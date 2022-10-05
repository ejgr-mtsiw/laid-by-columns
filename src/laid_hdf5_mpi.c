/*
 ============================================================================
 Name        : laid_hdf5_mpi.c
 Author      : Eduardo Ribeiro
 Description : OpenMPI implementation of the LAID algorithm in C + HDF5
 ============================================================================
 */

#include "dataset.h"
#include "dataset_hdf5.h"
#include "disjoint_matrix.h"
#include "jnsq.h"
#include "set_cover.h"
#include "types/best_attribute_t.h"
#include "types/dataset_hdf5_t.h"
#include "types/dataset_t.h"
#include "types/dm_t.h"
#include "types/steps_t.h"
#include "types/word_t.h"
#include "utils/bit.h"
#include "utils/block.h"
#include "utils/clargs.h"
#include "utils/math.h"
#include "utils/mpi_custom.h"
#include "utils/sort_r.h"
#include "utils/timing.h"

#include "hdf5.h"
#include "mpi.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

/**
 * In this mode we don't write the disjoint matrix (DM).
 * Everytime we need a line or column from the DM it's generated from the
 * dataset in memory.
 *
 * Each node root will open the original dataset file and read the contents to
 * memory. This allows us to save memory in each node, without sacrificing much
 * performance, because we're not having to send data across nodes.
 *
 * The node root(s) sort the dataset in memory, remove duplicates and  adds
 * jnsqs bits if necessary.
 *
 * Then they build a list of the steps needed to generate the disjoint matrix.
 * This list of steps allows us to generate any line or column of the disjoint
 * matrix.
 *
 * TLDR:
 * Each node root
 *  - Reads dataset attributes from hdf5 file
 *  - Read dataset
 *  - Sort dataset
 *  - Remove duplicates
 *  - Add jnsqs
 *  - Builds steps for matrix generation
 *
 * All processes
 *  - Apply set covering algorithm
 *
 * Global root
 *  - Show solution
 */
int main(int argc, char** argv)
{
	/**
	 * Command line arguments set by the user
	 */
	clargs_t args;

	// Parse command line arguments
	if (read_args(argc, argv, &args) == READ_CL_ARGS_NOK)
	{
		return EXIT_FAILURE;
	}

	/*
	 * Initialize MPI
	 */
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
	{
		printf("Error initializing MPI environment!\n");
		return EXIT_FAILURE;
	}

	/**
	 * Rank of process
	 */
	int rank;

	/**
	 * Number of processes
	 */
	int size;

	/**
	 * Global communicator group
	 */
	MPI_Comm comm = MPI_COMM_WORLD;

	/**
	 * Setup global rank and size
	 */
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	/**
	 * Node communicator group
	 */
	MPI_Comm node_comm = MPI_COMM_NULL;

	// Create node-local communicator
	// This communicator is used to share memory with processes intranode
	MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL,
						&node_comm);

	/**
	 * In-node rank of process
	 */
	int node_rank;

	/**
	 * In-node number of processes
	 */
	int node_size;

	/**
	 * Setup node rank and size
	 */
	MPI_Comm_size(node_comm, &node_size);
	MPI_Comm_rank(node_comm, &node_rank);

	/**
	 * The dataset
	 */
	dataset_t dataset;
	init_dataset(&dataset);

	/**
	 * The HDF5 dataset file
	 */
	dataset_hdf5_t hdf5_dset;

	/**
	 * The disjoint matrix info
	 */
	dm_t dm;

	/**
	 * Timing for the full operation
	 */
	struct timespec main_tick, main_tock;
	if (rank == 0)
	{
		clock_gettime(CLOCK_MONOTONIC_RAW, &main_tick);
	}

	/**
	 * Local timing structures
	 */
	SETUP_TIMING;
	TICK;

	// Only rank 0 on a node actually reads the dataset and allocates memory
	uint64_t dset_data_size = 0;

	// Open dataset file
	if (node_rank == 0)
	{
		if (hdf5_open_dataset(args.filename, args.datasetname, &hdf5_dset)
			== NOK)
		{
			return EXIT_FAILURE;
		}

		dataset.n_observations = hdf5_dset.dimensions[0];
		dataset.n_words		   = hdf5_dset.dimensions[1];

		dset_data_size = dataset.n_observations * dataset.n_words;
	}

	word_t* dset_data		= NULL;
	MPI_Win win_shared_dset = MPI_WIN_NULL;
	MPI_Win_allocate_shared(dset_data_size * sizeof(word_t), sizeof(word_t),
							MPI_INFO_NULL, node_comm, &dset_data,
							&win_shared_dset);

	// Set dataset data pointer
	if (node_rank == 0)
	{
		dataset.data = dset_data;
	}
	else
	{
		MPI_Aint win_size;
		int win_disp;
		MPI_Win_shared_query(win_shared_dset, 0, &win_size, &win_disp,
							 &dataset.data);
	}

	if (rank == 0)
	{
		fprintf(stdout, "- Finished MPI RMA Init ");
		TOCK(stdout);
	}
	// All table pointers should now point to copy on noderank 0

	// Setup dataset
	if (node_rank == 0)
	{
		TICK;

		fprintf(stdout, "- Loading dataset data\n - ");

		// Load dataset attributes
		hdf5_read_dataset_attributes(hdf5_dset.dataset_id, &dataset);

		// Load dataset data
		hdf5_read_dataset_data(hdf5_dset.dataset_id, dataset.data);

		print_dataset_details(stdout, &dataset);

		fprintf(stdout, " - Finished loading dataset data ");

		// We no longer need the dataset file
		hdf5_close_dataset(&hdf5_dset);

		TOCK(stdout);
		TICK;

		// Sort dataset
		fprintf(stdout, "- Sorting dataset\n");

		// We need to know the number of longs in each line of the dataset
		// so we can't use the standard qsort implementation
		sort_r(dataset.data, dataset.n_observations,
			   dataset.n_words * sizeof(word_t), compare_lines_extra,
			   &dataset.n_words);

		fprintf(stdout, " - Sorted dataset");
		TOCK(stdout);
		TICK;

		// Remove duplicates
		fprintf(stdout, "- Removing duplicates:\n");

		unsigned int duplicates = remove_duplicates(&dataset);

		fprintf(stdout, " - %d duplicate(s) removed ", duplicates);
		TOCK(stdout);
		TICK;

		// Fill class arrays
		fprintf(stdout, "- Checking classes: ");

		if (fill_class_arrays(&dataset) != OK)
		{
			return EXIT_FAILURE;
		}

		TOCK(stdout);

		for (unsigned int i = 0; i < dataset.n_classes; i++)
		{
			fprintf(stdout, " - class %d: %d item(s)\n", i,
					dataset.n_observations_per_class[i]);
		}

		TICK;

		// Set JNSQ
		fprintf(stdout, "- Setting up JNSQ attributes:\n");

		unsigned int max_jnsq = add_jnsqs(&dataset);

		fprintf(stdout, " - Max JNSQ: %d [%d bits] ", max_jnsq,
				dataset.n_bits_for_jnsqs);
		TOCK(stdout);
	}

	// End setup dataset
	// MPI_Barrier(node_comm);

	dm.n_matrix_lines = 0;
	if (node_rank == 0)
	{
		dm.n_matrix_lines = get_dm_n_lines(&dataset);
	}

	steps_t* localsteps		 = NULL;
	MPI_Win win_shared_steps = MPI_WIN_NULL;
	MPI_Win_allocate_shared(dm.n_matrix_lines * sizeof(steps_t),
							sizeof(steps_t), MPI_INFO_NULL, node_comm,
							&localsteps, &win_shared_steps);

	/**
	 * The steps
	 */
	steps_t* steps = NULL;

	// Set dataset data pointer
	if (node_rank == 0)
	{
		steps = localsteps;
	}
	else
	{
		MPI_Aint win_size;
		int win_disp;
		MPI_Win_shared_query(win_shared_steps, 0, &win_size, &win_disp, &steps);
	}
	// All table pointers should now point to copy on noderank 0

	// Setup steps
	if (node_rank == 0)
	{
		TICK;

		fprintf(stdout, "- Generating matrix steps\n");

		uint32_t nc	   = dataset.n_classes;
		uint32_t no	   = dataset.n_observations;
		uint32_t* opc  = dataset.observations_per_class;
		uint32_t* nopc = dataset.n_observations_per_class;

		// DO IT
		uint32_t cs = 0;

		for (uint32_t ca = 0; ca < nc - 1; ca++)
		{
			for (uint32_t ia = 0; ia < nopc[ca]; ia++)
			{
				for (uint32_t cb = ca + 1; cb < nc; cb++)
				{
					for (uint32_t ib = 0; ib < nopc[cb]; ib++)
					{
						steps[cs].indexA = opc[ca * no + ia];
						steps[cs].indexB = opc[cb * no + ib];

						cs++;
					}
				}
			}
		}

		free(dataset.n_observations_per_class);
		free(dataset.observations_per_class);

		dataset.n_observations_per_class = NULL;
		dataset.observations_per_class	 = NULL;

		fprintf(stdout, " - Finished generating %d matrix steps ",
				dm.n_matrix_lines);
		TOCK(stdout);
	}

	if (rank == 0)
	{
		fprintf(stdout, "- Broadcasting attributes\n");
	}

	uint64_t toshare[4];
	if (node_rank == 0)
	{
		toshare[0] = dataset.n_attributes;
		toshare[1] = dataset.n_observations;
		toshare[2] = dataset.n_words;
		toshare[3] = dm.n_matrix_lines;

		MPI_Bcast(&toshare, 4, MPI_UINT64_T, 0, node_comm);
	}
	else
	{
		MPI_Bcast(&toshare, 4, MPI_UINT64_T, 0, node_comm);

		dataset.n_attributes   = toshare[0];
		dataset.n_observations = toshare[1];
		dataset.n_words		   = toshare[2];
		dm.n_matrix_lines	   = toshare[3];
	}

	dm.steps = steps;

	// Distribute attributes by the processes
	// Keep them in multiples of N_WORDS_CACHE_LINE words
	// to maximize cache line usage
	dm.a_offset
		= BLOCK_LOW_MULTIPLE(rank, size, dataset.n_words, N_WORDS_CACHE_LINE);
	dm.a_size
		= BLOCK_SIZE_MULTIPLE(rank, size, dataset.n_words, N_WORDS_CACHE_LINE);

	if (rank == 0)
	{
		fprintf(stdout, " - Finished broadcasting attributes\n");
		TICK;
	}

	/**
	 * ALL:
	 *  - Setup line covered array -> 0
	 * ROOT:
	 *  - Setup attribute covered array -> 0
	 *
	 * ALL:
	 *  - Calculate attributes totals
	 *
	 *LOOP:
	 * All:
	 *  - Calculate best attribute
	 *  - MPI_Allreduce to select best global attribute
	 *
	 * ALL:
	 *  - if there are no more lines to blacklist (attribute == -1):
	 *  - goto SHOW_SOLUTION
	 *
	 * ROOT:
	 *  - Marks best attribute as selected
	 *
	 * NODE ROOTS:
	 *  - Generates cover line for the best atribute
	 *  - MPI_Bcast to everyone else on the same node
	 *
	 * PROCESSES WITHOUT THE BEST ATTRIBUTE:
	 *  - Wait for cover line
	 *
	 * ALL:
	 *  - Update atributes totals based on the received cover line
	 *  - Update line covered array
	 *
	 * ALL:
	 *  - goto LOOP
	 *
	 *SHOW_SOLUTION:
	 * ROOT:
	 *  - Display list of selected attributes
	 */

	if (rank == 0)
	{
		TICK;

		printf("- Applying set covering algorithm\n");
	}

	// The covered lines and covered attributes are bit arrays

	/**
	 * Number of words needed to store a column (attribute data)
	 */
	uint32_t n_words_in_column
		= dm.n_matrix_lines / WORD_BITS + (dm.n_matrix_lines % WORD_BITS != 0);

	/**
	 * Bit array with the information about the lines covered (1) or not (0) by
	 * the current best attribute
	 */
	word_t* best_column = (word_t*) calloc(n_words_in_column, sizeof(word_t));

	/**
	 * Bit array with the information about the lines already covered (1) or
	 * not (0) so far
	 */
	word_t* covered_lines = (word_t*) calloc(n_words_in_column, sizeof(word_t));

	/**
	 * The number of attributes is rounded so we can check all bits during
	 * the attribute totals calculation
	 */
	uint32_t* attribute_totals
		= (uint32_t*) calloc(dm.a_size * WORD_BITS, sizeof(uint32_t));

	/**
	 * Selected attributes aka the solution
	 */
	word_t* selected_attributes = NULL;

	/**
	 * Only root needs them
	 */
	if (rank == 0)
	{
		selected_attributes = (word_t*) calloc(dataset.n_words, sizeof(word_t));
	}

	//*********************************************************/
	// BUILD INITIAL TOTALS
	//*********************************************************/
	// calculate_initial_totals(&dm, &dataset, attribute_totals);
	for (uint32_t line = 0; line < dm.n_matrix_lines; line++)
	{
		word_t* la = dataset.data + dm.steps[line].indexA * dataset.n_words
			+ dm.a_offset;
		word_t* lb = dataset.data + dm.steps[line].indexB * dataset.n_words
			+ dm.a_offset;

		uint32_t c_attribute = 0;

		for (uint32_t w = 0; w < dm.a_size; w++)
		{
			word_t lxor = la[w] ^ lb[w];

			for (int8_t bit = WORD_BITS - 1; bit >= 0; bit--, c_attribute++)
			{
				attribute_totals[c_attribute] += BIT_CHECK(lxor, bit);
			}
		}
	}
	//*********************************************************/
	// END BUILD INITIAL TOTALS
	//*********************************************************/

	/**
	 * Build function to compare best atributes in MPI_Allreduce
	 * In C we don't have a MPI_LONG_LONG custom type so we cretate a new
	 * type and function to select the best attribute of all the best
	 * attributes selected by each process
	 */
	MPI_Op myOp;
	MPI_Datatype ctype;

	// explain to MPI how type best_attribute is defined
	MPI_Type_contiguous(2, MPI_LONG, &ctype);
	MPI_Type_commit(&ctype);

	// create the user-op
	MPI_Op_create(MPI_get_best_attribute, true, &myOp);

	while (true)
	{
		// What is my best attribute
		best_attribute_t local_max = get_best_attribute(
			attribute_totals,
			MIN(dm.a_size * WORD_BITS,
				dataset.n_attributes - (dm.a_offset * WORD_BITS)));

		// My best attribute translates to which global attribute?
		if (local_max.attribute != -1)
		{
			local_max.attribute += dm.a_offset * WORD_BITS;
		}

		//	for (int r = 0; r < size; r++)
		//	{
		//		if (r == rank)
		//		{
		//			printf("[%d] max:%ld, attr: %ld\n", rank, local_max.max,
		//				   local_max.attribute);
		//		}
		//		sleep(1);
		//	}

		// Reset best
		best_attribute_t best_max = { .max = 0, .attribute = -1 };

		// At this point, the answer resides on best_max
		MPI_Allreduce(&local_max, &best_max, 1, ctype, myOp, comm);

		/**
		 * If the best attribute is -1 we're done here
		 */
		if (best_max.attribute == -1)
		{
			goto show_solution;
		}

		/**
		 * Mark attributed as selected
		 */
		if (rank == 0)
		{
			printf(" - Selected attribute #%ld, covers %ld lines ",
				   best_max.attribute, best_max.max);

			TOCK(stdout);
			TICK;

			// Which word has the best attribute
			uint32_t best_word = best_max.attribute / WORD_BITS;

			// Which bit?
			uint32_t best_bit = WORD_BITS - best_max.attribute % WORD_BITS - 1;

			// Mark best attribute as selected
			BIT_SET(selected_attributes[best_word], best_bit);
		}

		/**
		 * As all processes have access to the full dataset each node root will
		 * generate the cover line for it's node
		 */
		if (node_rank == 0)
		{
			// Which word has the best attribute
			uint32_t best_word = best_max.attribute / WORD_BITS;

			// Which bit?
			uint32_t best_bit = WORD_BITS - best_max.attribute % WORD_BITS - 1;

			//***********************************************************/
			// BUILD BEST COLUMN
			//***********************************************************/
			for (uint32_t line = 0; line < dm.n_matrix_lines; line++)
			{
				word_t* la = dataset.data
					+ dm.steps[line].indexA * dataset.n_words + best_word;
				word_t* lb = dataset.data
					+ dm.steps[line].indexB * dataset.n_words + best_word;

				word_t lxor = *la ^ *lb;

				if (BIT_CHECK(lxor, best_bit))
				{
					// Where to save it
					uint32_t current_word = line / WORD_BITS;

					// Which bit?
					uint32_t current_bit = WORD_BITS - line % WORD_BITS - 1;

					BIT_SET(best_column[current_word], current_bit);
				}
			}
			//***********************************************************/
			// END BUILD BEST COLUMN
			//***********************************************************/
		}

		MPI_Bcast(best_column, n_words_in_column, MPI_UINT64_T, 0, node_comm);

		//***************************************************************/
		// UPDATE ATTRIBUTES TOTALS
		//***************************************************************/

		/**
		 * Subtract the contribution of the lines covered by the best attribute
		 * from the totals
		 */
		for (uint32_t line = 0; line < dm.n_matrix_lines; line++)
		{
			// Is this line already covered?
			// Yes: skip
			// No. Is it covered by the best attribute?
			// Yes: add
			// No: skip
			uint32_t current_word = line / WORD_BITS;
			uint8_t current_bit	  = WORD_BITS - line % WORD_BITS - 1;

			// Is this line already covered?
			if (BIT_CHECK(covered_lines[current_word], current_bit))
			{
				// This line is already covered: skip
				continue;
			}

			// Is this line covered by the best attribute?
			if (!BIT_CHECK(best_column[current_word], current_bit))
			{
				// This line is NOT covered: skip
				continue;
			}

			// This line was uncovered, but is covered now
			// Calculate attributes totals
			word_t* la = dataset.data + dm.steps[line].indexA * dataset.n_words
				+ dm.a_offset;
			word_t* lb = dataset.data + dm.steps[line].indexB * dataset.n_words
				+ dm.a_offset;

			uint32_t c_attribute = 0;

			for (uint32_t w = 0; w < dm.a_size; w++)
			{
				word_t lxor = la[w] ^ lb[w];

				for (int8_t bit = WORD_BITS - 1; bit >= 0; bit--, c_attribute++)
				{
					attribute_totals[c_attribute] -= BIT_CHECK(lxor, bit);
				}

				//				__m256i lxor_ = _mm256_set1_epi64x(lxor);
				//				__m256i zeros = _mm256_set1_epi64x(0);
				//
				//				for (int8_t bit = WORD_BITS, pos = 0; bit > 0;
				//					 bit -= 4, pos += 4, c_attribute += 4)
				//				{
				//					__m256i mask = _mm256_set_epi64x(
				//						AND_MASK_TABLE[bit - 1],
				// AND_MASK_TABLE[bit
				//- 2], 						AND_MASK_TABLE[bit - 3],
				// AND_MASK_TABLE[bit - 4]);
				//					__m256i attr_totals = _mm256_set_epi64x(
				//						attribute_totals[pos],
				// attribute_totals[pos
				//+ 1], 						attribute_totals[pos + 2],
				// attribute_totals[pos + 3]);
				//
				//					__m256i a	 = _mm256_and_si256(lxor_,
				// mask);
				//					__m256i ones = _mm256_cmpgt_epi64(a, zeros);
				//
				//					attr_totals = _mm256_sub_epi64(attr_totals,
				// ones);
				//
				//					attribute_totals[c_attribute + 0] =
				// attr_totals[0]; 					attribute_totals[c_attribute
				// + 1]
				// =
				// attr_totals[1]; 					attribute_totals[c_attribute
				// + 2]
				// =
				// attr_totals[2]; 					attribute_totals[c_attribute
				// + 3] = attr_totals[3];
				//				}
			}
		}
		//***************************************************************/
		// END UPDATE ATTRIBUTES TOTALS
		//***************************************************************/

		//***************************************************************/
		// UPDATE COVERED LINES
		//***************************************************************/
		for (uint32_t w = 0; w < n_words_in_column; w++)
		{
			covered_lines[w] |= best_column[w];
		}

		//***************************************************************/
		// END UPDATE COVERED LINES
		//***************************************************************/
	}

show_solution:
	if (rank == 0)
	{
		printf("Solution: { ");
		uint32_t current_attribute = 0;
		for (uint32_t w = 0; w < dataset.n_words; w++)
		{
			for (int8_t bit = WORD_BITS - 1; bit >= 0;
				 bit--, current_attribute++)
			{
				if (BIT_CHECK(selected_attributes[w], bit))
				{
					printf("%d ", current_attribute);
				}
			}
		}
		printf("}\n");
	}

	if (rank == 0)
	{
		fprintf(stdout, "All done! ");

		clock_gettime(CLOCK_MONOTONIC_RAW, &main_tock);
		fprintf(stdout, "[%0.3fs]\n",
				(main_tock.tv_nsec - main_tick.tv_nsec) / 1000000000.0
					+ (main_tock.tv_sec - main_tick.tv_sec));
	}

	//  wait for everyone
	MPI_Barrier(comm);

	MPI_Win_free(&win_shared_steps);
	MPI_Win_free(&win_shared_dset);

	dataset.data = NULL;
	dm.steps	 = NULL;

	free_dataset(&dataset);
	free_dm(&dm);

	free(attribute_totals);
	free(best_column);
	free(covered_lines);
	free(selected_attributes);

	/* shut down MPI */
	MPI_Finalize();

	return EXIT_SUCCESS;
}
