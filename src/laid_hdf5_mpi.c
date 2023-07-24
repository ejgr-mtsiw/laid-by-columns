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
#include "types/word_t.h"
#include "utils/bit.h"
#include "utils/block.h"
#include "utils/clargs.h"
#include "utils/math.h"
#include "utils/mpi_custom.h"
#include "utils/output.h"
#include "utils/ranks.h"
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

	/**
	 * Parse command line arguments
	 */
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

	/**
	 * Create node-local communicator
	 * This communicator is used to share memory with processes intranode
	 */
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
	time_t main_tick = 0, main_tock = 0;

	if (rank == ROOT_RANK)
	{
		main_tick = time(0);
	}

	/**
	 * Local timing structures
	 */
	time_t tick = 0, tock = 0;

	ROOT_SHOWS("Using dataset '%s'\n", args.filename);
	ROOT_SHOWS("Using %d processes\n\n", size);

	/**
	 * Only rank 0 on a node actually reads the dataset and allocates memory
	 */
	uint64_t dset_data_size = 0;

	/**
	 * Open dataset file
	 */
	ROOT_SAYS("Initializing MPI RMA: ");
	TICK;

	if (node_rank == LOCAL_ROOT_RANK)
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

	/**
	 * All table pointers should now point to copy on noderank 0
	 */

	TOCK;

	/**
	 * Setup dataset
	 */
	if (node_rank == LOCAL_ROOT_RANK)
	{
		ROOT_SAYS("Reading dataset: ");
		TICK;

		/**
		 * Load dataset attributes
		 */
		hdf5_read_dataset_attributes(hdf5_dset.dataset_id, &dataset);

		/**
		 * Load dataset data
		 */
		hdf5_read_dataset_data(hdf5_dset.dataset_id, dataset.data);

		TOCK;
		/**
		 * Print dataset details
		 */
		ROOT_SHOWS("  Classes = %d", dataset.n_classes);
		ROOT_SHOWS(" [%d bits]\n", dataset.n_bits_for_class);
		ROOT_SHOWS("  Attributes = %d \n", dataset.n_attributes);
		ROOT_SHOWS("  Observations = %d \n", dataset.n_observations);

		/**
		 * We no longer need the dataset file
		 */
		hdf5_close_dataset(&hdf5_dset);

		TICK;

		/**Sort dataset
		 *
		 */
		ROOT_SAYS("Sorting dataset: ");

		/**
		 * We need to know the number of longs in each line of the dataset
		 * so we can't use the standard qsort implementation
		 */
		sort_r(dataset.data, dataset.n_observations,
			   dataset.n_words * sizeof(word_t), compare_lines_extra,
			   &dataset.n_words);

		TOCK;
		TICK;

		/**
		 * Remove duplicates
		 */
		ROOT_SAYS("Removing duplicates: ");

		unsigned int duplicates = remove_duplicates(&dataset);

		TOCK;
		ROOT_SHOWS("  %d duplicate(s) removed\n", duplicates);
		TICK;

		/**
		 * Fill class arrays
		 */
		ROOT_SAYS("Checking classes: ");

		// Number of classes
		uint32_t nc = dataset.n_classes;

		// Number of observations
		uint32_t no = dataset.n_observations;

		/**
		 * Array that stores the number of observations for each class
		 */
		uint32_t* n_opc = (uint32_t*) calloc(nc, sizeof(uint32_t));
		if (n_opc == NULL)
		{
			fprintf(stderr, "Error allocating n_observations_per_class\n");
			return NOK;
		}

		/**
		 * Matrix that stores the list of observations per class
		 */
		uint32_t* opc = (uint32_t*) calloc(nc * no, sizeof(uint32_t*));
		if (opc == NULL)
		{
			fprintf(stderr, "Error allocating observations_per_class\n");
			return NOK;
		}

		dataset.n_observations_per_class = n_opc;
		dataset.observations_per_class	 = opc;

		if (fill_class_arrays(&dataset) != OK)
		{
			return EXIT_FAILURE;
		}

		TOCK;

		for (unsigned int i = 0; i < dataset.n_classes; i++)
		{
			ROOT_SHOWS("  Class %d: ", i);
			ROOT_SHOWS("%d item(s)\n", dataset.n_observations_per_class[i]);
		}

		TICK;

		/**
		 * Set JNSQ
		 */
		ROOT_SAYS("Setting up JNSQ attributes: ");

		unsigned int max_jnsq = add_jnsqs(&dataset);

		TOCK;
		ROOT_SHOWS("  Max JNSQ: %d", max_jnsq);
		ROOT_SHOWS(" [%d bits]\n", dataset.n_bits_for_jnsqs);

		/**
		 * Displayes the estimated disjoint matrix size
		 * This is the size it would take if it was saved on disk
		 */
		dm.n_matrix_lines = get_dm_n_lines(&dataset);

		double matrix_size
			= ((double) dm.n_matrix_lines
			   * (dataset.n_attributes + dataset.n_bits_for_class))
			/ (1024 * 1024 * 8);

		ROOT_SHOWS("  %d matrix steps generated\n", dm.n_matrix_lines);
		ROOT_SHOWS("  Estimated disjoint matrix size: %3.2fMB\n", matrix_size);

		TICK;
		ROOT_SAYS("Sorting dataset by class: ");
		// Sort by class
		sort_dataset_by_class(&dataset);

		free(dataset.observations_per_class);
		dataset.observations_per_class = NULL;

		TOCK;
	}

	ROOT_SAYS("Broadcasting attributes: ");
	TICK;

	uint32_t toshare[5];
	if (node_rank == 0)
	{
		toshare[0] = dataset.n_classes;
		toshare[1] = dataset.n_attributes;
		toshare[2] = dataset.n_observations;
		toshare[3] = dataset.n_words;
		toshare[4] = dm.n_matrix_lines;

		MPI_Bcast(&toshare, 5, MPI_UINT32_T, 0, node_comm);
	}
	else
	{
		MPI_Bcast(&toshare, 5, MPI_UINT32_T, 0, node_comm);

		dataset.n_classes	   = toshare[0];
		dataset.n_attributes   = toshare[1];
		dataset.n_observations = toshare[2];
		dataset.n_words		   = toshare[3];
		dm.n_matrix_lines	   = toshare[4];
	}

	TOCK;

	ROOT_SAYS("Broadcasting observations per class: ");
	TICK;

	if (node_rank == 0)
	{
		MPI_Bcast(dataset.n_observations_per_class, dataset.n_classes,
				  MPI_UINT32_T, 0, node_comm);
	}
	else
	{
		dataset.n_observations_per_class
			= (uint32_t*) malloc(dataset.n_classes * sizeof(uint32_t));
		MPI_Bcast(dataset.n_observations_per_class, dataset.n_classes,
				  MPI_UINT32_T, 0, node_comm);
	}

	TOCK;

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

	ROOT_SAYS("Applying set covering algorithm:\n");
	TICK;

	// Distribute attributes by the processes
	dm.a_offset = BLOCK_LOW(rank, size, dataset.n_words);
	dm.a_size	= BLOCK_SIZE(rank, size, dataset.n_words);

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
	 * the attribute totals calculation, avoiding one if in the inner cycle
	 */
	uint32_t* attribute_totals
		= (uint32_t*) calloc(dm.a_size * WORD_BITS, sizeof(uint32_t));

	/**
	 * Selected attributes aka the solution, only root needs them
	 */
	word_t* selected_attributes = NULL;
	if (rank == 0)
	{
		selected_attributes = (word_t*) calloc(dataset.n_words, sizeof(word_t));
	}

	/***
	 * Calculate class offsets
	 */
	uint32_t* class_offsets
		= (uint32_t*) malloc(dataset.n_classes * sizeof(uint32_t));
	uint32_t running_total = 0;
	for (uint32_t i = 0; i < dataset.n_classes; i++)
	{
		class_offsets[i] = running_total;
		running_total += dataset.n_observations_per_class[i];
	}

	//*********************************************************/
	// BUILD INITIAL TOTALS
	//*********************************************************/
	for (uint32_t c = 0; c < dataset.n_classes - 1; c++)
	{
		word_t* start_class_a
			= dataset.data + class_offsets[c] * dataset.n_words + dm.a_offset;
		word_t* end_class_a = start_class_a
			+ dataset.n_observations_per_class[c] * dataset.n_words;

		word_t* start_class_b = end_class_a;
		word_t* end_class_b	  = dataset.data
			+ (dataset.n_observations) * dataset.n_words + dm.a_offset;

		//	printf("\n");
		//	printf("[%d] start_a: %ld, end_a: %ld, start_b= %ld, end_b= %ld",
		// rank, 		   start_class_a - dataset.data, end_class_a -
		// dataset.data, 		   start_class_b - dataset.data, end_class_b -
		// dataset.data); 	printf("[%d] c: %d, offset: %d, size= %d", rank, c,
		// dm.a_offset, 		   dm.a_size); 	printf("\n");

		for (word_t* la = start_class_a; la < end_class_a;
			 la += dataset.n_words)
		{
			for (word_t* lb = start_class_b; lb < end_class_b;
				 lb += dataset.n_words)
			{
				uint32_t c_attribute = 0;

				//	printf("[%d] %ld x %ld\n", rank, la - dataset.data, lb -
				// dataset.data);

				for (uint32_t w = 0; w < dm.a_size; w++)
				{
					word_t lxor = la[w] ^ lb[w];

					for (int8_t bit = WORD_BITS - 1; bit >= 0;
						 bit--, c_attribute++)
					{
						attribute_totals[c_attribute] += BIT_CHECK(lxor, bit);
					}
				}
			}
		}
	}

	//	printf("FINAL ATTR TOTALS:\n");
	//	for (unsigned long i = 0; i < dataset.n_attributes; i++)
	//	{
	//		printf("i: %ld, %d\n", i, attribute_totals[i]);
	//	}

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

	/**
	 * Number of covered lines so far
	 */
	uint32_t n_covered_lines = 0;

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
		best_attribute_t best_max = { .n_covered_lines = 0, .attribute = -1 };

		MPI_Allreduce(&local_max, &best_max, 1, ctype, myOp, comm);
		/**
		 * At this point, the answer resides on best_max
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
			ROOT_SHOWS("  Selected attribute #%ld, ", best_max.attribute);
			ROOT_SHOWS("covers %lu lines ", best_max.n_covered_lines);
			TOCK;
			TICK;

			// Which word has the best attribute
			uint32_t best_word = best_max.attribute / WORD_BITS;

			// Which bit?
			uint32_t best_bit = WORD_BITS - best_max.attribute % WORD_BITS - 1;

			// Mark best attribute as selected
			BIT_SET(selected_attributes[best_word], best_bit);
		}

		// Update number of covered lines
		n_covered_lines += best_max.n_covered_lines;

		// If we have everything covered goto show_solution
		if (n_covered_lines == dm.n_matrix_lines)
		{
			goto show_solution;
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
			uint32_t line = 0;
			for (uint32_t c = 0; c < dataset.n_classes - 1; c++)
			{
				word_t* start_class_a = dataset.data
					+ class_offsets[c] * dataset.n_words + best_word;
				word_t* end_class_a = start_class_a
					+ dataset.n_observations_per_class[c] * dataset.n_words;

				word_t* start_class_b = end_class_a;
				word_t* end_class_b	  = dataset.data
					+ (dataset.n_observations) * dataset.n_words + best_word;

				for (word_t* la = start_class_a; la < end_class_a;
					 la += dataset.n_words)
				{
					for (word_t* lb = start_class_b; lb < end_class_b;
						 lb += dataset.n_words)
					{
						word_t lxor = *la ^ *lb;

						if (BIT_CHECK(lxor, best_bit))
						{
							// Where to save it
							uint32_t current_word = line / WORD_BITS;

							// Which bit?
							uint32_t current_bit
								= WORD_BITS - line % WORD_BITS - 1;

							BIT_SET(best_column[current_word], current_bit);
						}

						line++;
					}
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
		uint32_t line = 0;
		for (uint32_t c = 0; c < dataset.n_classes - 1; c++)
		{
			word_t* start_class_a = dataset.data
				+ class_offsets[c] * dataset.n_words + dm.a_offset;
			word_t* end_class_a = start_class_a
				+ dataset.n_observations_per_class[c] * dataset.n_words;

			word_t* start_class_b = end_class_a;
			word_t* end_class_b	  = dataset.data
				+ (dataset.n_observations) * dataset.n_words + dm.a_offset;

			for (word_t* la = start_class_a; la < end_class_a;
				 la += dataset.n_words)
			{
				for (word_t* lb = start_class_b; lb < end_class_b;
					 lb += dataset.n_words, line++)
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

					uint32_t c_attribute = 0;

					//	printf("[%d] %ld x %ld\n", rank, la - dataset.data, lb -
					// dataset.data);

					for (uint32_t w = 0; w < dm.a_size; w++)
					{
						word_t lxor = la[w] ^ lb[w];

						for (int8_t bit = WORD_BITS - 1; bit >= 0;
							 bit--, c_attribute++)
						{
							attribute_totals[c_attribute]
								-= BIT_CHECK(lxor, bit);
						}
					}
				}
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
	if (rank == ROOT_RANK)
	{
		fprintf(stdout, "Solution: { ");
		uint32_t current_attribute = 0;
		for (uint32_t w = 0; w < dataset.n_words; w++)
		{
			for (int8_t bit = WORD_BITS - 1; bit >= 0;
				 bit--, current_attribute++)
			{
				if (BIT_CHECK(selected_attributes[w], bit))
				{
					fprintf(stdout, "%d ", current_attribute);
				}
			}
		}
		fprintf(stdout, "}\nAll done! ");

		main_tock = time(0);
		fprintf(stdout, "[%lds]\n", main_tock - main_tick);
	}

	//  wait for everyone
	MPI_Barrier(comm);

	MPI_Win_free(&win_shared_dset);

	dataset.data = NULL;
	free_dataset(&dataset);

	free(attribute_totals);
	free(best_column);
	free(covered_lines);
	free(selected_attributes);

	/* shut down MPI */
	MPI_Finalize();

	return EXIT_SUCCESS;
}
