/*
 ============================================================================
 Name        : laid_by_columns.c
 Author      : Eduardo Ribeiro
 Description : OpenMPI implementation of the LAID algorithm in C + HDF5
 ============================================================================
 */

#include "dataset.h"
#include "dataset_hdf5.h"
#include "disjoint_matrix.h"
#include "disjoint_matrix_mpi.h"
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

#include <assert.h>
#include <math.h>
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
	 * Timing for the full operation
	 */
	SETUP_TIMING_GLOBAL;

	/**
	 * Local timing structures
	 */
	SETUP_TIMING;

	/**
	 * The dataset
	 */
	dataset_t dataset;
	init_dataset(&dataset);

	/**
	 * The HDF5 dataset file
	 */
	dataset_hdf5_t hdf5_dset;

	// Open dataset file
	ROOT_SHOWS("Using dataset '%s'\n", args.filename);
	ROOT_SHOWS("Using %d processes\n\n", size);

	/**
	 * Dataset data size
	 * Only rank 0 on a node actually reads the dataset and allocates memory
	 */
	uint64_t shared_data_size = 0;

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

		shared_data_size = dataset.n_observations * dataset.n_words;
	}

	word_t* dset_data		= NULL;
	MPI_Win win_shared_dset = MPI_WIN_NULL;
	MPI_Win_allocate_shared(shared_data_size * sizeof(word_t), sizeof(word_t),
							MPI_INFO_NULL, node_comm, &dset_data,
							&win_shared_dset);

	// Set dataset data pointer
	if (node_rank == LOCAL_ROOT_RANK)
	{
		dataset.data = dset_data;
	}
	else
	{
		MPI_Aint win_size;
		int win_disp;
		MPI_Win_shared_query(win_shared_dset, LOCAL_ROOT_RANK, &win_size,
							 &win_disp, &dataset.data);
	}
	// All dataset.data pointers should now point to copy on noderank 0

	TOCK;

	// Setup dataset

	ROOT_SAYS("Reading dataset: ");
	TICK;

	if (node_rank == LOCAL_ROOT_RANK)
	{
		// Load dataset attributes
		hdf5_read_dataset_attributes(hdf5_dset.dataset_id, &dataset);

		// Load dataset data
		hdf5_read_dataset_data(hdf5_dset.dataset_id, dataset.data);

		TOCK;

		// Print dataset details
		ROOT_SHOWS("  Classes = %d", dataset.n_classes);
		ROOT_SHOWS(" [%d bits]\n", dataset.n_bits_for_class);
		ROOT_SHOWS("  Attributes = %d \n", dataset.n_attributes);
		ROOT_SHOWS("  Observations = %d \n", dataset.n_observations);

		// We no longer need the dataset file
		hdf5_close_dataset(&hdf5_dset);

		// Sort dataset
		ROOT_SAYS("Sorting dataset: ");
		TICK;

		/**
		 * We need to know the number of longs in each line of the dataset
		 * so we can't use the standard qsort implementation
		 */
		sort_r(dataset.data, dataset.n_observations,
			   dataset.n_words * sizeof(word_t), compare_lines_extra,
			   &dataset.n_words);

		TOCK;

		// Remove duplicates
		ROOT_SAYS("Removing duplicates: ");
		TICK;

		unsigned int duplicates = remove_duplicates(&dataset);

		TOCK;
		ROOT_SHOWS("  %d duplicate(s) removed\n", duplicates);
	}

	// Share current dataset attributes
	MPI_Bcast(&(dataset.n_attributes), 1, MPI_UINT32_T, LOCAL_ROOT_RANK,
			  node_comm);
	MPI_Bcast(&(dataset.n_observations), 1, MPI_UINT32_T, LOCAL_ROOT_RANK,
			  node_comm);
	MPI_Bcast(&(dataset.n_classes), 1, MPI_UINT32_T, LOCAL_ROOT_RANK,
			  node_comm);
	MPI_Bcast(&(dataset.n_bits_for_class), 1, MPI_UINT8_T, LOCAL_ROOT_RANK,
			  node_comm);
	MPI_Bcast(&(dataset.n_words), 1, MPI_UINT32_T, LOCAL_ROOT_RANK, node_comm);

	// Fill class arrays
	ROOT_SAYS("Checking classes: ");
	TICK;

	/**
	 * Array that stores the number of observations for each class
	 */
	dataset.n_observations_per_class
		= (uint32_t*) calloc(dataset.n_classes, sizeof(uint32_t));
	assert(dataset.n_observations_per_class != NULL);

	/**
	 * Matrix that stores the list of observations per class
	 */
	dataset.observations_per_class = (word_t**) calloc(
		dataset.n_classes * dataset.n_observations, sizeof(word_t*));
	assert(dataset.observations_per_class != NULL);

	// Must make sure the dataset is filled before proceeding
	MPI_Barrier(node_comm);

	if (fill_class_arrays(&dataset) != OK)
	{
		return EXIT_FAILURE;
	}

	TOCK;

	if (rank == ROOT_RANK)
	{
		for (unsigned int i = 0; i < dataset.n_classes; i++)
		{
			ROOT_SHOWS("  Class %d: ", i);
			ROOT_SHOWS("%d item(s)\n", dataset.n_observations_per_class[i]);
		}
	}

	// Must make sure everyone has finished before changing the dataset
	MPI_Barrier(node_comm);

	// Set JNSQ
	if (node_rank == LOCAL_ROOT_RANK)
	{
		ROOT_SAYS("Setting up JNSQ attributes: ");
		TICK;

		uint32_t max_inconsistency = add_jnsqs(&dataset);

		// Update number of bits needed for jnsqs
		if (max_inconsistency > 0)
		{
			// How many bits are needed for jnsq attributes
			uint8_t n_bits_for_jnsq = ceil(log2(max_inconsistency + 1));

			dataset.n_bits_for_jnsqs = n_bits_for_jnsq;
		}

		TOCK;
		ROOT_SHOWS("  Max JNSQ: %d", max_inconsistency);
		ROOT_SHOWS(" [%d bits]\n", dataset.n_bits_for_jnsqs);
	}

	// Update dataset data because only node_root knows if we added jnsqs
	MPI_Bcast(&(dataset.n_bits_for_jnsqs), 1, MPI_UINT8_T, LOCAL_ROOT_RANK,
			  node_comm);

	// JNSQ attributes are treated just like all the other attributes from
	// this point forward
	dataset.n_attributes += dataset.n_bits_for_jnsqs;

	// n_words may have changed?
	// If we have 5 classes (3 bits) and only one bit is in the last word
	// If we only use 2 bits for jnsqs then we need one less n_words
	// This could impact all the calculations that assume n_words - 1
	// is the last word!
	dataset.n_words = dataset.n_attributes / WORD_BITS
		+ (dataset.n_attributes % WORD_BITS != 0);

	// End setup dataset

	MPI_Barrier(node_comm);

	// Setup disjoint matrix

	/**
	 * The disjoint matrix info
	 */
	dm_t dm;

	ROOT_SAYS("Calculating disjoint matrix lines to generate: ");
	TICK;

	// Calculate the number of disjoint matrix lines
	dm.n_matrix_lines = get_dm_n_lines(&dataset);

	// Distribute attributes by the processes
	// To avoid having to share words between processes we distribute
	// blocks of WORD_BITS attributes for each process
	dm.a_offset = BLOCK_LOW(rank, size, dataset.n_words);
	dm.a_size	= BLOCK_SIZE(rank, size, dataset.n_words);

	dm.first_attribute = dm.a_offset * WORD_BITS;

	dm.n_attributes = MIN(dm.a_size * WORD_BITS, dataset.n_attributes);

	if (dm.n_attributes > 0)
	{
		dm.last_attribute = dm.first_attribute + dm.n_attributes - 1;
	}
	else
	{
		dm.last_attribute = dm.first_attribute;
	}

	TOCK;

	if (rank == ROOT_RANK)
	{
		ROOT_SHOWS("  Number of lines in the disjoint matrix: %lu\n",				   dm.n_matrix_lines);

		double matrix_size = dm.n_matrix_lines * dataset.n_attributes;
		matrix_size /= (1024.0 * 1024.0 * 1024.0 * 8.0);
		ROOT_SHOWS("  Estimated disjoint matrix size: %0.2lfGB\n", matrix_size);

		for (int r = 0; r < size; r++)
		{
			uint32_t a_offset = BLOCK_LOW(r, size, dataset.n_words);
			uint32_t a_size	  = BLOCK_SIZE(r, size, dataset.n_words);

			uint32_t first_attribute = a_offset * WORD_BITS;
			uint32_t n_attributes
				= MIN(a_size * WORD_BITS, dataset.n_attributes);
			uint32_t last_attribute = dm.first_attribute;
			if (n_attributes > 0)
			{
				last_attribute = first_attribute + n_attributes - 1;
			}

			if (a_size > 0)
			{
				fprintf(stdout,
						"    Process %d will manage %u attributes [%u -> %u]\n",
						r, n_attributes, first_attribute, last_attribute);
			}
			else
			{
				fprintf(stdout, "    Process %d will manage 0 attributes\n", r);
			}
		}
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

	ROOT_SAYS("Applying set covering algorithm:\n");
	TICK;

	/**
	 * Number of words needed to store a column (attribute data)
	 */
	dm.n_words_in_a_column
		= dm.n_matrix_lines / WORD_BITS + (dm.n_matrix_lines % WORD_BITS != 0);

	/**
	 * Bit array with the information about the lines covered (1) or not (0) by
	 * the current best attribute
	 */
	word_t* best_column
		= (word_t*) calloc(dm.n_words_in_a_column, sizeof(word_t));

	/**
	 * The number of attributes is rounded so we can check all bits during
	 * the attribute totals calculation, avoiding one if in the inner cycle
	 */
	uint32_t* attribute_totals
		= (uint32_t*) calloc(dm.a_size * WORD_BITS, sizeof(uint32_t));

	// The covered lines and selected attributes are bit arrays
	/**
	 * Bit array with the information about the lines already covered (1) or
	 * not (0) so far
	 */
	word_t* covered_lines
		= (word_t*) calloc(dm.n_words_in_a_column, sizeof(word_t));

	/**
	 * Selected attributes aka the solution, only root needs them
	 */
	word_t* selected_attributes = NULL;
	if (rank == ROOT_RANK)
	{
		selected_attributes = (word_t*) calloc(dataset.n_words, sizeof(word_t));
	}

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
	 * Number of uncovered lines remaining
	 */
	uint32_t n_uncovered_lines = dm.n_matrix_lines;

	// Calculate initial values
	calculate_attribute_totals_add(&dataset, &dm, covered_lines,
								   attribute_totals);

	while (true)
	{
		// What is my best attribute?
		best_attribute_t local_max
			= get_best_attribute(attribute_totals, dm.n_attributes);

		// My best attribute translates to which global attribute?
		if (local_max.attribute != -1)
		{
			local_max.attribute += dm.first_attribute;
		}

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

		// Mark attributed as selected
		if (rank == ROOT_RANK)
		{
			ROOT_SHOWS("  Selected attribute #%ld, ", best_max.attribute);
			ROOT_SHOWS("covers %ld lines ", best_max.n_covered_lines);
			TOCK;
			TICK;

			// Mark best attribute as selected
			mark_attribute_as_selected(best_max.attribute, selected_attributes);
		}

		// Update number of uncovered lines
		n_uncovered_lines -= best_max.n_covered_lines;

		// If we have everything covered goto show_solution
		if (n_uncovered_lines == 0)
		{
			goto show_solution;
		}

		// Get the data for the best attribute
		get_column(&dataset, &dm, best_max.attribute, best_column);

		/**
		 * If this attribute covers more lines than what remains
		 * to be covered we calculate the sum of the remaining lines.
		 * If the attribute covers only a few lines we remove the
		 * contribution of the covered lines.
		 * The objetive is to reduce the number of lines generated.
		 */
		if (n_uncovered_lines < best_max.n_covered_lines)
		{

			// Update covered lines
			update_covered_lines(best_column, dm.n_words_in_a_column,
								 covered_lines);

			// Update attributes totals
			calculate_attribute_totals_add(&dataset, &dm, covered_lines,
										   attribute_totals);
		}
		else
		{
			// Update best attribute data, leaving only the new covered lines
			for (uint32_t w = 0; w < dm.n_words_in_a_column; w++)
			{
				best_column[w] &= ~covered_lines[w];
			}

			// Update attributes totals
			calculate_attribute_totals_sub(&dataset, &dm, best_column,
										   attribute_totals);

			// Update covered lines
			update_covered_lines(best_column, dm.n_words_in_a_column,
								 covered_lines);
		}
	}

show_solution:
	if (rank == ROOT_RANK)
	{
		printf("Solution: { ");

		uint32_t current_attribute = 0;
		uint32_t solution_size	   = 0;

		for (uint32_t w = 0; w < dataset.n_words; w++)
		{
			for (int8_t bit = WORD_BITS - 1;
				 bit >= 0 && current_attribute < dataset.n_attributes;
				 bit--, current_attribute++)
			{
				if (selected_attributes[w] & AND_MASK_TABLE[bit])
				{
					// This attribute is set so it's part of the solution
					printf("%d ", current_attribute);
					solution_size++;
				}
			}
		}

		printf("}\nSolution has %d attributes: %d / %d = %3.4f%%\n",
			   solution_size, solution_size, dataset.n_attributes,
			   ((float) solution_size / (float) dataset.n_attributes) * 100);

		fprintf(stdout, "}\nAll done! ");
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

	PRINT_TIMING_GLOBAL;

	/* shut down MPI */
	MPI_Finalize();

	return EXIT_SUCCESS;
}
