/*
 ============================================================================
 Name        : laid-hdf5-mpi.c
 Author      : Eduardo Ribeiro
 Description : OpenMPI implementation of the LAID algorithm in C + HDF5
 ============================================================================
 */

#include "clargs.h"
#include "dataset.h"
#include "disjoint_matrix.h"
#include "jnsq.h"
#include "mpi_disjoint_matrix.h"
#include "mpi_hdf5_dataset.h"
#include "sort_r.h"
#include "timing.h"
#include "word_t.h"

#include "mpi.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

/**
 * Reads dataset attributes from hdf5 file
 * Read dataset
 * Sort dataset
 * Remove duplicates
 * Add jnsqs
 * Write disjoint matrix
 * Apply set covering algorithm
 * Show solution
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

	// rank of process
	int rank;
	// number of processes
	int size;

	// Global communicator group
	MPI_Comm comm = MPI_COMM_WORLD;

	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	MPI_Comm node_comm = MPI_COMM_NULL;

	// Create node-local communicator
	MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL,
						&node_comm);

	int node_size, node_rank;

	MPI_Comm_size(node_comm, &node_size);
	MPI_Comm_rank(node_comm, &node_rank);

	// Open dataset
	hdf5_dataset_t hdf5_dset;
	mpi_hdf5_open_dataset(args.filename, args.datasetname, comm, MPI_INFO_NULL,
						  &hdf5_dset);

	uint32_t n_observations = hdf5_dset.dimensions[0];
	uint32_t n_words		= hdf5_dset.dimensions[1];

	struct timespec main_tick, main_tock;
	if (rank == 0)
	{
		// Timing for the full operation
		clock_gettime(CLOCK_MONOTONIC_RAW, &main_tick);
	}

	// We can jump straight to the set covering algorithm
	// if we already have the matrix in the hdf5 dataset
	// Only root checks to avoid having all processes trying
	// to open and read the file at the same time
	uint8_t skip_matrix_creation = 0;
	if (rank == 0)
	{
		skip_matrix_creation = is_matrix_created(args.filename);

		printf("- Disjoint matrix dataset %sfound!.\n",
			   skip_matrix_creation ? "" : "not ");
	}

	MPI_Bcast(&skip_matrix_creation, 1, MPI_INT8_T, 0, comm);

	if (!skip_matrix_creation)
	{
		// SETUP_TIMING
		// TICK;

		// Only rank 0 on a node actually allocates memory
		uint64_t localtablesize = 0;
		if (node_rank == 0)
		{
			localtablesize = n_observations * n_words;
		}

		char node_name[MPI_MAX_PROCESSOR_NAME];
		int node_str_len = 0;
		MPI_Get_processor_name(node_name, &node_str_len);

		// debug info
		//		printf(
		//			"Rank %d of %d, rank %d of %d in node <%s>, localtablesize
		//%lu\n", 			rank, size, node_rank, node_size, node_name,
		// localtablesize);

		word_t* localtable		= NULL;
		MPI_Win win_shared_dset = MPI_WIN_NULL;
		MPI_Win_allocate_shared(localtablesize * sizeof(word_t), sizeof(word_t),
								MPI_INFO_NULL, node_comm, &localtable,
								&win_shared_dset);

		/**
		 * The dataset
		 */
		dataset_t dataset;

		init_dataset(&dataset);

		// Set dataset data pointer
		if (node_rank == 0)
		{
			dataset.data = localtable;
		}
		else
		{
			MPI_Aint win_size;
			int win_disp;
			MPI_Win_shared_query(win_shared_dset, 0, &win_size, &win_disp,
								 &dataset.data);
		}

		SETUP_TIMING

		// All table pointers should now point to copy on noderank 0
		// Initialise table on rank 0 with appropriate synchronisation
		MPI_Win_fence(0, win_shared_dset);

		if (node_rank == 0)
		{
			TICK;

			fprintf(stdout, "- Loading dataset data\n - ");

			// Load dataset attributes
			hdf5_read_dataset_attributes(hdf5_dset.dataset_id, &dataset);
			// Load dataset data
			hdf5_read_data(hdf5_dset.dataset_id, &dataset);

			print_dataset_details(stdout, &dataset);

			fprintf(stdout, "- Finished loading dataset data ");

			TOCK(stdout)
			TICK;

			// Sort dataset
			fprintf(stdout, "- Sorting dataset\n");

			// We need to know the number of longs in each line of the dataset
			// so we can't use the standard qsort implementation
			sort_r(dataset.data, dataset.n_observations,
				   dataset.n_words * sizeof(word_t), compare_lines_extra,
				   &dataset.n_words);

			fprintf(stdout, "- Sorted dataset");
			TOCK(stdout)
			TICK;

			// Remove duplicates
			fprintf(stdout, "- Removing duplicates:\n");

			unsigned int duplicates = remove_duplicates(&dataset);

			fprintf(stdout, " - %d duplicate(s) removed ", duplicates);
			TOCK(stdout)
			TICK;

			// Fill class arrays
			fprintf(stdout, "- Checking classes: ");

			if (fill_class_arrays(&dataset) != OK)
			{
				return EXIT_FAILURE;
			}

			TOCK(stdout)

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
			TOCK(stdout)
		}

		MPI_Win_fence(0, win_shared_dset);

		// Only rank 0 on a node actually allocates memory
		uint64_t n_matrix_lines = 0;
		if (node_rank == 0)
		{
			n_matrix_lines = get_dm_n_lines(&dataset);
		}

		// debug info
		//		printf("Rank %d of %d, rank %d of %d in node <%s>, local_dm_size
		//%lu\n", 			   rank, size, node_rank, node_size, node_name,
		// n_matrix_lines);

		steps_t* localsteps		 = NULL;
		MPI_Win win_shared_steps = MPI_WIN_NULL;
		MPI_Win_allocate_shared(n_matrix_lines * sizeof(steps_t),
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
			MPI_Win_shared_query(win_shared_steps, 0, &win_size, &win_disp,
								 &steps);
		}

		// All table pointers should now point to copy on noderank 0
		// Initialise table on rank 0 with appropriate synchronisation
		MPI_Win_fence(0, win_shared_steps);

		if (node_rank == 0)
		{
			SETUP_TIMING
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

			//			for (uint64_t i=0;i<n_matrix_lines;i++){
			//				printf("[%lu]: %d ^ %d\n", i, steps[i].indexA+1,
			// steps[i].indexB+1);
			//			}

			fprintf(stdout, "- Finished generating matrix steps ");
			TOCK(stdout)
		}

		MPI_Win_fence(0, win_shared_steps);

		if (rank == 0)
		{
			fprintf(stdout, "- Broadcasting number of matrix lines\n");
		}

		MPI_Bcast(&n_matrix_lines, 1, MPI_UINT64_T, 0, comm);

		if (rank == 0)
		{
			fprintf(stdout, "- Finished broadcasting number of matrix lines\n");
		}

		if (rank == 0)
		{
			fprintf(stdout, "- Building disjoint matrix\n");
		}

		if (rank == 0)
		{
			TICK;
		}

		// Build part of the disjoint matrix and store it in the hdf5 file
		mpi_create_line_dataset(hdf5_dset, dataset.data, n_matrix_lines,
								n_words, steps, rank, size);

		MPI_Barrier(comm);
		if (rank == 0)
		{
			fprintf(stdout, "Finished building disjoint matrix ");
			TOCK(stdout)
		}

		MPI_Win_free(&win_shared_dset);
		MPI_Win_free(&win_shared_steps);

		dataset.data = NULL;
		free_dataset(&dataset);

		hdf5_close_dataset(&hdf5_dset);

		//		unsigned long matrix_lines = get_dm_n_lines(&dataset);
		//
		//		double matrixsize
		//			= (matrix_lines * (dataset.n_attributes +
		// dataset.n_bits_for_class)) 			/ (1024 * 1024 * 8);
		//
		//		fprintf(stdout, "\nBuilding disjoint matrix.\n");
		//		fprintf(stdout,
		//				"Estimated disjoint matrix size: %lu lines
		//[%0.2fMB]\n", 				matrix_lines, matrixsize);

		// Sync everyone before starting building the matrix dataset
		// MPI_Barrier(comm);

		// Build part of the disjoint matrix and store it in the hdf5 file
		//		mpi_create_disjoint_matrix(args.filename, &dataset, rank, size,
		// comm, 								   MPI_INFO_NULL);

		//		fprintf(stdout, "Finished building disjoint matrix ");
		//		TOCK(stdout)

		/**
		 * From this point forward we no longer need the dataset
		 */
		// free_dataset(&dataset);
	}
	//	//
	//	//	cover_t cover;
	//	//	init_cover(&cover);
	//	//
	//	//	TICK
	//	//	fprintf(stdout, "\nApplying set covering algorithm.\n");
	//	//	if (calculate_solution(args.filename, &cover) != OK) {
	//	//		return EXIT_FAILURE;
	//	//	}
	//	//
	//	//	print_solution(stdout, &cover);
	//	//
	//	//	free_cover(&cover);
	//	//	TOCK(stdout)

	// the_end:
	//  Wait for everyone
	MPI_Barrier(comm);

	if (rank == 0)
	{
		fprintf(stdout, "All done! ");

		clock_gettime(CLOCK_MONOTONIC_RAW, &main_tock);
		fprintf(stdout, "[%0.3fs]\n",
				(main_tock.tv_nsec - main_tick.tv_nsec) / 1000000000.0
					+ (main_tock.tv_sec - main_tick.tv_sec));
	}

	/* shut down MPI */
	MPI_Finalize();

	return EXIT_SUCCESS;
}