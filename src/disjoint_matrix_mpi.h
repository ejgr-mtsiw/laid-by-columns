/*
 ============================================================================
 Name        : disjoint_matrix_mpi.h
 Author      : Eduardo Ribeiro
 Description : Structures and functions to manage the disjoint matrix in MPIIO
 ============================================================================
 */

#ifndef MPI_DISJOINT_MATRIX_H
#define MPI_DISJOINT_MATRIX_H

#include "types/dataset_t.h"
#include "types/dm_t.h"
#include "types/oknok_t.h"

oknok_t get_column(const dataset_t* dataset, const dm_t* dm,
				   const uint32_t index, word_t* column);

#endif // MPI_DISJOINT_MATRIX_H
