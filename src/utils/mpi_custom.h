/*
 ============================================================================
 Name        : utils/mpi.h
 Author      : Eduardo Ribeiro
 Description : Custom function to select best global attribute in MPI_Allreduce
 ============================================================================
 */

#ifndef UTILS_MPI_CUSTOM_H
#define UTILS_MPI_CUSTOM_H

#include "mpi.h"

/**
 * Custom function to select best global attribute in MPI_Allreduce
 */
void MPI_get_best_attribute(void* in, void* inout, int* len, MPI_Datatype* dptr);

#endif // UTILS_MPI_CUSTOM_H
