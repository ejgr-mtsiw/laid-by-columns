/*
 ============================================================================
 Name        : mpi_custom.c
 Author      : Eduardo Ribeiro
 Description : Custom function to select best global attribute in MPI_Allreduce
 ============================================================================
 */

#include "types/best_attribute_t.h"
#include "utils/mpi_custom.h"

void MPI_get_best_attribute(void* in, void* inout, int* len, MPI_Datatype* dptr)
{
	// Suppress unused warning/error
	(void) len;
	(void) dptr;

	best_attribute_t* ina	 = (best_attribute_t*) in;
	best_attribute_t* inouta = (best_attribute_t*) inout;

	if (ina->max > inouta->max)
	{
		inouta->max		  = ina->max;
		inouta->attribute = ina->attribute;
	}
}
