/*
 ============================================================================
 Name        : mpi_custom.c
 Author      : Eduardo Ribeiro
 Description : Custom function to select best global attribute in MPI_Allreduce
 ============================================================================
 */

#include "utils/mpi_custom.h"

#include "types/best_attribute_t.h"

#include "mpi.h"

void MPI_get_best_attribute(void* in, void* inout, int* len, MPI_Datatype* dptr)
{
	// Suppress unused warning/error
	(void) len;
	(void) dptr;

	best_attribute_t* ina	 = (best_attribute_t*) in;
	best_attribute_t* inouta = (best_attribute_t*) inout;

	if (ina->n_covered_lines > inouta->n_covered_lines)
	{
		inouta->n_covered_lines = ina->n_covered_lines;
		inouta->attribute		= ina->attribute;
	}
	//// Right now the program selects the first attribute sent by the
	//// first process to terminate, uncomment next lines to revert to
	//// standard greedy behaviour selecting the lowest attribute
	//
	 else if (ina->n_covered_lines == inouta->n_covered_lines
			 && inouta->attribute > ina->attribute)
	{
		inouta->attribute = ina->attribute;
	}
}
