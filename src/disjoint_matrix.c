/*
 ============================================================================
 Name        : disjoint_matrix.c
 Author      : Eduardo Ribeiro
 Description : Structures and functions to manage the disjoint matrix
 ============================================================================
 */

#include "disjoint_matrix.h"

#include "types/dataset_t.h"
#include "types/dm_t.h"

#include <stdint.h>
#include <stdlib.h>

uint64_t get_dm_n_lines(const dataset_t* dataset)
{
	// Calculate number of lines for the matrix
	uint64_t n = 0;

	uint64_t n_classes	  = dataset->n_classes;
	uint64_t* n_class_obs = dataset->n_observations_per_class;

	for (uint64_t class_a = 0; class_a < n_classes - 1; class_a++)
	{
		for (uint64_t class_b = class_a + 1; class_b < n_classes; class_b++)
		{
			n += (uint64_t)n_class_obs[class_a] * (uint64_t)n_class_obs[class_b];
		}
	}

	return n;
}
