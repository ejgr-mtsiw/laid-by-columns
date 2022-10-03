/*
 ============================================================================
 Name        : set_cover.c
 Author      : Eduardo Ribeiro
 Description : Structures and functions to apply the set cover algorithm
 ============================================================================
 */

#include "set_cover.h"

#include <stdint.h>

int64_t get_best_attribute_index(const uint32_t* totals,
								 const uint32_t n_attributes)
{
	uint32_t max_total	  = 0;
	int64_t max_attribute = -1;

	for (uint32_t i = 0; i < n_attributes; i++)
	{
		if (totals[i] > max_total)
		{
			max_total	  = totals[i];
			max_attribute = i;
		}
	}

	return max_attribute;
}
