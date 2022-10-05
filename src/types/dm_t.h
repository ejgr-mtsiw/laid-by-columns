/*
 ============================================================================
 Name        : dm_t.h
 Author      : Eduardo Ribeiro
 Description : Datatype representing one disjoint matriz or part of one.
			   Each column corresponds to one attribute
 ============================================================================
 */

#ifndef DM_T_H
#define DM_T_H

#include "../types/steps_t.h"

#include <stdint.h>

typedef struct dm_t
{
	/**
	 * The number of lines of the full matrix
	 */
	uint32_t n_matrix_lines;

	/**
	 * The first attribute of this process
	 */
	uint32_t a_offset;

	/**
	 * Number of attributes
	 */
	uint32_t a_size;

	/**
	 * Steps to generate the lines of the disjoint matrix
	 */
	steps_t* steps;
} dm_t;

#endif // DM_T_H
