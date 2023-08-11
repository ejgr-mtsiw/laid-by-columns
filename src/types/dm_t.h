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

#include <stdint.h>

typedef struct dm_t
{
	/**
	 * The number of lines of the full matrix
	 */
	uint64_t n_matrix_lines;

	/**
	 * The number of words to store a full column of the matrix
	 */
	uint64_t n_words_in_a_column;

	/**
	 * The number of attributes in the partial disjoint matrix
	 */
	uint32_t n_attributes;

	/**
	 * The first attribute of this process
	 */
	uint32_t first_attribute;

	/**
	 * The last attribute os this process
	 */
	uint32_t last_attribute;

	/**
	 * The first attribute we can generate
	 */
	uint32_t a_offset;

	/**
	 * Number of attributes we can generate
	 */
	uint32_t a_size;
} dm_t;

#endif // DM_T_H
