/*
 ============================================================================
 Name        : types/best_attribute_t.h
 Author      : Eduardo Ribeiro
 Description : Datatype to store best attributes
 ============================================================================
 */

#ifndef TYPES_BEST_ATTRIBUTE_T_H
#define TYPES_BEST_ATTRIBUTE_T_H

#include <stdint.h>

typedef struct best_attribute_t
{
	/**
	 * Attribute
	 */
	long attribute;

	/**
	 * Number of covered lines
	 */
	uint64_t n_covered_lines;
} best_attribute_t;

#endif // TYPES_BEST_ATTRIBUTE_T_H
