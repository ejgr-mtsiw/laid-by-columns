/*
 ============================================================================
 Name        : types/best_attribute_t.h
 Author      : Eduardo Ribeiro
 Description : Datatype to store best attributes
 ============================================================================
 */

#ifndef TYPES_BEST_ATTRIBUTE_T_H
#define TYPES_BEST_ATTRIBUTE_T_H

typedef struct best_attribute_t
{
	/**
	 * Attribute
	 */
	long attribute;

	/**
	 * Number of covered lines
	 */
	long max;
} best_attribute_t;

#endif // TYPES_BEST_ATTRIBUTE_T_H
