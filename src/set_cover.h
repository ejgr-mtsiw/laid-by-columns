/*
 ============================================================================
 Name        : set_cover.h
 Author      : Eduardo Ribeiro
 Description : Structures and functions to apply the set cover algorithm
 ============================================================================
 */

#ifndef SET_COVER_H
#define SET_COVER_H

#include <stdint.h>

/**
 * Searches the attribute totals array for the highest score and returns the
 * correspondent attribute index.
 * Returns -1 if there are no more attributes available.
 */
int64_t get_best_attribute_index(const uint32_t* totals,
								 const uint32_t n_attributes);

#endif
