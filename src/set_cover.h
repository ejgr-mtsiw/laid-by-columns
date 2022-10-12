/*
 ============================================================================
 Name        : set_cover.h
 Author      : Eduardo Ribeiro
 Description : Structures and functions to apply the set cover algorithm
 ============================================================================
 */

#ifndef SET_COVER_H
#define SET_COVER_H

#include "types/best_attribute_t.h"
#include "types/dataset_t.h"
#include "types/dm_t.h"

#include <stdint.h>

/**
 * Searches the attribute totals array for the highest score and returns the
 * correspondent best attribute.
 * Return has index -1 if there are no more attributes available.
 */
best_attribute_t get_best_attribute(const uint32_t* totals,
									const uint32_t n_attributes);

///**
// * Calculates the initial totals for all attributes
// */
//void calculate_initial_totals(dm_t* dm, dataset_t* dataset, uint32_t* totals);
//
///**
// * Adds the contribution from lxor to the atribute totals starting at attribute
// * start
// */
//void add_to_totals(uint32_t* totals, uint32_t* start, word_t lxor);
//
///**
// * Removes the contribution of the covered lines from the totals array
// */
//void update_totals(dm_t* dm, dataset_t* dataset, uint32_t* totals,
//				   word_t* covered_lines, word_t* best_column);
//
///**
// * Removes contribution from lxor to the attribute totals starting at attribute
// * start
// */
//void sub_from_totals(uint32_t* totals, uint32_t* start, word_t lxor);

#endif
