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
#include "types/oknok_t.h"
#include "types/word_t.h"

#include <stdint.h>

/**
 * Number of words to process per disjoint matrix cycle
 */
#define N_WORDS_PER_CYCLE 8

/**
 * Searches the attribute totals array for the highest score and returns the
 * correspondent best attribute.
 * Return has index -1 if there are no more attributes available.
 */
best_attribute_t get_best_attribute(const uint32_t* totals,
									const uint32_t n_attributes);

oknok_t mark_attribute_as_selected(const int64_t attribute,
								   word_t* selected_attributes);

/**
 * Calculates the current attributes totals
 */
oknok_t calculate_attribute_totals_add(const dataset_t* dataset, const dm_t* dm,
									   const word_t* covered_lines,
									   uint32_t* attribute_totals);

oknok_t calculate_attribute_totals_sub(const dataset_t* dataset, const dm_t* dm,
									   const word_t* covered_lines,
									   uint32_t* attribute_totals);

oknok_t update_covered_lines(const word_t* best_column,
							 const uint32_t n_words_in_a_column,
							 word_t* covered_lines);
#endif
