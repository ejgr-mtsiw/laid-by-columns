/*
 ============================================================================
 Name        : jnsq.c
 Author      : Eduardo Ribeiro
 Description : Structures and functions to manage JNSQ
 ============================================================================
 */

#include "jnsq.h"

#include "dataset.h"
#include "types/dataset_t.h"
#include "types/word_t.h"
#include "utils/bit.h"

#include <math.h>
#include <stdint.h>

#include <stdio.h>

void set_jnsq_bits(word_t* line, uint64_t inconsistency,
				   const uint64_t n_attributes, const uint64_t n_words,
				   const uint8_t n_bits_for_class)
{

	// Words that have attributes
	uint64_t true_n_words =  n_attributes/ WORD_BITS + (n_attributes% WORD_BITS != 0);

	// The first word with jnsq
	uint64_t word_with_jnsq = true_n_words-1+n_words*0;

	// How many attributes remain on last word with attributes
	uint8_t remaining_bits = n_attributes % WORD_BITS;

	// n jnsq bits needed
	uint8_t n_bits_needed = n_bits_for_class;

	//printf("nw[%lu] tnw[%lu] wwjnsq[%lu] rb[%d] nbn[%d] i[%lu]\n", n_words, true_n_words, word_with_jnsq, remaining_bits, n_bits_needed, inconsistency);

	if (remaining_bits+n_bits_needed >WORD_BITS) {
		// We must split the jnsq attributes

		// n jnsq bits on first word
		uint8_t n_bits = WORD_BITS - remaining_bits;

		// Invert consistency
		inconsistency = invert_n_bits((word_t) inconsistency, n_bits);

		line[word_with_jnsq]=set_bits(line[word_with_jnsq], inconsistency, 0, n_bits);

		// There's no more attributes
		remaining_bits=0;

		// Remove used bits from inconsistency
		inconsistency >>= n_bits;

		// n jnsq bits on next word
		n_bits_needed -= n_bits;

		// next word_with_jnsq
		word_with_jnsq++;
	}

	// All remaining jnsq bits are in the same word
	uint8_t jnsq_start = WORD_BITS - remaining_bits - n_bits_needed;

	// Invert consistency
	inconsistency = invert_n_bits((word_t) inconsistency, n_bits_needed);

	line[word_with_jnsq] = set_bits(line[word_with_jnsq], inconsistency, jnsq_start, n_bits_needed);
}

//void update_jnsq(word_t* to_update, const word_t* to_compare,
//				 uint64_t* inconsistency, const uint64_t n_attributes,
//				 const uint64_t n_words, const uint8_t n_bits_for_class)
//{
//	// Set the line JNSQ
//	set_jnsq_bits(to_update, (*inconsistency), n_attributes, n_words,
//				  n_bits_for_class);
//
//	if (has_same_attributes(to_update, to_compare, n_attributes))
//	{
//		// Inconsistency!
//		(*inconsistency)++; // Because observations are sorted by class
//	}
//	else
//	{
//		// Differente attributes - reset JNSQ
//		(*inconsistency) = 0;
//	}
//}

/**
 * Adds the JNSQs attributes to the dataset
 */
uint64_t add_jnsqs(dataset_t* dataset)
{
	// Current line
	word_t* current = dataset->data;

	// Previous line
	word_t* prev = current;

	// Number of attributes
	uint64_t n_attributes = dataset->n_attributes;

	// Number of longs in a line
	uint64_t n_words = dataset->n_words;

	// Number of observations in the dataset
	uint64_t n_observations = dataset->n_observations;

	// Number of bits needed to store class
	uint8_t n_bits_for_class = dataset->n_bits_for_class;

	// Last line
	word_t* last = GET_LAST_LINE(dataset->data, n_observations, n_words);

	// Inconsistency
	uint64_t inconsistency = 0;

	// Max inconsistency found
	uint64_t max_inconsistency = 0;

	// first line has jnsq=0
	set_jnsq_bits(current, 0, n_attributes, n_words, n_bits_for_class);

	// Now do the others
	for (prev = dataset->data; prev < last; NEXT_LINE(prev, n_words))
	{
		NEXT_LINE(current, n_words);

		if (has_same_attributes(current, prev, n_attributes))
		{
			// Inconsistency!
			// Because observations with the same attributes end up sorted by class
			inconsistency++;

			// Update max
			if (inconsistency > max_inconsistency)
			{
				max_inconsistency = inconsistency;
			}
		}
		else
		{
			// Differente attributes - reset JNSQ
			inconsistency = 0;
		}

		// Set the line JNSQ
		set_jnsq_bits(current, inconsistency, n_attributes, n_words, n_bits_for_class);
	}

	return max_inconsistency;
}
