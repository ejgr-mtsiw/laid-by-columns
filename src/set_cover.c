/*
 ============================================================================
 Name        : set_cover.c
 Author      : Eduardo Ribeiro
 Description : Structures and functions to apply the set cover algorithm
 ============================================================================
 */

#include "set_cover.h"

#include "types/best_attribute_t.h"
#include "types/dataset_t.h"
#include "types/dm_t.h"
#include "types/oknok_t.h"
#include "utils/bit.h"

#include <stdint.h>
#include <string.h>
#include <stdio.h>

best_attribute_t get_best_attribute(const uint64_t* totals,
									const uint64_t n_attributes)
{
	uint64_t max_total	  = 0;
	int64_t max_attribute = -1;

	for (uint64_t i = 0; i < n_attributes; i++)
	{
		if (totals[i] > max_total)
		{
			max_total	  = totals[i];
			max_attribute = i;
		}
	}

	best_attribute_t best
		= { .attribute = max_attribute, .n_covered_lines = max_total };
	return best;
}

oknok_t mark_attribute_as_selected(const int64_t attribute,
								   word_t* selected_attributes)
{

	// Which word has the best attribute
	uint64_t best_word = attribute / WORD_BITS;

	// Which bit?
	uint8_t best_bit = WORD_BITS - attribute % WORD_BITS - 1;

	// Mark best attribute as selected
	BIT_SET(selected_attributes[best_word], best_bit);

	return OK;
}

oknok_t calculate_attribute_totals_add(const dataset_t* dataset, const dm_t* dm,
									   const word_t* covered_lines,
									   uint64_t* attribute_totals)
{

	uint64_t nw = dataset->n_words;
	uint64_t nc	   = dataset->n_classes;
	line_class_t*classes = dataset->classes;

	uint64_t cl = 0;
	uint64_t ca = 0;
	uint64_t ia = 0;
	uint64_t cb = 0;
	uint64_t ib = 0;

//	printf("nw %lu, nc %lu, as %lu, ao %lu\n", nw, nc, dm->a_size, dm->a_offset);

	// Reset totals
	memset(attribute_totals, 0, dm->a_size * WORD_BITS * sizeof(uint64_t));

	for (uint64_t ww = dm->a_offset; ww < dm->a_offset + dm->a_size; ww+=8)
	{
		cl=0;
		for (ca = 0; ca < nc - 1; ca++)
		{
			for (cb = ca + 1; cb < nc; cb++)
			{
				for (ia = 0; ia < classes[ca].n_observations; ia++)
				{
					for (ib = 0; ib < classes[cb].n_observations; ib++, cl++)
					{

						uint64_t w = cl / WORD_BITS;
						uint8_t b  = WORD_BITS - (cl % WORD_BITS) - 1;

						if (!BIT_CHECK(covered_lines[w], b))
						{
							// This line is not yet covered
							// Generate the partial line
							word_t* la = classes[ca].first_observation_address + ia*nw;
							word_t* lb = classes[cb].first_observation_address + ib*nw;

							uint64_t c_attribute = (ww-dm->a_offset)*WORD_BITS;
							//printf ("cl %lu, ww %lu, ca %lu\n", cl, ww, c_attribute);

							for (uint64_t www = ww; www < ww+8; www++)
							{
								word_t attributes = la[www] ^ lb[www];

								for (int8_t bit = WORD_BITS - 1; bit >= 0; bit--, c_attribute++)
								{
									attribute_totals[c_attribute] += BIT_CHECK(attributes, bit);
								}
							}
						}
					}
				}
			}
		}
	}

	return OK;
}

oknok_t calculate_attribute_totals_sub(const dataset_t* dataset, const dm_t* dm,
									   const word_t* covered_lines,
									   uint64_t* attribute_totals)
{

	uint64_t nw = dataset->n_words;
	uint64_t nc	   = dataset->n_classes;
	line_class_t*classes = dataset->classes;

	uint64_t cl = 0;
	uint64_t ca = 0;
	uint64_t ia = 0;
	uint64_t cb = 0;
	uint64_t ib = 0;

	for (uint64_t ww = dm->a_offset; ww < dm->a_offset + dm->a_size; ww+=8)
		{
			cl=0;
			for (ca = 0; ca < nc - 1; ca++)
			{
				for (cb = ca + 1; cb < nc; cb++)
				{
					for (ia = 0; ia < classes[ca].n_observations; ia++)
					{
						for (ib = 0; ib < classes[cb].n_observations; ib++, cl++)
						{

							uint64_t w = cl / WORD_BITS;
							uint8_t b  = WORD_BITS - (cl % WORD_BITS) - 1;

							if (!BIT_CHECK(covered_lines[w], b))
							{
								// This line is not yet covered
								// Generate the partial line
								word_t* la = classes[ca].first_observation_address + ia*nw;
								word_t* lb = classes[cb].first_observation_address + ib*nw;

								uint64_t c_attribute = (ww-dm->a_offset)*WORD_BITS;
								//printf ("cl %lu, ww %lu, ca %lu\n", cl, ww, c_attribute);

								for (uint64_t www = ww; www < ww+8; www++)
								{
									word_t attributes = la[www] ^ lb[www];

									for (int8_t bit = WORD_BITS - 1; bit >= 0; bit--, c_attribute++)
									{
										attribute_totals[c_attribute] -= BIT_CHECK(attributes, bit);
									}
								}
							}
						}
					}
				}
			}
		}

	return OK;
}

oknok_t update_covered_lines(const word_t* best_column,
							 const uint64_t n_words_in_a_column,
							 word_t* covered_lines)
{
	for (uint64_t w = 0; w < n_words_in_a_column; w++)
	{
		covered_lines[w] |= best_column[w];
	}

	return OK;
}
