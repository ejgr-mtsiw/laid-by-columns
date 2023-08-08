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

best_attribute_t get_best_attribute(const uint32_t* totals,
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

	best_attribute_t best
		= { .attribute = max_attribute, .n_covered_lines = max_total };
	return best;
}

oknok_t mark_attribute_as_selected(const int64_t attribute,
								   word_t* selected_attributes)
{

	// Which word has the best attribute
	uint32_t best_word = attribute / WORD_BITS;

	// Which bit?
	uint8_t best_bit = WORD_BITS - attribute % WORD_BITS - 1;

	// Mark best attribute as selected
	BIT_SET(selected_attributes[best_word], best_bit);

	return OK;
}

oknok_t calculate_attribute_totals_add(const dataset_t* dataset, const dm_t* dm,
									   const word_t* covered_lines,
									   uint32_t* attribute_totals)
{

	uint32_t nc	   = dataset->n_classes;
	uint32_t nobs  = dataset->n_observations;
	word_t** opc   = dataset->observations_per_class;
	uint32_t* nopc = dataset->n_observations_per_class;

	uint32_t cs = 0;
	uint32_t ca = 0;
	uint32_t ia = 0;
	uint32_t cb = 0;
	uint32_t ib = 0;

	// Reset totals
	memset(attribute_totals, 0, dm->a_size * WORD_BITS * sizeof(uint32_t));

	for (ca = 0; ca < nc - 1; ca++)
	{
		for (cb = ca + 1; cb < nc; cb++)
		{
			for (ia = 0; ia < nopc[ca]; ia++)
			{
				for (ib = 0; ib < nopc[cb]; ib++, cs++)
				{
					uint32_t w = cs / WORD_BITS;
					uint8_t b  = WORD_BITS - (cs % WORD_BITS) - 1;

					if (!BIT_CHECK(covered_lines[w], b))
					{
						// This line is not yet covered

						uint32_t c_attribute = 0;
						for (uint32_t ww = dm->a_offset;
							 ww < dm->a_offset + dm->a_size; ww++)
						{
							// Generate next line
							word_t** bla = opc + ca * nobs;
							word_t** blb = opc + cb * nobs;

							word_t* la = *(bla + ia);
							word_t* lb = *(blb + ib);

							word_t attributes = la[ww] ^ lb[ww];

							for (int8_t bit = WORD_BITS - 1; bit >= 0;
								 bit--, c_attribute++)
							{
								attribute_totals[c_attribute]
									+= BIT_CHECK(attributes, bit);
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
									   uint32_t* attribute_totals)
{

	uint32_t nc	   = dataset->n_classes;
	uint32_t nobs  = dataset->n_observations;
	word_t** opc   = dataset->observations_per_class;
	uint32_t* nopc = dataset->n_observations_per_class;

	uint32_t cs = 0;
	uint32_t ca = 0;
	uint32_t ia = 0;
	uint32_t cb = 0;
	uint32_t ib = 0;

	for (ca = 0; ca < nc - 1; ca++)
	{
		for (cb = ca + 1; cb < nc; cb++)
		{
			for (ia = 0; ia < nopc[ca]; ia++)
			{
				for (ib = 0; ib < nopc[cb]; ib++, cs++)
				{
					uint32_t w = cs / WORD_BITS;
					uint8_t b  = WORD_BITS - (cs % WORD_BITS) - 1;

					if (BIT_CHECK(covered_lines[w], b))
					{
						// This line wasn't covered but it is now

						uint32_t c_attribute = 0;
						for (uint32_t ww = dm->a_offset;
							 ww < dm->a_offset + dm->a_size; ww++)
						{
							// Generate next line
							word_t** bla = opc + ca * nobs;
							word_t** blb = opc + cb * nobs;

							word_t* la = *(bla + ia);
							word_t* lb = *(blb + ib);

							word_t attributes = la[ww] ^ lb[ww];

							for (int8_t bit = WORD_BITS - 1; bit >= 0;
								 bit--, c_attribute++)
							{
								attribute_totals[c_attribute]
									-= BIT_CHECK(attributes, bit);
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
							 const uint32_t n_words_in_a_column,
							 word_t* covered_lines)
{
	for (uint32_t w = 0; w < n_words_in_a_column; w++)
	{
		covered_lines[w] |= best_column[w];
	}

	return OK;
}
