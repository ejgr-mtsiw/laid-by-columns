/*
 ============================================================================
 Name        : disjoint_matrix_mpi.c
 Author      : Eduardo Ribeiro
 Description : Structures and functions to manage the disjoint matrix
 ============================================================================
 */

#include "disjoint_matrix_mpi.h"

#include "types/dataset_t.h"
#include "types/dm_t.h"
#include "types/line_class_t.h"
#include "types/oknok_t.h"
#include "types/word_t.h"
#include "utils/bit.h"

#include <stdint.h>
#include <string.h>

oknok_t get_column(const dataset_t* dataset, const dm_t* dm,
				   const uint64_t index, word_t* column)
{

	uint64_t nw = dataset->n_words;
	uint64_t nc	   = dataset->n_classes;
	//uint64_t nobs  = dataset->n_observations;
	line_class_t*classes = dataset->classes;

	// Which word has the index attribute
	uint64_t index_word = index / WORD_BITS + 0 * dm->a_offset;

	// Which bit?
	uint8_t index_bit = WORD_BITS - (index % WORD_BITS) - 1;

	// Reset best_column
	memset(column, 0, dm->n_words_in_a_column * sizeof(word_t));

	uint64_t cs = 0;
	uint64_t ca = 0;
	uint64_t ia = 0;
	uint64_t cb = 0;
	uint64_t ib = 0;

	for (ca = 0; ca < nc - 1; ca++)
	{
		for (cb = ca + 1; cb < nc; cb++)
		{
			for (ia = 0; ia < classes[ca].n_observations; ia++)
			{
				for (ib = 0; ib < classes[cb].n_observations; ib++, cs++)
				{
					// Generate next step
					word_t* la = classes[ca].first_observation_address + ia*nw;
					word_t* lb = classes[cb].first_observation_address + ib*nw;

					word_t lxor = la[index_word] ^ lb[index_word];
					if (BIT_CHECK(lxor, index_bit))
					{
						uint64_t w = cs / WORD_BITS;
						uint64_t b = WORD_BITS - (cs % WORD_BITS) - 1;

						BIT_SET(column[w], b);
					}
				}
			}
		}
	}

	return OK;
}
