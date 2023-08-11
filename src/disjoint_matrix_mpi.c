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
#include "types/oknok_t.h"
#include "types/word_t.h"
#include "utils/bit.h"

#include <stdint.h>
#include <string.h>

oknok_t get_column(const dataset_t* dataset, const dm_t* dm,
				   const uint64_t index, word_t* column)
{

	uint64_t nc	   = dataset->n_classes;
	uint64_t nobs  = dataset->n_observations;
	word_t** opc   = dataset->observations_per_class;
	uint64_t* nopc = dataset->n_observations_per_class;

	// Which word has the index attribute
	uint64_t index_word = index / WORD_BITS + 0 * dm->a_offset;

	// Which bit?
	uint8_t index_bit = WORD_BITS - (index % WORD_BITS) - 1;

	// Reset best_column
	memset(column, 0, dm->n_words_in_a_column * sizeof(word_t));

	// TODO: I think we can optimize the order of the lines
	// to optimize cache usage and get faster results
	// when calculating the attributes totals
	uint64_t cs = 0;
	uint64_t ca = 0;
	uint64_t ia = 0;
	uint64_t cb = 0;
	uint64_t ib = 0;

	// TODO: is there a better way?

	for (ca = 0; ca < nc - 1; ca++)
	{
		for (cb = ca + 1; cb < nc; cb++)
		{
			for (ia = 0; ia < nopc[ca]; ia++)
			{
				for (ib = 0; ib < nopc[cb]; ib++, cs++)
				{
					// Generate next step
					word_t** bla = opc + ca * nobs;
					word_t** blb = opc + cb * nobs;

					word_t* la = *(bla + ia);
					word_t* lb = *(blb + ib);

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
