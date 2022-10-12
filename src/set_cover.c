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
#include "utils/bit.h"

#include <stdint.h>

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

	best_attribute_t best = { .attribute = max_attribute, .max = max_total };
	return best;
}

//void calculate_initial_totals(dm_t* dm, dataset_t* dataset, uint32_t* totals)
//{
//	for (uint32_t line = 0; line < dm->n_matrix_lines; line++)
//	{
//		word_t* la = dataset->data + dm->steps[line].indexA * dataset->n_words
//			+ dm->a_offset;
//		word_t* lb = dataset->data + dm->steps[line].indexB * dataset->n_words
//			+ dm->a_offset;
//
//		uint32_t c_attribute = 0;
//
//		for (uint32_t w = 0; w < dm->a_size; w++)
//		{
//			word_t lxor = la[w] ^ lb[w];
//
//			add_to_totals(totals, &c_attribute, lxor);
//		}
//	}
//}
//
//#ifdef __GNUC__
///*
// * GCC has trouble optimizing this loop se we unrolled it manually
// * CLANG can do it faster using AVX
// */
//void add_to_totals(uint32_t* totals, uint32_t* start, const word_t lxor)
//{
//	totals[(*start) + 0] += !!(lxor & (1UL << 63));
//	totals[(*start) + 1] += !!(lxor & (1UL << 62));
//	totals[(*start) + 2] += !!(lxor & (1UL << 61));
//	totals[(*start) + 3] += !!(lxor & (1UL << 60));
//	totals[(*start) + 4] += !!(lxor & (1UL << 59));
//	totals[(*start) + 5] += !!(lxor & (1UL << 58));
//	totals[(*start) + 6] += !!(lxor & (1UL << 57));
//	totals[(*start) + 7] += !!(lxor & (1UL << 56));
//	totals[(*start) + 8] += !!(lxor & (1UL << 55));
//	totals[(*start) + 9] += !!(lxor & (1UL << 54));
//
//	totals[(*start) + 10] += !!(lxor & (1UL << 53));
//	totals[(*start) + 11] += !!(lxor & (1UL << 52));
//	totals[(*start) + 12] += !!(lxor & (1UL << 51));
//	totals[(*start) + 13] += !!(lxor & (1UL << 50));
//	totals[(*start) + 14] += !!(lxor & (1UL << 49));
//	totals[(*start) + 15] += !!(lxor & (1UL << 48));
//	totals[(*start) + 16] += !!(lxor & (1UL << 47));
//	totals[(*start) + 17] += !!(lxor & (1UL << 46));
//	totals[(*start) + 18] += !!(lxor & (1UL << 45));
//	totals[(*start) + 19] += !!(lxor & (1UL << 44));
//
//	totals[(*start) + 20] += !!(lxor & (1UL << 43));
//	totals[(*start) + 21] += !!(lxor & (1UL << 42));
//	totals[(*start) + 22] += !!(lxor & (1UL << 41));
//	totals[(*start) + 23] += !!(lxor & (1UL << 40));
//	totals[(*start) + 24] += !!(lxor & (1UL << 39));
//	totals[(*start) + 25] += !!(lxor & (1UL << 38));
//	totals[(*start) + 26] += !!(lxor & (1UL << 37));
//	totals[(*start) + 27] += !!(lxor & (1UL << 36));
//	totals[(*start) + 28] += !!(lxor & (1UL << 35));
//	totals[(*start) + 29] += !!(lxor & (1UL << 34));
//
//	totals[(*start) + 30] += !!(lxor & (1UL << 33));
//	totals[(*start) + 31] += !!(lxor & (1UL << 32));
//	totals[(*start) + 32] += !!(lxor & (1UL << 31));
//	totals[(*start) + 33] += !!(lxor & (1UL << 30));
//	totals[(*start) + 34] += !!(lxor & (1UL << 29));
//	totals[(*start) + 35] += !!(lxor & (1UL << 28));
//	totals[(*start) + 36] += !!(lxor & (1UL << 27));
//	totals[(*start) + 37] += !!(lxor & (1UL << 26));
//	totals[(*start) + 38] += !!(lxor & (1UL << 25));
//	totals[(*start) + 39] += !!(lxor & (1UL << 24));
//
//	totals[(*start) + 40] += !!(lxor & (1UL << 23));
//	totals[(*start) + 41] += !!(lxor & (1UL << 22));
//	totals[(*start) + 42] += !!(lxor & (1UL << 21));
//	totals[(*start) + 43] += !!(lxor & (1UL << 20));
//	totals[(*start) + 44] += !!(lxor & (1UL << 19));
//	totals[(*start) + 45] += !!(lxor & (1UL << 18));
//	totals[(*start) + 46] += !!(lxor & (1UL << 17));
//	totals[(*start) + 47] += !!(lxor & (1UL << 16));
//	totals[(*start) + 48] += !!(lxor & (1UL << 15));
//	totals[(*start) + 49] += !!(lxor & (1UL << 14));
//
//	totals[(*start) + 50] += !!(lxor & (1UL << 13));
//	totals[(*start) + 51] += !!(lxor & (1UL << 12));
//	totals[(*start) + 52] += !!(lxor & (1UL << 11));
//	totals[(*start) + 53] += !!(lxor & (1UL << 10));
//	totals[(*start) + 54] += !!(lxor & (1UL << 9));
//	totals[(*start) + 55] += !!(lxor & (1UL << 8));
//	totals[(*start) + 56] += !!(lxor & (1UL << 7));
//	totals[(*start) + 57] += !!(lxor & (1UL << 6));
//	totals[(*start) + 58] += !!(lxor & (1UL << 5));
//	totals[(*start) + 59] += !!(lxor & (1UL << 4));
//
//	totals[(*start) + 60] += !!(lxor & (1UL << 3));
//	totals[(*start) + 61] += !!(lxor & (1UL << 2));
//	totals[(*start) + 62] += !!(lxor & (1UL << 1));
//	totals[(*start) + 63] += !!(lxor & (1UL << 0));
//
//	(*start) += 64;
//}
//#else // __GNUC__
//void add_to_totals(uint32_t* totals, uint32_t* start, const word_t lxor)
//{
//	for (int8_t bit = WORD_BITS - 1; bit >= 0; bit--, (*start)++)
//	{
//		totals[*start] += BIT_CHECK(lxor, bit);
//	}
//}
//#endif // __GNUC__
//
//void update_totals(dm_t* dm, dataset_t* dataset, uint32_t* totals,
//				   word_t* covered_lines, word_t* best_column)
//{
//	for (uint32_t line = 0; line < dm->n_matrix_lines; line++)
//	{
//		// Is this line covered?
//		// Yes: skip
//		// No. Is it covered by the best attribute?
//		// Yes: add
//		// No: skip
//		uint32_t current_word = line / WORD_BITS;
//		uint8_t current_bit	  = WORD_BITS - line % WORD_BITS - 1;
//
//		// Is this line already covered?
//		if (BIT_CHECK(covered_lines[current_word], current_bit))
//		{
//			// This line is already covered: skip
//			continue;
//		}
//
//		// Is this line covered by the best attribute?
//		if (!BIT_CHECK(best_column[current_word], current_bit))
//		{
//			// This line is NOT covered: skip
//			continue;
//		}
//
//		// This line was uncovered, but is covered now
//		// Calculate attributes totals
//		word_t* la = dataset->data + dm->steps[line].indexA * dataset->n_words
//			+ dm->a_offset;
//		word_t* lb = dataset->data + dm->steps[line].indexB * dataset->n_words
//			+ dm->a_offset;
//
//		uint32_t c_attribute = 0;
//
//		for (uint32_t w = 0; w < dm->a_size; w++)
//		{
//			word_t lxor = la[w] ^ lb[w];
//
//			sub_from_totals(totals, &c_attribute, lxor);
//
//			//				__m256i lxor_ = _mm256_set1_epi64x(lxor);
//			//				__m256i zeros = _mm256_set1_epi64x(0);
//			//
//			//				for (int8_t bit = WORD_BITS, pos = 0; bit > 0;
//			//					 bit -= 4, pos += 4, c_attribute += 4)
//			//				{
//			//					__m256i mask = _mm256_set_epi64x(
//			//						AND_MASK_TABLE[bit - 1], AND_MASK_TABLE[bit
//			//- 2], 						AND_MASK_TABLE[bit - 3],
//			// AND_MASK_TABLE[bit - 4]);
//			//					__m256i attr_totals = _mm256_set_epi64x(
//			//						attribute_totals[pos], attribute_totals[pos
//			//+ 1], 						attribute_totals[pos + 2],
//			// attribute_totals[pos + 3]);
//			//
//			//					__m256i a	 = _mm256_and_si256(lxor_, mask);
//			//					__m256i ones = _mm256_cmpgt_epi64(a, zeros);
//			//
//			//					attr_totals = _mm256_sub_epi64(attr_totals,
//			// ones);
//			//
//			//					attribute_totals[c_attribute + 0] =
//			// attr_totals[0]; 					attribute_totals[c_attribute +
//			// 1]
//			// =
//			// attr_totals[1]; 					attribute_totals[c_attribute +
//			// 2]
//			// =
//			// attr_totals[2]; 					attribute_totals[c_attribute +
//			// 3] = attr_totals[3];
//			//				}
//		}
//	}
//}
//
//#ifdef __GNUC__
///*
// * GCC has trouble optimizing this loop se we unrolled it manually
// * CLANG can do it faster using AVX
// */
//void sub_from_totals(uint32_t* totals, uint32_t* start, const word_t lxor)
//{
//	totals[(*start) + 0] -= !!(lxor & (1UL<<63));
//	totals[(*start) + 1] -= !!(lxor & (1UL<<62));
//	totals[(*start) + 2] -= !!(lxor & (1UL<<61));
//	totals[(*start) + 3] -= !!(lxor & (1UL<<60));
//	totals[(*start) + 4] -= !!(lxor & (1UL<<59));
//	totals[(*start) + 5] -= !!(lxor & (1UL<<58));
//	totals[(*start) + 6] -= !!(lxor & (1UL<<57));
//	totals[(*start) + 7] -= !!(lxor & (1UL<<56));
//	totals[(*start) + 8] -= !!(lxor & (1UL<<55));
//	totals[(*start) + 9] -= !!(lxor & (1UL<<54));
//
//	totals[(*start) + 10] -= !!(lxor & (1UL<<53));
//	totals[(*start) + 11] -= !!(lxor & (1UL<<52));
//	totals[(*start) + 12] -= !!(lxor & (1UL<<51));
//	totals[(*start) + 13] -= !!(lxor & (1UL<<50));
//	totals[(*start) + 14] -= !!(lxor & (1UL<<49));
//	totals[(*start) + 15] -= !!(lxor & (1UL<<48));
//	totals[(*start) + 16] -= !!(lxor & (1UL<<47));
//	totals[(*start) + 17] -= !!(lxor & (1UL<<46));
//	totals[(*start) + 18] -= !!(lxor & (1UL<<45));
//	totals[(*start) + 19] -= !!(lxor & (1UL<<44));
//
//	totals[(*start) + 20] -= !!(lxor & (1UL<<43));
//	totals[(*start) + 21] -= !!(lxor & (1UL<<42));
//	totals[(*start) + 22] -= !!(lxor & (1UL<<41));
//	totals[(*start) + 23] -= !!(lxor & (1UL<<40));
//	totals[(*start) + 24] -= !!(lxor & (1UL<<39));
//	totals[(*start) + 25] -= !!(lxor & (1UL<<38));
//	totals[(*start) + 26] -= !!(lxor & (1UL<<37));
//	totals[(*start) + 27] -= !!(lxor & (1UL<<36));
//	totals[(*start) + 28] -= !!(lxor & (1UL<<35));
//	totals[(*start) + 29] -= !!(lxor & (1UL<<34));
//
//	totals[(*start) + 30] -= !!(lxor & (1UL<<33));
//	totals[(*start) + 31] -= !!(lxor & (1UL<<32));
//	totals[(*start) + 32] -= !!(lxor & (1UL<<31));
//	totals[(*start) + 33] -= !!(lxor & (1UL<<30));
//	totals[(*start) + 34] -= !!(lxor & (1UL<<29));
//	totals[(*start) + 35] -= !!(lxor & (1UL<<28));
//	totals[(*start) + 36] -= !!(lxor & (1UL<<27));
//	totals[(*start) + 37] -= !!(lxor & (1UL<<26));
//	totals[(*start) + 38] -= !!(lxor & (1UL<<25));
//	totals[(*start) + 39] -= !!(lxor & (1UL<<24));
//
//	totals[(*start) + 40] -= !!(lxor & (1UL<<23));
//	totals[(*start) + 41] -= !!(lxor & (1UL<<22));
//	totals[(*start) + 42] -= !!(lxor & (1UL<<21));
//	totals[(*start) + 43] -= !!(lxor & (1UL<<20));
//	totals[(*start) + 44] -= !!(lxor & (1UL<<19));
//	totals[(*start) + 45] -= !!(lxor & (1UL<<18));
//	totals[(*start) + 46] -= !!(lxor & (1UL<<17));
//	totals[(*start) + 47] -= !!(lxor & (1UL<<16));
//	totals[(*start) + 48] -= !!(lxor & (1UL<<15));
//	totals[(*start) + 49] -= !!(lxor & (1UL<<14));
//
//	totals[(*start) + 50] -= !!(lxor & (1UL<<13));
//	totals[(*start) + 51] -= !!(lxor & (1UL<<12));
//	totals[(*start) + 52] -= !!(lxor & (1UL<<11));
//	totals[(*start) + 53] -= !!(lxor & (1UL<<10));
//	totals[(*start) + 54] -= !!(lxor & (1UL<<9));
//	totals[(*start) + 55] -= !!(lxor & (1UL<<8));
//	totals[(*start) + 56] -= !!(lxor & (1UL<<7));
//	totals[(*start) + 57] -= !!(lxor & (1UL<<6));
//	totals[(*start) + 58] -= !!(lxor & (1UL<<5));
//	totals[(*start) + 59] -= !!(lxor & (1UL<<4));
//
//	totals[(*start) + 60] -= !!(lxor & (1UL<<3));
//	totals[(*start) + 61] -= !!(lxor & (1UL<<2));
//	totals[(*start) + 62] -= !!(lxor & (1UL<<1));
//	totals[(*start) + 63] -= !!(lxor & (1UL<<0));
//
//	(*start) += 64;
//}
//#else // __GNUC__
//void sub_from_totals(uint32_t* totals, uint32_t* start, const word_t lxor)
//{
//	for (int8_t bit = WORD_BITS - 1; bit >= 0; bit--, (*start)++)
//	{
//		totals[*start] -= BIT_CHECK(lxor, bit);
//	}
//}
//#endif // __GNUC__
