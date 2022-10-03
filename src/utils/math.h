/*
 ============================================================================
 Name        : utils/math.h
 Author      : Eduardo Ribeiro
 Description : Mathematical helpers
 ============================================================================
 */

#ifndef UTILS_MATH_H
#define UTILS_MATH_H

#include <stdint.h>

#define MAX(a, b) ((a) > (b) ? a : b)
#define MIN(a, b) ((a) < (b) ? a : b)

uint32_t roundUp(uint32_t numToRound, uint32_t multiple);

#endif // UTILS_MATH_H
