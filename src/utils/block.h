/*
 ============================================================================
 Name        : utils/block.h
 Author      : Eduardo Ribeiro
 Description : Defines expressions to calculate intervals for multiprocessing
 ============================================================================
 */

#ifndef UTILS_BLOCK_H
#define UTILS_BLOCK_H

/**
 * Number of words we can fit in a cache line
 */
#define N_WORDS_CACHE_LINE 8

// First element controlled by process id out of p processes, array length n
#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))

// Last element controlled by process id out of p processes, array length n
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id + 1), p, n) - 1)

// Size of the block controlled by process id out of p processes, array length n
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW((id + 1), p, n) - BLOCK_LOW(id, p, n))

// Process that controls item index from array with length n, p processes
#define BLOCK_OWNER(index, p, n) (((p) * ((index) + 1) - 1) / (n))

// First element controlled by process id out of p processes, array length n,
// size is multiple of m
#define BLOCK_LOW_MULTIPLE(id, p, n, m) ((id) * (n / m + (n % m != 0)) / (p) *m)

// Last element controlled by process id out of p processes, array length n
// size is multiple of m
#define BLOCK_HIGH_MULTIPLE(id, p, n, m)                                       \
	(BLOCK_LOW_MULTIPLE((id + 1), p, n, m) - 1)

// Size of the block controlled by process id out of p processes, array length n
// size is multiple of m
#define BLOCK_SIZE_MULTIPLE(id, p, n, m)                                       \
	(BLOCK_LOW_MULTIPLE((id + 1), p, n, m) - BLOCK_LOW_MULTIPLE(id, p, n, m))

#endif // UTILS_BLOCK_H
