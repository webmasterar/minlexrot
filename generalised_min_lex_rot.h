/*
* Copyright (c) 2025 Ahmad Retha; MIT License.
*/

#include <math.h>
#include <stdint.h>
#include <inttypes.h>

// Check Windows
#if _WIN32 || _WIN64
#if _WIN64
#define NUM_BITS 64
#define WORD_FMT PRIu64
typedef uint64_t WORD;
#else
#define NUM_BITS 32
#define WORD_FMT PRIu32
typedef uint32_t WORD;
#endif
#endif

// Check GCC
#if __GNUC__
#if __x86_64__ || __ppc64__
#define NUM_BITS 64
#define WORD_FMT PRIu64
typedef uint64_t WORD;
#else
#define NUM_BITS 32
#define WORD_FMT PRIu32
typedef uint32_t WORD;
#endif
#endif

// Algorithm constants
#define ALGO_BW 0
#define ALGO_KR 1
#define ROLL_DIR_BACKWARD 0
#define ROLL_DIR_FORWARD  1

// Algorithm support functions
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define LOG2(x) (log(x) / log(2.))

uint32_t min_lex_rot(char* T, uint32_t n, char* A, uint32_t a, uint32_t algo, uint32_t roll_direction);
