/*
* Copyright (c) 2025 Ahmad Retha; MIT License.
*
* This is the generalised C implementation of both the bitwise algorithm (min_lex_rot_bw) and the algorithm
* inspired by Karp-Rabin (min_lex_rot_kr), where there is no limit to the length of the string. It uses more computer
* words to store the hash of the string.
*
* This is an alternative to Booth's algorithm (see https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation),
* but suited to shorter strings.
*
* See the Python implementation (generalised_min_lex_rot.py) for an analysis of the algorithm.
*/

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "generalised_min_lex_rot.h"


uint32_t calculate_how_many_letters_fit_into_word(size_t a, uint32_t algo) {
	if (algo == ALGO_KR) {
		return (uint32_t) floor(NUM_BITS / LOG2(a));
    } else {
		return (uint32_t) floor(NUM_BITS / ceil(LOG2(a)));
    }
}


uint32_t calculate_how_many_words_for_text(size_t n, uint32_t letters_per_word) {
    return (uint32_t) ceil((double)n / letters_per_word);
}


WORD get_letter_idx(char letter, char* A) {
    char *e = strchr(A, letter);
    return (e == NULL) ? SIZE_MAX : (WORD) (e - A);
}


WORD* build_power_table(size_t a, uint32_t letters_per_word) {
    WORD* P = (uint64_t*) malloc(letters_per_word * sizeof(WORD));
    WORD i = 1;
    P[0] = i;
    for (; i < letters_per_word; i++) {
        P[i] = (WORD) pow((double)a, (double)i);
    }
    return P;
}


void copy_hash(WORD* H_dest, WORD* H_src, uint32_t m) {
    uint32_t i;
    for (i = 0; i < m; i++) {
        H_dest[i] = H_src[i];
    }
}


int compare_hashes(WORD* H_a, WORD* H_b, uint32_t m) {
    uint32_t i;
    for (i = 0; i < m; i++) {
        if (H_a[i] < H_b[i]) {
            return -1;
        } else if (H_a[i] > H_b[i]) {
            return 1;
        }
    }
    return 0;
}


WORD* build_initial_hash(char* T, size_t n, uint32_t m, char* A, size_t a, uint32_t letters_per_word, uint32_t algo) {
    WORD* H = (WORD*) calloc(m, sizeof(WORD));  //m = num_words
    uint32_t i, word_idx = 0;
    WORD ltr_idx;

    if (algo == ALGO_BW) {
        uint32_t shift_amt = (uint32_t) ceil(LOG2(a));
        int32_t pos_in_word = ((n-1) % letters_per_word) * shift_amt;
        for (i = 0; i < n; i++) {
            ltr_idx = get_letter_idx(T[i], A);
            H[word_idx] = H[word_idx] | (ltr_idx << pos_in_word);
            pos_in_word = pos_in_word - shift_amt;
            if (pos_in_word < 0) {
                pos_in_word = (letters_per_word - 1) * shift_amt;
                word_idx++;
            }
        }
    } else {
        WORD* P = build_power_table(a, letters_per_word);
        int32_t j = (n-1) % letters_per_word;
        for (i = 0; i < n; i++) {
            ltr_idx = get_letter_idx(T[i], A);
            H[word_idx] = H[word_idx] + P[j] * ltr_idx;
            j--;
            if (j < 0) {
                j = letters_per_word - 1;
                word_idx++;
            }
        }
        free(P);
    }

    return H;
}


// ACDE EACD DEAC CDEA forward (right to left)
void roll_hash_forward(WORD* H, uint32_t m, uint32_t n, uint32_t a, uint32_t letters_per_word, uint32_t algo) {
    uint32_t word_idx = m - 1;
    WORD prev_letter, last_letter;

    if (algo == ALGO_BW) {
        uint32_t shift_amt = (uint32_t) ceil(LOG2(a));
        WORD mask = ((WORD)1 << shift_amt) - 1;
        last_letter = H[word_idx] & mask;
        uint32_t pos_in_word = (letters_per_word-1) * shift_amt;
        while (word_idx > 0) {
            prev_letter = H[word_idx-1] & mask;
            H[word_idx] = H[word_idx] >> shift_amt;
            H[word_idx] = H[word_idx] | (prev_letter << pos_in_word);
            word_idx--;
        }
        H[0] = H[0] >> shift_amt;
        pos_in_word = ((n-1) % letters_per_word) * shift_amt;
        H[0] = H[0] | (last_letter << pos_in_word);
    } else {
        WORD top_pow = (WORD) pow((double)a, (double)(letters_per_word - 1));
        last_letter = H[word_idx] % a;
        while (word_idx > 0) {
            prev_letter = H[word_idx-1] % a;
            H[word_idx] = H[word_idx] / a;
            H[word_idx] = H[word_idx] + top_pow * prev_letter;
            word_idx--;
        }
        H[0] = H[0] / a;
        top_pow = (WORD) pow((double)a, (double)((n-1) % letters_per_word));
        H[0] = H[0] + top_pow * last_letter;
    }
}


// ACDE CDEA DEAC EACD backward (left to right), [ACDEACD]
void roll_hash_backward(WORD* H, uint32_t m, uint32_t n, uint32_t a, uint32_t letters_per_word, uint32_t algo) {
    uint32_t word_idx = 1;

    if (algo == ALGO_BW) {
        uint32_t shift_amt = (uint32_t) ceil(LOG2(a));
        uint32_t pos_in_word = ((n-1) % letters_per_word) * shift_amt;
        WORD one = 1;
        uint32_t min_shift = MIN(NUM_BITS, pos_in_word + shift_amt);
        WORD mask = (min_shift == NUM_BITS) ? ~((WORD)0) : (one <<  min_shift) - one;
        WORD first_letter = H[0] >> pos_in_word;
        H[0] = H[0] << shift_amt;
        H[0] = H[0] & mask;
        pos_in_word = (letters_per_word-1) * shift_amt;
        min_shift = MIN(NUM_BITS, pos_in_word + shift_amt);
        mask = (min_shift == NUM_BITS) ? ~((WORD)0) : (one << min_shift) - one;
        while (word_idx < m) {
            WORD curr_top = H[word_idx] >> pos_in_word;
            H[word_idx-1] = H[word_idx-1] | curr_top;
            H[word_idx] = H[word_idx] << shift_amt;
            H[word_idx] = H[word_idx] & mask;
            word_idx++;
        }
        H[m-1] = H[m-1] | first_letter;
    } else {
        uint32_t num_letters_in_first_word = (n-1) % letters_per_word;
        WORD top_pow = (WORD) pow((double)a, (double)num_letters_in_first_word);
        WORD first_letter = H[0] / top_pow;
        H[0] = H[0] - (first_letter * top_pow);
        top_pow = (WORD) pow((double)a, (double)(letters_per_word - 1));
        while (word_idx < m) {
            WORD curr_top = H[word_idx] / top_pow;
            H[word_idx-1] = H[word_idx-1] * a + curr_top;
            H[word_idx] = H[word_idx] - (curr_top * top_pow);
            word_idx++;
        }
        H[m-1] = H[m-1] * a + first_letter;
    }
}


uint32_t min_lex_rot(char* T, uint32_t n, char* A, uint32_t a, uint32_t algo, uint32_t roll_direction) {
    uint32_t letters_per_word = calculate_how_many_letters_fit_into_word(a, algo);
    uint32_t m = calculate_how_many_words_for_text(n, letters_per_word);

    WORD* H = build_initial_hash(T, n, m, A, a, letters_per_word, algo);  // Hash of T[i=0]

    WORD* M = (WORD*) malloc(m * sizeof(WORD));
    copy_hash(M, H, m);
    uint32_t M_idx = 0;

    void (*roll_words_fn)(WORD*, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t);
    roll_words_fn = (roll_direction == ROLL_DIR_BACKWARD) ? &roll_hash_backward : &roll_hash_forward;

    uint32_t i;
    for (i = 1; i < n; i++) {
        roll_words_fn(H, m, n, a, letters_per_word, algo);
        if (compare_hashes(H, M, m) == -1) {
            copy_hash(M, H, m);
            M_idx = i;
        }
    }

    free(H);
    free(M);

    return (roll_direction == ROLL_DIR_BACKWARD) ? M_idx : (n - M_idx) % n;
}
