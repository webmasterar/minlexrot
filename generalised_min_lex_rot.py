# Copyright (c) 2025 Ahmad Retha; MIT License.
#
# This is the generalised Python implementation of both the bitwise algorithm (min_lex_rot_bw) and the algorithm
# inspired by Karp-Rabin (min_lex_rot_kr), where there is no limit to the length of the string. It uses more computer
# words to store the hash of the string.
#
# This is an alternative to Booth's algorithm (see https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation),
# but suited to shorter strings.
#
# Analysis of the algorithm:
#
# Let T be a string of length n = |T| letters from a constant alphabet A of size a = |A|. Let w be the size of a
# computer word, which is w=64 bits on most modern computers.
# 
# We can determine how many letters fit into a computer word, g, using these formulas:
#   - For BW, g = floor(w / ceil(log2(a)));
#   - For KR, g = floor(w / log2(a)).
#
# Where w is divisible by the size of the alphabet without a remainder, the number of letters that fit into a computer
# word is the same, otherwise the BW algorithm uses more space. To calculate how many words are needed, i.e. the space
# complexity of the algorithm, we simply divide the number of letters in T by the number of letters that fit into a
# computer word:
#   - ceil(n / g).
#
# For each character of T, we update the ceil(n / g) word(s) encoding the text and check that it is the Lexicographically
# Minimal String Rotation (LMSR). We update the word hashes using operations that require constant time (O(1)) per word,
# either using bitwise operations (BW) or mathematical operations (KR). Afterwards, to determine the LMSR, we compare
# O(ceil(n / g)) words against the current LMSR hash.
#
# This means that for n letters, we perform O(n * ceil(n/g)) operations, giving us a worst case time complexity of O(n^2 / g)
# and a space complexity of O(ceil(n / g)).
#
# A quadratic factor in the time complexity is not an attractive attribute for any algorithm, however, when ceil(n / g) = 1,
# the algorithm runs in linear time O(n) and uses one computer word of space (8 bytes on a 64 bit machine).
#
# The size of g is determined by w and a, the size of a computer word and the size of the alphabet. The algorithm uses
# up to ceil(log_2(a)) bits to store each letter in a computer word (i.e. 64-bits). For a text T over the DNA alphabet
# A = (A,C,G,T), the text can be up to n=32 bases long to fit in a single computer word. With a protein alphabet,
# A = (A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y), a = 20, it can be up to n = 12 residues long for the BW algorithm, or
# n = 14 for the KR algorithm.
#
# I recommend using the single word implementations (provided in Python) of the algorithms in this repository for
# especially short strings:
#   - min_lex_rot_bw.py
#   - min_lex_rot_kr.py
#
# But for longer strings, this generalised implementation (provided in both Python and C) should be used instead:
#   - generalised_min_lex_rot.py
#   - generalised_min_lex_rot.c
#

import argparse
import math


ALGO_BW = 0
ALGO_KR = 1
ROLL_DIR_BACKWARD = 'BW'
ROLL_DIR_FORWARD  = 'FW'
SIZEOF_WORD = 64


def calculate_how_many_letters_fit_into_word(A, algo) -> int:
    len_A = len(A)

    if algo == ALGO_KR:
        return int( math.floor( SIZEOF_WORD / math.log(len_A, 2) ) )
    else:
        return int( math.floor( SIZEOF_WORD / math.ceil( math.log(len_A, 2) ) ) )


def calculate_how_many_words_for_text(T, letters_per_word) -> int:
    return int(math.ceil( len(T) / letters_per_word ))


def get_letter_idx(letter, A) -> int:
    return A.index(letter)


def build_power_table(a, letters_per_word):
    P = [1] * letters_per_word
    for i in range(1, letters_per_word):
        P[i] = a**i
    return P


def build_initial_hash(T, A, letters_per_word, algo) -> list[int]:
    word_idx = 0
    n = len(T)
    a = len(A)
    num_words = calculate_how_many_words_for_text(T, letters_per_word)
    H = [0] * num_words

    if algo == ALGO_BW:
        shift_amt = int( math.ceil( math.log(a, 2) ) )
        pos_in_word = ((n-1) % letters_per_word) * shift_amt
        for ltr in T:
            ltr_idx = get_letter_idx(ltr, A)
            H[word_idx] = H[word_idx] | (ltr_idx << pos_in_word)
            pos_in_word -= shift_amt
            if pos_in_word < 0:
                pos_in_word = (letters_per_word-1) * shift_amt
                word_idx += 1
    else:
        P = build_power_table(a, letters_per_word)
        i = (n-1) % letters_per_word
        for ltr in T:
            ltr_idx = get_letter_idx(ltr, A)
            H[word_idx] = H[word_idx] + P[i] * ltr_idx
            i -= 1
            if i < 0:
                i = letters_per_word - 1
                word_idx += 1

    return H


# ACDE EACD DEAC CDEA forward (right to left)
def roll_hash_forward(n, H, A, letters_per_word, algo) -> list[int]:
    m = len(H)
    word_idx = m - 1
    a = len(A)

    if algo == ALGO_BW:
        shift_amt = int( math.ceil( math.log(a, 2) ) )
        mask = (1 << shift_amt) - 1
        last_letter = H[word_idx] & mask
        pos_in_word = (letters_per_word-1) * shift_amt
        while word_idx > 0:
            prev_bot = H[word_idx-1] & mask
            H[word_idx] = H[word_idx] >> shift_amt
            H[word_idx] = H[word_idx] | (prev_bot << pos_in_word)
            word_idx -= 1
        H[0] = H[0] >> shift_amt
        pos_in_word = ((n-1) % letters_per_word) * shift_amt
        H[0] = H[0] | (last_letter << pos_in_word)
    else:
        top_pow = a**(letters_per_word - 1)
        last_letter = H[word_idx] % a
        while word_idx > 0:
            prev_bot = H[word_idx-1] % a
            H[word_idx] = H[word_idx] // a
            H[word_idx] = H[word_idx] + top_pow * prev_bot
            word_idx -= 1
        H[0] = H[0] // a
        top_pow = a**((n-1) % letters_per_word)
        H[0] = H[0] + top_pow * last_letter

    return H


# ACDE CDEA DEAC EACD backward (left to right), [ACDEACD]
def roll_hash_backward(n, H, A, letters_per_word, algo) -> list[int]:
    m = len(H)
    word_idx = 1
    a = len(A)

    if algo == ALGO_BW:
        shift_amt = int( math.ceil( math.log(a, 2) ) )
        pos_in_word = ((n-1) % letters_per_word) * shift_amt
        mask = (1 << min(SIZEOF_WORD, pos_in_word + shift_amt)) - 1
        first_letter = H[0] >> pos_in_word
        H[0] = H[0] << shift_amt
        H[0] = H[0] & mask
        pos_in_word = (letters_per_word-1) * shift_amt
        mask = (1 << min(SIZEOF_WORD, pos_in_word + shift_amt)) - 1
        while word_idx < m:
            curr_top = H[word_idx] >> pos_in_word
            H[word_idx-1] = H[word_idx-1] | curr_top
            H[word_idx] = H[word_idx] << shift_amt
            H[word_idx] = H[word_idx] & mask
            word_idx += 1
        H[m-1] = H[m-1] | first_letter
    else:
        num_letters_in_first_word = (n-1) % letters_per_word
        top_pow = a**num_letters_in_first_word
        first_letter = H[0] // top_pow
        H[0] = H[0] - (first_letter * top_pow)
        top_pow = a**(letters_per_word - 1)
        while word_idx < m:
            curr_top = H[word_idx] // top_pow
            H[word_idx-1] = H[word_idx-1] * a + curr_top
            H[word_idx] = H[word_idx] - (curr_top * top_pow)
            word_idx += 1
        H[m-1] = H[m-1] * a + first_letter

    return H


def compare_hashes(H_a, H_b) -> int:
    for i in range(len(H_a)):
        if H_a[i] < H_b[i]:
            return -1
        elif H_a[i] > H_b[i]:
            return 1
    return 0


def copy_hash(H_dest, H_src):
    for i in range(len(H_src)):
        H_dest[i] = H_src[i]


def min_lex_rot(T, A, algo, roll_direction) -> int:
    n = len(T)
    letters_per_word = calculate_how_many_letters_fit_into_word(A, algo)
    H = build_initial_hash(T, A, letters_per_word, algo)  # Hash of T[i=0]
    M = H[:]
    M_idx = 0
    if roll_direction == ROLL_DIR_BACKWARD:
        roll_method = roll_hash_backward
    else:
        roll_method = roll_hash_forward
    
    for i in range(1, n):
        H = roll_method(n, H, A, letters_per_word, algo)
        if compare_hashes(H, M) == -1:
            copy_hash(M, H)
            M_idx = i

    if roll_direction == ROLL_DIR_BACKWARD:
        return M_idx
    return (n - M_idx) % n


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A generalised implementation of the Lexicographically Minimal String Rotation using algorithm BW or KR')
    parser.add_argument('-a', '--alphabet', type=str, required=True, help='The alphabet, e.g. "ACGT"')
    parser.add_argument('-g', '--algorithm', type=str, required=True, choices=['BW', 'KR'], help='The algorithm to use (BW=Bitwise or KR=Karp-Rabin).')
    parser.add_argument('-d', '--direction',type=str, required=True, default='BW', choices=['FW', 'BW'], help='The direction of rotation (FW=Forward or BW=Backward).')
    parser.add_argument('-t', '--text', type=str, required=True, help='The input string')
    args = parser.parse_args()

    alphabet = args.alphabet
    algo = args.algorithm
    direction = args.direction
    text = args.text

    rot = min_lex_rot(text, alphabet, algo, direction)
    n = len(text)
    print(str(rot) + ' ' + (text*2)[rot:rot+n])
