# Copyright (c) 2025 Ahmad Retha; MIT License.
#
# Bitwise encoding of a string is the process of representing each letter as a binary code and packing it onto a
# computer word. We can use bitwise operations on a computer word to solve the Lexicographically Minimal String
# Rotation problem.
#
# This is an alternative to Booth's algorithm (see: https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation),
# suited to short strings.
#
# For a text T of length n, made up of letters from a constant alphabet A, we can find the LMSR in time O(n) and
# constant extra space = O(1), provided T is short enough it can be packed into a single computer word.
#
# The limitation of this implementation is that the string T needs to be short and it should be based on a small
# alphabet (A) because the algorithm uses ceil(log_2(|A|)) bits to store each letter in a computer word (i.e. 64-bits).
# For a text T over the DNA alphabet A = (A,C,G,T), the text can be up to n=32 bases long. With a protein alphabet,
# A = (A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y), |A|=20, it can be up to n=12 residues long.
#
# Note: The implementation of get_letter_index() below represents a constant factor in the time complexity so it
# is best suited for small alphabets. You can write an alternative implementation for a large alphabet.

import argparse
import math


def get_letter_index(A, c):
    return A.find(c)


def initial_hash(text, n, A, shift_amount):
    hs = 0
    for i in range(n):
        hs = hs << shift_amount
        letter_code = get_letter_index(A, text[i])
        hs = hs | letter_code
    return hs


def rolling_hash(old_hash, mask, shift_amount, n):
    top_letter = old_hash >> (shift_amount * n - shift_amount)
    hs = old_hash << shift_amount
    hs = hs | top_letter
    hs = hs & mask
    return hs


def min_lex_rot_bw(alphabet, text):
    n = len(text)
    shift_amount = int(math.ceil(math.log2(len(alphabet))))

    if shift_amount > 64 or n > math.floor(64 / shift_amount):
        exit('Aborting... Text too long to fit into a computer word.')

    mask = (1 << (shift_amount * n)) - 1
    hs = initial_hash(text, n, alphabet, shift_amount)

    min_hs = hs
    min_i = 0

    for i in range(1, n):
        hs = rolling_hash(hs, mask, shift_amount, n)
        if hs < min_hs:
            min_hs = hs
            min_i = i

    return min_i


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate the Lexicographically Minimal String Rotation using the BW algorithm')
    parser.add_argument('-a', '--alphabet', type=str, required=True, help='The alphabet, e.g. "ACGT"')
    parser.add_argument('-t', '--text', type=str, required=True, help='The input string')
    args = parser.parse_args()

    alphabet = args.alphabet
    text = args.text

    # text = 'MISSISSIPPI'
    # alphabet = 'IMPS'
    rot = min_lex_rot_bw(alphabet, text)
    print(str(rot) + ' ' + (text*2)[rot:rot+len(text)])  # 10 (IMISSISSIPP)
