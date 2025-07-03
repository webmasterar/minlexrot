# Copyright (c) 2025 Ahmad Retha; MIT License.
#
# Using a fundamental concept found in the Karp-Rabin algorithm, we can encode the letters of a string as numbers and
# store them in a computer word to solve the Lexicographically Minimal String Rotation problem.
#
# This is an alternative to Booth's algorithm: https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation
# suited to short strings.
#
# For a text T of length n, made up of letters from a constant alphabet A, we can find the LMSR in time O(n) and
# constant extra space = O(1), provided T is short enough it can be represented as a number in a single computer word.
#
# The limitation of this implementation is that the string T needs to be short and it should be based on a small
# alphabet (A) because the algorithm uses log_2(|A|) bits to store each letter in a computer word (i.e. 64-bits).
# For a text T over the DNA alphabet A = (A,C,G,T), the text can be up to n=32 bases long. With a protein alphabet,
# A = (A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y), |A|=20, it can be up to n=14 residues long. Observe that this
# algorithm offers slightly greater capacity than the Bitwise algorithm which can only store up to n=12 residues.
#
# Note: The implementation of get_letter_index() below represents a constant factor in the time complexity so it
# is best suited for small alphabets. You can write an alternative implementation for a large alphabet.

import argparse
import math


def get_letter_index(A, c):
    return A.find(c)


def initial_hash(T, n, A):
    a = len(A)
    d = 0
    i = 0
    while i < n:
        ltr_idx = get_letter_index(A, T[i])
        d += a**(n-i-1) * ltr_idx
        i += 1
    return d


def rolling_hash(old_hash, top_pow, a):
    top_num = old_hash // top_pow
    d = old_hash - (top_pow * top_num)
    d *= a
    d += top_num
    return d


def min_lex_rot_kr(alphabet, text):
    n = len(text)
    
    if n > math.floor(64 / math.log(len(alphabet), 2)):
        exit('Aborting... Text too long to fit into a computer word.')
    
    hs = initial_hash(text, n, alphabet)

    a = len(alphabet)
    top_pow = a**(n-1)
    min_hs = hs
    min_i = 0

    for i in range(1, n):
        hs = rolling_hash(hs, top_pow, a)
        if hs < min_hs:
            min_hs = hs
            min_i = i

    return min_i


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate the Lexicographically Minimal String Rotation using the KR algorithm')
    parser.add_argument('-a', '--alphabet', type=str, required=True, help='The alphabet, e.g. "ACGT"')
    parser.add_argument('-t', '--text', type=str, required=True, help='The input string')
    args = parser.parse_args()

    alphabet = args.alphabet
    text = args.text

    # text = 'MISSISSIPPI'
    # alphabet = 'IMPS'
    rot = min_lex_rot_kr(alphabet, text)
    print(str(rot) + ' ' + (text*2)[rot:rot+len(text)])  # 10 (IMISSISSIPP)

