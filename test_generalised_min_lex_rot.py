# Copyright (c) 2025 Ahmad Retha; MIT License.

import pytest

from generalised_min_lex_rot import *
from random import choices, randint


DNA_ALPHABET     = ['A', 'C', 'G', 'T']
PROTEIN_ALPHABET = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


@pytest.mark.parametrize('algo,A,expected', [
	(ALGO_BW, DNA_ALPHABET,     32),
	(ALGO_KR, DNA_ALPHABET,     32),
	(ALGO_BW, PROTEIN_ALPHABET, 12),
	(ALGO_KR, PROTEIN_ALPHABET, 14),
]) 
def test_calculate_how_many_letters_fit_into_word(algo, A, expected):
	num_letters = calculate_how_many_letters_fit_into_word(A, algo)
	assert num_letters == expected
	#log 20 / log 2 = 4.32
	#64 / 4.32 = 14.8 = floor(14.8) = 14 residues per word using kr filling
	#64 / ceil(4.32) = 64/5 = 12.8 = floor(12.8) = 12 residues per word using bw filling


@pytest.mark.parametrize('n,A,algo,expected', [
	(1,  DNA_ALPHABET,     ALGO_BW, 1),
	(32, DNA_ALPHABET,     ALGO_BW, 1),
	(33, DNA_ALPHABET,     ALGO_BW, 2),
	(1,  DNA_ALPHABET,     ALGO_KR, 1),
	(32, DNA_ALPHABET,     ALGO_KR, 1),
	(33, DNA_ALPHABET,     ALGO_KR, 2),
	(12, PROTEIN_ALPHABET, ALGO_BW, 1),
	(13, PROTEIN_ALPHABET, ALGO_BW, 2),
	(14, PROTEIN_ALPHABET, ALGO_KR, 1),
	(15, PROTEIN_ALPHABET, ALGO_KR, 2),
])
def test_calculate_how_many_words_for_text(n, A, algo, expected):
	T = ['A'] * n
	letters_per_word = calculate_how_many_letters_fit_into_word(A, algo)
	assert calculate_how_many_words_for_text(T, letters_per_word) == expected


@pytest.mark.parametrize('A', [DNA_ALPHABET, PROTEIN_ALPHABET])
def test_get_letter_idx(A):
	for i, a in enumerate(A):
		assert get_letter_idx(a, A) == i


@pytest.mark.parametrize('algo,T,A,expected', [
	(ALGO_BW, 'C', DNA_ALPHABET, [0b1]),
	(ALGO_BW, 'CA', DNA_ALPHABET, [0b100]),
	(ALGO_BW, 'ACD', PROTEIN_ALPHABET, [0b100010]),
	(ALGO_BW, 'ACDE', PROTEIN_ALPHABET, [0b10001000011]),
	(ALGO_BW, 'ACDEF', PROTEIN_ALPHABET, [0b1000100001100100]),
	(ALGO_BW, 'ACDEFG', PROTEIN_ALPHABET, [0b100010000110010000101]),
	(ALGO_BW, 'ACDEFGH', PROTEIN_ALPHABET, [0b10001000011001000010100110]),
    (ALGO_BW, 'ACDEFGHI', PROTEIN_ALPHABET, [0b1000100001100100001010011000111]),
    (ALGO_BW, 'ACDEFGHIK', PROTEIN_ALPHABET, [0b100010000110010000101001100011101000]),
	(ALGO_BW, PROTEIN_ALPHABET[0:12], PROTEIN_ALPHABET, [0b100010000110010000101001100011101000010010101001011]),
	(ALGO_BW, PROTEIN_ALPHABET[0:13], PROTEIN_ALPHABET, [0, 0b10001000011001000010100110001110100001001010100101101100]),
	(ALGO_BW, PROTEIN_ALPHABET[0:12] * 2, PROTEIN_ALPHABET, [0b100010000110010000101001100011101000010010101001011, 0b100010000110010000101001100011101000010010101001011]),
    (ALGO_KR, 'C', DNA_ALPHABET, [1]),
    (ALGO_KR, 'CA', DNA_ALPHABET, [4]),
    (ALGO_KR, 'GCA', DNA_ALPHABET, [36]),
    (ALGO_KR, 'TGCA', DNA_ALPHABET, [228]),
    (ALGO_KR, 'A' * 32, DNA_ALPHABET, [0]),
    (ALGO_KR, 'Y' * 14, PROTEIN_ALPHABET, [1638399999999999999]),
    (ALGO_KR, 'Y' * 15, PROTEIN_ALPHABET, [19, 1638399999999999999]),
    (ALGO_KR, 'Y' * 16, PROTEIN_ALPHABET, [399, 1638399999999999999]),
])
def test_build_initial_hash(algo, T, A, expected: list[int]):
	letters_per_word = calculate_how_many_letters_fit_into_word(A, algo)
	H = build_initial_hash(T, A, letters_per_word, algo)
	assert H == expected


@pytest.mark.parametrize('algo,T,A,expected', [
    (ALGO_BW, 'ACDE', PROTEIN_ALPHABET, [0b11000000000100010]),
    (ALGO_BW, 'EACD', PROTEIN_ALPHABET, [0b10000110000000001]),
    (ALGO_BW, 'DEAC', PROTEIN_ALPHABET, [0b1000100001100000]),
    (ALGO_BW, 'CDEA', PROTEIN_ALPHABET, [0b10001000011]),
    (ALGO_BW, 'ACDEACDEACDEA', PROTEIN_ALPHABET, [0, 0b100010000110000000001000100001100000000010001000011]),
    (ALGO_BW, 'AACDEACDEACDE', PROTEIN_ALPHABET, [0b11, 0b1000100001100000000010001000011000000000100010]),
    (ALGO_BW, 'EAACDEACDEACD', PROTEIN_ALPHABET, [0b10, 0b000110000000000000010001000011000000000100010000110000000001]),
    (ALGO_BW, 'DEAACDEACDEAC', PROTEIN_ALPHABET, [0b1, 0b100001100000000000000100010000110000000001000100001100000]),
    (ALGO_BW, 'DDEAACDEACDEAC', PROTEIN_ALPHABET, [0b100010, 0b100001100000000000000100010000110000000001000100001100000]),
    (ALGO_KR, 'Y', PROTEIN_ALPHABET, [19]),
    (ALGO_KR, 'YA', PROTEIN_ALPHABET, [19]),
    (ALGO_KR, 'Y' * 13 + 'D', PROTEIN_ALPHABET, [245759999999999999]),
    (ALGO_KR, 'Y' * 14, PROTEIN_ALPHABET, [1638399999999999999]),
    (ALGO_KR, 'Y' * 14 + 'D', PROTEIN_ALPHABET, [2, 1638399999999999999]),
    (ALGO_KR, 'Y' * 15, PROTEIN_ALPHABET, [19, 1638399999999999999]),
    (ALGO_KR, 'A' + 'Y' * 14 + 'D', PROTEIN_ALPHABET, [40, 1638399999999999999]),
])
def test_roll_hash_forward(algo, T, A, expected: list[int]):
    letters_per_word = calculate_how_many_letters_fit_into_word(A, algo)
    prev_H = build_initial_hash(T, A, letters_per_word, algo)
    n = len(T)
    new_H = roll_hash_forward(n, prev_H, A, letters_per_word, algo)
    assert new_H == expected


@pytest.mark.parametrize('algo,T,A,expected', [
    (ALGO_BW, 'ACDE', PROTEIN_ALPHABET, [0b1000100001100000]),
    (ALGO_BW, 'CDEA', PROTEIN_ALPHABET, [0b10000110000000001]),
    (ALGO_BW, 'DEAC', PROTEIN_ALPHABET, [0b11000000000100010]),
    (ALGO_BW, 'EACD', PROTEIN_ALPHABET, [0b10001000011]),
    (ALGO_BW, 'ACDEACDEACDEA', PROTEIN_ALPHABET, [0b1, 0b100001100000000010001000011000000000100010000110000000000]),
    (ALGO_BW, 'DDEAACDEACDEAC', PROTEIN_ALPHABET, [0b1000011, 0b1000100001100000000010001000011000000000100010]),
    (ALGO_KR, 'Y', PROTEIN_ALPHABET, [19]),
    (ALGO_KR, 'YA', PROTEIN_ALPHABET, [19]),
    (ALGO_KR, 'DY' + 'A' * 12, PROTEIN_ALPHABET, [1556480000000000002]),
    (ALGO_KR, 'Y' * 14, PROTEIN_ALPHABET, [1638399999999999999]),
    (ALGO_KR, 'D' + 'Y' * 14, PROTEIN_ALPHABET, [19, 1638399999999999982]),
    (ALGO_KR, 'Y' * 15, PROTEIN_ALPHABET, [19, 1638399999999999999]),
    (ALGO_KR, 'DA' + 'Y' * 14, PROTEIN_ALPHABET, [19, 1638399999999999982]),
])
def test_roll_hash_backward(algo, T, A, expected: list[int]):
    letters_per_word = calculate_how_many_letters_fit_into_word(A, algo)
    prev_H = build_initial_hash(T, A, letters_per_word, algo)
    n = len(T)
    new_H = roll_hash_backward(n, prev_H, A, letters_per_word, algo)
    assert new_H == expected


@pytest.mark.parametrize('a,b,expected', [
    ([0], [1], -1),
    ([1], [0],  1),
    ([0], [0],  0),
    ([1], [1],  0),
    ([0, 1], [1, 0], -1),
    ([1, 0], [0, 1],  1),
    ([0, 0], [0, 0],  0),
    ([1, 1], [1, 1],  0),
])
def test_compare_hashes(a: list[int], b: list[int], expected: int):
    assert compare_hashes(a, b) == expected


@pytest.mark.parametrize('a,letters_per_word', [
    (4,  32),
    (20, 14)
])
def test_build_power_table(a, letters_per_word):
    P = build_power_table(a, letters_per_word)
    assert P[letters_per_word-1] == a**(letters_per_word-1)


@pytest.mark.parametrize('algo,direction,A,T,expected', [
    (ALGO_BW, ROLL_DIR_BACKWARD, DNA_ALPHABET, 'AAA', 0),
    (ALGO_BW, ROLL_DIR_BACKWARD, DNA_ALPHABET, 'ATA', 2),
    (ALGO_BW, ROLL_DIR_BACKWARD, DNA_ALPHABET, 'AT'*33, 0),
    (ALGO_BW, ROLL_DIR_BACKWARD, DNA_ALPHABET, 'AAA'*64, 0),
    (ALGO_BW, ROLL_DIR_BACKWARD, PROTEIN_ALPHABET, 'AAAACCCCDDDDEEEE', 0),
    (ALGO_BW, ROLL_DIR_BACKWARD, PROTEIN_ALPHABET, 'ACCCCDDDDEEEEAAA', 13),
    (ALGO_BW, ROLL_DIR_BACKWARD, PROTEIN_ALPHABET, 'ACCCCAAAAEEEEAA', 5),
    (ALGO_BW, ROLL_DIR_FORWARD, DNA_ALPHABET, 'AAA', 0),
    (ALGO_BW, ROLL_DIR_FORWARD, DNA_ALPHABET, 'ATA', 2),
    (ALGO_BW, ROLL_DIR_FORWARD, PROTEIN_ALPHABET, 'AAAACCCCDDDDEEEE', 0),
    (ALGO_BW, ROLL_DIR_FORWARD, PROTEIN_ALPHABET, 'ACCCCDDDDEEEEAAA', 13),
    (ALGO_BW, ROLL_DIR_FORWARD, PROTEIN_ALPHABET, 'ACCCCAAAAEEEEAA', 5),
    (ALGO_KR, ROLL_DIR_FORWARD, DNA_ALPHABET, 'AAA', 0),
    (ALGO_KR, ROLL_DIR_FORWARD, DNA_ALPHABET, 'ATA', 2),
    (ALGO_KR, ROLL_DIR_FORWARD, DNA_ALPHABET, 'AT'*33, 0),
    (ALGO_KR, ROLL_DIR_FORWARD, DNA_ALPHABET, 'AAA'*64, 0),
    (ALGO_KR, ROLL_DIR_FORWARD, PROTEIN_ALPHABET, 'AAAACCCCDDDDEEEE', 0),
    (ALGO_KR, ROLL_DIR_FORWARD, PROTEIN_ALPHABET, 'ACCCCDDDDEEEEAAA', 13),
    (ALGO_KR, ROLL_DIR_FORWARD, PROTEIN_ALPHABET, 'ACCCCAAAAEEEEAA', 5),
    (ALGO_KR, ROLL_DIR_BACKWARD, DNA_ALPHABET, 'AAA', 0),
    (ALGO_KR, ROLL_DIR_BACKWARD, DNA_ALPHABET, 'ATA', 2),
    (ALGO_KR, ROLL_DIR_BACKWARD, DNA_ALPHABET, 'AT'*33, 0),
    (ALGO_KR, ROLL_DIR_BACKWARD, DNA_ALPHABET, 'AAA'*64, 0),
    (ALGO_KR, ROLL_DIR_BACKWARD, PROTEIN_ALPHABET, 'AAAACCCCDDDDEEEE', 0),
    (ALGO_KR, ROLL_DIR_BACKWARD, PROTEIN_ALPHABET, 'ACCCCDDDDEEEEAAA', 13),
    (ALGO_KR, ROLL_DIR_BACKWARD, PROTEIN_ALPHABET, 'ACCCCAAAAEEEEAA', 5),
])
def test_min_lex_rot(algo, direction, A, T, expected):
    n = len(T)
    TT = T + T
    i = min_lex_rot(T, A, algo, direction)
    assert TT[i:i+n] == TT[expected:expected+n]


def test_min_lex_rot_randstrings():
    for algo in [ALGO_BW, ALGO_KR]:
        for direction in [ROLL_DIR_BACKWARD, ROLL_DIR_FORWARD]:
            for A in [DNA_ALPHABET, PROTEIN_ALPHABET]:
                for _ in range(100):
                    k = randint(500, 1500)
                    T_rand = choices(A, k=k)
                    TT = T_rand + T_rand
                    rot = min_lex_rot(T_rand, A, algo, direction)
                    T_mlr = TT[rot:rot+k]
                    for i in range(k):
                        assert TT[i:i+k] >= T_mlr
