# Copyright (c) 2025 Ahmad Retha; MIT License.
#
# This is an algorithm I am calling min_lex_rot_count. It solves the Lexicographically
# Minimal String Rotation (LMSR) problem in O(n) time and O(a) extra space, where
# n is the length of a text T and a is the size of a fixed-size alphabet A. It
# works by finding the longest contiguous stretch of characters that is lexicographically
# ascending and confirms that it is lexicographically minimal in the whole string.
#
# This is an alternative to Booth's algorithm (see https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation).
# 
# How it works:
# 
# The intuition is simply that if we listed all the rotations of T and sorted
# them in ascending order -- well, the first string would begin with the longest
# contiguous stretch of characters that is both minimal and lexicographically
# ascending.
#
# So all we have to do is find the longest contiguous stretch in any rotation
# of T. We can find this by effectively concatenating T with itself, creating a
# string T' of length 2*n-1 (where n = |T|). We look at each character T'[i] from
# left to right, checking it is lexicographically equal to or greater than the
# previous character and incrementing the count of this character in a table C
# of size a = |A|. If we find that character T'[i] < T'[i-1], we end the stretch
# and check to confirm that it is minimal so far in the text, then we reset the
# current counter, begin a new stretch and continue looking through the rest of
# the string.
#
# Analysis of the algorithm:
#
# Before we scan the text we need to create a table C for the counts of the alphabet
# in the current stretch as well as another table M for the minimal stretch. The
# size of each table is a = |A| and require O(a) time to create. We keep tabs on
# the starting and ending position of the stretch using two pairs of variables
# and this requires constant space and time to update.
#
# We scan through a sequence of 2n-1 characters in a text of length n. This takes
# time O(n). With every character we update its count in a table C and the position
# of the stretch in constant time.
# 
# There are two worst case string combination scenarios - a lexicographically
# descending sequence and one where at every alternating position (i+1) the character
# is lexicographically smaller. We then have to compare table C against M and
# if the counts are greater then copy it to M. Finally, we have to reset C.
# These actions take time O(a).
#
# So, for every position in T or every other position in T, we could be doing
# O(a) operations, so the worst-case time complexity is O(na); Since the
# alphabet is constant in a fixed-size alphabet (a = O(1)), the time complexity
# could be said to be linear, O(n).
# 
# The space complexity of the algorithm is O(a) since we require to store counts
# in tables C and M.
#

import argparse


def _get_letter_idx(c, A):
    return A.index(c)


def _check_counts_better(min_counts, curr_counts, a):
    j = 0
    while j < a:
        if curr_counts[j] > min_counts[j]:
            return True
        elif curr_counts[j] < min_counts[j]:
            return False
        j += 1
    return False


def min_lex_rot_count(seq, A):
    n = len(seq)
    a = len(A)
    start_pos = 0
    end_pos = 0
    min_start_pos = 0
    min_end_pos = 0
    min_counts = [0] * a  # M
    curr_counts = [0] * a  # C

    # first letter
    ltr_idx = _get_letter_idx(seq[0], A)
    curr_counts[ltr_idx] = 1
    min_counts[ltr_idx] = 1

    i = 1
    while i < 2*n-1:
        prev_ltr = ltr_idx
        ltr_idx = _get_letter_idx(seq[i%n], A)

        # if curr_letter is not the same/bigger than previous letter
        if prev_ltr > ltr_idx:

            # check last stretch is greater than current one and if so copy it
            if _check_counts_better(min_counts, curr_counts, a):
                j = 0
                while j < a:
                    min_counts[j] = curr_counts[j]
                    j += 1
                min_start_pos = start_pos
                min_end_pos = end_pos

            # reset counts
            j = 0
            while j < a:
                curr_counts[j] = 0
                j += 1
            curr_counts[ltr_idx] = 1
            start_pos = i

        else:

            # count letter
            curr_counts[ltr_idx] += 1

        end_pos = i

        i += 1

    # check last stretch is greater than current one and if so update the min positions
    if _check_counts_better(min_counts, curr_counts, a):
        min_start_pos = start_pos
        min_end_pos = end_pos

    return min_start_pos, min_end_pos


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='min_lex_rot_count algorithm for finding the Lexicographically Minimal String Rotation')
    parser.add_argument('-a', '--alphabet', type=str, required=True, help='The alphabet, e.g. "ACGT"')
    parser.add_argument('-t', '--text', type=str, required=True, help='The input string')
    args = parser.parse_args()
    
    T = args.text
    A = args.alphabet
    n = len(T)
    start_pos, end_pos = min_lex_rot_count(T, A)
    # print(T)
    # print(start_pos, end_pos)
    TT = T + T
    print(TT[start_pos:start_pos+n])
