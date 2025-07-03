# Copyright (c) 2025 Ahmad Retha; MIT License.
#
# A naive implementation of the lexicographically minimal string rotation algorithm. It requires O(n^2) time and O(n)
# space in this implementation.


def min_lex_rot_naive(text):
    n = len(text)
    text2 = text + text
    min_i = 0

    for i in range(1, n):
        if text2[i:i+n] < text2[min_i:min_i+n]:
            min_i = i

    return min_i


if __name__ == '__main__':
    text = 'MISSISSIPPI'
    print(min_lex_rot_naive(text))  # 10 (IMISSISSIPP)
