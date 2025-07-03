# MIN_LEX_ROT: Lexicographically Minimal String Rotation of short strings

The naive algorithm for solving the [Lexicographically Minimal String Rotation](https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation)
(LMSR) problem involves iterating through successive rotations of a string and
comparing each rotation to the previously identified lexicographically minimal
rotation and updating the index accordingly. It requires quadratic time to run
in the length of the string, O(n^2).

An algorithm for solving the LMSR problem was devised by Booth in 1980 and was
based on the _KMP_ search algorithm. It requires linear time and space. In 1983
Duval devised an algorithm based on _Lyndon_ words that required linear time and
constant space.

Here you will find two new LMSR algorithms (KR and BW) -- one inspired by the
_Karp-Rabin_ algorithm and the other based on _bitwise_ packing of a string into
a computer word (typically, a 64-bit int). The premise of both algorithms is to
first pack each letter of a text into a computer word, then rotate it O(n) times,
calculating a lexicographical score each time. The score is the numerical value
of the integer stored in the computer word and can be checked in constant time.
Rotating the computer word can also be done in constant time. Thus, the algorithms
require linear time and constant space when the text fits into a single computer
word, but in practice this means it is limited to short strings using a small
constant alphabet.

In this repository you will find the KR and BW single-word implementations
written in Python:

	- min_lex_rot_kr.py
	- min_lex_rot_bw.py

Also in this repository, you will find a generalised implementation, which
can handle longer strings by using more memory, but consequently _doesn't_ do
it in linear time and constant space. Please refer to the comments of the
Python implementation for an analysis of the time and space complexity of that
algorithm. There is also an implementation in C for greater performance:

	- generalised_min_lex_rot.py
	- generalised_min_lex_rot.c

For long strings, I recommend using Booth's or Duval's algorithm, but for short
ones, try the algorithms provided here.

The code in this repository is provided only as a proof of concept with no
guarantees. Please refer to the license.

Copyright (c) 2025 Ahmad Retha; MIT License.
