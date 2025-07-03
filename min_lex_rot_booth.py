# Booth's lexicographically minimal string rotation algorithm as found on Wikipedia:
# https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation

def min_lex_rot_booth(text):
    n = len(text)
    f = [-1] * (2 * n)
    k = 0
    for j in range(1, 2 * n):
        i = f[j - k - 1]
        while i != -1 and text[j % n] != text[(k + i + 1) % n]:
            if text[j % n] < text[(k + i + 1) % n]:
                k = j - i - 1
            i = f[i]
        if i == -1 and text[j % n] != text[(k + i + 1) % n]:
            if text[j % n] < text[(k + i + 1) % n]:
                k = j
            f[j - k] = -1
        else:
            f[j - k] = i + 1
    return k


if __name__ == '__main__':
    text = 'MISSISSIPPI'
    print(min_lex_rot_booth(text))  # 10 (IMISSISSIPP)
