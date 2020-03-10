import itertools


def hamming(s1, s2):
    """Calculate the Hamming distance between two bit strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


max_len = 0
perms = itertools.product('ACGT', repeat=3)
for p in perms:
    s = set()
    s.add(''.join(p))

    for word in itertools.product('ACGT', repeat=3):
        if all(x == word[0] for x in word):
            continue
        word = ''.join(word)
        distances = [hamming(word, item) >= 2 for item in s]
        if all(distances):
            s.add(word)

    print(len(s))
    if len(s) > max_len:
        max_len = len(s)
        winner = s
