from typing import Tuple
import itertools


def dna_sequence_generator(sequence_len=12, symbols=('A', 'C', 'G', 'T')) -> Tuple[str]:
    barcodes = itertools.product(symbols, repeat=sequence_len)
    while True:
        try:
            yield next(barcodes)
        except StopIteration:
            return