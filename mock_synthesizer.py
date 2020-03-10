from pathlib import Path
import random
from typing import Union, Dict, List


class Synthesizer:
    def __init__(self, input_file: Union[Path, str],
                 results_file: Union[Path, str],
                 number_of_oligos_per_barcode: int,
                 letter_error_ratio: float,
                 k_mer_representative_to_z: Dict):
        self.input_file = input_file
        self.results_file = open(results_file, 'w+')
        self.number_of_oligos_per_barcode = number_of_oligos_per_barcode
        self.letter_error_ratio = letter_error_ratio
        self.k_mer_representative_to_z = k_mer_representative_to_z

    def synthesize(self):
        with open(self.input_file, 'r') as file:
            for line in file:
                line_list = line.strip('\n').split(',')
                barcode, payload = line_list[0], line_list[1:]
                x_list = self.get_x_list(payload=payload)
                a = 3
                for x_tuple in x_list:
                    number_of_nuc = random.randint(self.number_of_oligos_per_barcode*0.7, self.number_of_oligos_per_barcode*1.3)
                    samples = self.constrained_sum_sample_pos(n=2, total=number_of_nuc)
                    a = 3

    def get_x_list(self, payload):
        x_list = []
        for z in payload:
            for key, val, in self.k_mer_representative_to_z.items():
                if val == z:
                    x_list.append(key)
        return x_list

    def constrained_sum_sample_pos(self, n, total):
        """Return a randomly chosen list of n positive integers summing to total.
        Each such list is equally likely to occur."""

        dividers = sorted(random.sample(range(1, total), n - 1))
        return [a - b for a, b in zip(dividers + [total], [0] + dividers)]