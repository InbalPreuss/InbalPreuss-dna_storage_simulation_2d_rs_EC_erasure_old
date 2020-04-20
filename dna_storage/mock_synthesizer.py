from pathlib import Path
import random
from typing import Union, Dict, List

import numpy as np


class Synthesizer:
    def __init__(self, input_file: Union[Path, str],
                 results_file: Union[Path, str],
                 synthesis_config: Dict,
                 barcode_total_len: int,
                 subset_size: int,
                 k_mer_representative_to_z: Dict,
                 k_mer_to_dna: Dict,
                 mode: str):
        self.input_file = input_file
        self.results_file = results_file
        open(self.results_file, 'w').close()
        self.synthesis_config = synthesis_config
        self.barcode_total_len = barcode_total_len
        self.subset_size = subset_size
        self.k_mer_representative_to_z = k_mer_representative_to_z
        self.k_mer_to_dna = k_mer_to_dna
        self.mode = mode

    def synthesize(self):
        if self.mode == 'test':
            np.random.seed(self.synthesis_config['seed'])
            random.seed(self.synthesis_config['seed'])
        with open(self.input_file, 'r', encoding='utf-8') as input_file, open(self.results_file, 'w+', encoding='utf-8') as results_file:
            for line in input_file:
                line_list = line.strip('\n').split(',')
                barcode, payload = line_list[0], line_list[1:]
                x_list = self.get_x_list(payload=payload)
                number_of_nuc = max(1, int(round(np.random.normal(self.synthesis_config['number_of_oligos_per_barcode'], 2))))
                x_mat = np.empty([number_of_nuc, self.barcode_total_len], dtype=np.dtype(('U', 5)))
                x_mat[:] = np.array(tuple(barcode))
                for idx, x_tuple in enumerate(x_list, 1):
                    vec = np.random.choice(x_tuple, size=(number_of_nuc,))
                    if self.mode == 'test':
                        while True:
                            if len(np.unique(vec)) == self.subset_size:
                                break
                            vec = np.random.choice(x_tuple, size=(number_of_nuc,))
                    col = np.array([tuple(self.k_mer_to_dna[k_mer]) for k_mer in vec])
                    x_mat = np.hstack((x_mat, col))

                dna_list = [''.join(row) for row in x_mat]
                dna_list = self.add_remove_replace(dna_list)
                results_file.write('\n'.join(dna_list) + '\n')

    def add_remove_replace(self, dna_list: List[str]):
        for row_idx, oligo in enumerate(dna_list):
            remove = [np.random.binomial(1, self.synthesis_config['letter_remove_error_ratio']) for i in range(len(oligo))]
            oligo = ''.join([char if remove[idx] == 0 else '' for idx, char in enumerate(oligo)])
            add_idx = [np.random.binomial(1, self.synthesis_config['letter_add_error_ratio']) for i in range(len(oligo))]
            add = [random.choice('ACGT') if i == 1 else '' for i in add_idx]
            oligo = ''.join(''.join(x) for x in zip(oligo, add))
            replace_idx = [np.random.binomial(1, self.synthesis_config['letter_replace_error_ratio']) for i in
                           range(len(oligo))]
            for letter_idx, letter in enumerate(oligo):
                if replace_idx == 1:
                    diff = {'A', 'C', 'G', 'T'} - set(letter)
                    oligo[letter_idx] = random.choice(''.join(diff))
            dna_list[row_idx] = oligo
        return dna_list

    def get_x_list(self, payload: List[str]):
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
