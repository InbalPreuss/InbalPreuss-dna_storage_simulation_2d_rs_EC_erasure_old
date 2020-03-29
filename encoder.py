import itertools
from textwrap import wrap
from collections import Counter
from typing import Union, Dict
from pathlib import Path

from compoiste_algorithm import CompositeAlgorithm
from k_mer_algorithm import KMerAlgorithm


#################################################################
# @ Class: Encoder
# @ Description: Retrieve the oligo to the oligo that was written
#                in originally
#################################################################
class Encoder:
    def __init__(self, number_of_barcode_letters: int,
                 oligo_length: int,
                 binary_file_name: str,
                 algorithm: Union[CompositeAlgorithm, KMerAlgorithm],
                 shrink_dict: Dict,
                 k_mer: int,
                 k_mer_representative_to_z: Dict,
                 binary_to_z: Dict,
                 bits_per_z: int,
                 results_file: Union[Path, str],
                 ):
        self.file_name = binary_file_name
        self.number_of_barcode_letters = number_of_barcode_letters
        self.oligo_length = oligo_length
        self.algorithm = algorithm
        self.shrink_dict = shrink_dict
        self.k_mer = k_mer
        self.k_mer_representative_to_z = k_mer_representative_to_z
        self.binary_to_z_dict = binary_to_z
        self.bits_per_z = bits_per_z
        self.results_file = results_file
        open(self.results_file, 'w').close()
        self.barcode_generator = self.get_barcode_generator()

    def run(self):
        with open(self.file_name, 'r', encoding='utf-8') as file:
            z_list = []
            z_list_len = int(self.oligo_length / self.k_mer)
            for line in file:
                line = line.strip('\n')
                for binary_to_transform in wrap(line, self.bits_per_z):
                    z = self.binary_to_z(binary=binary_to_transform)
                    z_list.append(z)
                    if len(z_list) == z_list_len:
                        oligo = self.z_to_oligo(z_list)
                        self.save_oligo(oligo=oligo)
                        z_list = []

            if len(z_list) > 0:
                oligo = self.z_to_oligo(z_list)
                self.save_oligo(oligo=oligo)

    def binary_to_z(self, binary):
        binary_tuple = tuple([int(b) for b in binary])
        return self.binary_to_z_dict[binary_tuple]

    def z_to_oligo(self, z_list):
        oligo = self.error_correction(payload=z_list)
        barcode = next(self.barcode_generator)
        barcode = self.error_correction(payload=barcode)
        barcode = "".join(barcode)
        oligo.insert(0, barcode)
        return ",".join(oligo)

    def get_barcode_generator(self):
        barcodes = itertools.product(['A', 'C', 'G', 'T'], repeat=self.number_of_barcode_letters)
        while True:
            try:
                yield next(barcodes)
            except StopIteration:
                return

    def error_correction(self, payload):
        return payload

    def save_oligo(self, oligo):
        with open(self.results_file, 'a+', encoding='utf-8') as f:
            f.write(oligo + '\n')


