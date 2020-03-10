import itertools
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
                 z_to_binary: Dict,
                 results_file: Union[Path, str],
                 ):
        self.file_name = binary_file_name
        self.number_of_barcode_letters = number_of_barcode_letters
        self.oligo_length = oligo_length
        self.algorithm = algorithm
        self.shrink_dict = shrink_dict
        self.k_mer = k_mer
        self.k_mer_representative_to_z = k_mer_representative_to_z
        self.z_to_binary = z_to_binary
        self.results_file = open(results_file, 'w+')
        self.barcode_generator = self.get_barcode_generator()

    def __del__(self):
        self.results_file.close()

    def run(self):
        with open(self.file_name, 'r') as file:
            z_list = []
            for line in file:
                line = line.strip('\n')
                if len(line) < 12:
                    pass  # TODO
                if len(z_list) < int(self.oligo_length / self.k_mer):
                    z = self.binary_to_z(binary=line)
                    z_list.append(z)
                else:
                    oligo = self.z_to_oligo(z_list)
                    self.save_oligo(oligo=oligo)
                    z_list = []
            if len(z_list) > 0:
                oligo = self.z_to_oligo(z_list)
                self.save_oligo(oligo=oligo)

    def binary_to_z(self, binary):
        binary_tuple = tuple([int(b) for b in binary])
        for key, val in self.z_to_binary.items():
            if val == binary_tuple:
                return key

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
        self.results_file.write(oligo + '\n')


