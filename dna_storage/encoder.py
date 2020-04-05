import itertools
from textwrap import wrap
from typing import Union, Dict, List, Tuple
from pathlib import Path

from dna_storage.reedsolomon import barcode_rs_encode


#################################################################
# @ Class: Encoder
# @ Description: Retrieve the oligo to the oligo that was written
#                in originally
#################################################################

class Encoder:
    def __init__(self, barcode_len: int,
                 barcode_rs_len: int,
                 payload_len: int,
                 payload_rs_len: int,
                 binary_file_name: str,
                 shrink_dict: Dict,
                 k_mer: int,
                 k_mer_representative_to_z: Dict,
                 binary_to_z: Dict,
                 subset_size: int,
                 bits_per_z: int,
                 rs_encoders: Dict,
                 results_file: Union[Path, str],
                 ):
        self.file_name = binary_file_name
        self.barcode_len = barcode_len
        self.barcode_rs_len = barcode_rs_len
        self.payload_len = payload_len
        self.payload_rs_len = payload_rs_len
        self.shrink_dict = shrink_dict
        self.k_mer = k_mer
        self.k_mer_representative_to_z = k_mer_representative_to_z
        self.binary_to_z_dict = binary_to_z
        self.subset_size = subset_size
        self.bits_per_z = bits_per_z
        self.rs_encoders = rs_encoders
        self.results_file = results_file
        open(self.results_file, 'w').close()
        self.barcode_generator = self.get_barcode_generator()

    def run(self):
        with open(self.file_name, 'r', encoding='utf-8') as file:
            z_list = []
            z_list_len = int(self.payload_len)
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

    def binary_to_z(self, binary: str) -> str:
        binary_tuple = tuple([int(b) for b in binary])
        return self.binary_to_z_dict[binary_tuple]

    def z_to_oligo(self, z_list: List[str]) -> str:
        oligo = self.add_payload_rs_symbols_for_error_correction(payload=z_list)
        barcode = next(self.barcode_generator)
        barcode = self.add_barcode_rs_symbols_for_error_correction(barcode=barcode)
        barcode = "".join(barcode)
        oligo.insert(0, barcode)
        return ",".join(oligo)

    def get_barcode_generator(self):
        barcodes = itertools.product(['A', 'C', 'G', 'T'], repeat=self.barcode_len)
        while True:
            try:
                yield next(barcodes)
            except StopIteration:
                return

    def add_payload_rs_symbols_for_error_correction(self, payload: Union[str, List[str]]) -> List[str]:
        if isinstance(payload, str):
            payload = [c for c in payload]
        try:
            encoder = self.select_encoder()
            payload_encoded = encoder(payload)
        except:
            payload_encoded = payload + ['Z1'] * self.payload_rs_len
            # TODO: check possible failure reasons
        return payload_encoded

    def add_barcode_rs_symbols_for_error_correction(self, barcode: Tuple[str]) -> List[str]:
        barcode = list(barcode)
        try:
            barcode_encoded = barcode_rs_encode(barcode)
        except:
            barcode_encoded = barcode + ['A'] * self.barcode_rs_len

        return barcode_encoded

    def select_encoder(self):
        try:
            return self.rs_encoders[self.subset_size]
        except KeyError:
            raise NotImplementedError('No Reed-Solomon is implemented for this subset size')

    def save_oligo(self, oligo: str) -> None:
        with open(self.results_file, 'a+', encoding='utf-8') as f:
            f.write(oligo + '\n')


