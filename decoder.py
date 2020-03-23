from collections import Counter
import re
from typing import Union, Dict
from pathlib import Path

from compoiste_algorithm import CompositeAlgorithm
from k_mer_algorithm import KMerAlgorithm


#################################################################
# @ Class: Decoder
# @ Description: Retrieve the oligo to the oligo that was written
#                in originally
#################################################################
class Decoder:
    def __init__(self, number_of_barcode_letters: int,
                 oligo_length: int,
                 input_file: str,
                 algorithm: Union[CompositeAlgorithm, KMerAlgorithm],
                 shrink_dict: Dict,
                 k_mer: int,
                 k_mer_representative_to_z: Dict,
                 z_to_binary: Dict,
                 results_file: Union[Path, str],
                 ):
        self.input_file = input_file
        self.number_of_barcode_letters = number_of_barcode_letters
        self.oligo_length = oligo_length
        self.algorithm = algorithm
        self.shrink_dict = shrink_dict
        self.k_mer = k_mer
        self.k_mer_representative_to_z = k_mer_representative_to_z
        self.z_to_binary = z_to_binary
        self.results_file = results_file
        open(self.results_file, 'w').close()

    def run(self):
        barcode_prev = ''
        payload_accumulation = []
        with open(self.input_file, 'r', encoding='utf-8') as file:
            for line in file:
                barcode_and_payload = line.split(sep=' ')[0].rstrip()
                barcode, payload = barcode_and_payload[:self.number_of_barcode_letters], barcode_and_payload[
                                                                                         self.number_of_barcode_letters:]
                if barcode != barcode_prev:
                    if len(payload_accumulation) != 0:
                        binary = self.dna_to_binary(payload_accumulation=payload_accumulation)
                        self.save_binary(binary=binary, barcode_prev=barcode_prev)
                    payload_accumulation = [payload]
                    barcode_prev = barcode
                else:
                    payload_accumulation.append(payload)
            binary = self.dna_to_binary(payload_accumulation=payload_accumulation)
            self.save_binary(binary=binary, barcode_prev=barcode_prev)

    def dna_to_binary(self, payload_accumulation):
        shrunk_payload = self.shrink_payload(payload_accumulation=payload_accumulation)
        shrunk_payload_histogram = self.payload_histogram(payload=shrunk_payload)
        unique_payload = self.payload_histogram_to_payload(payload_histogram=shrunk_payload_histogram)
        unique_payload = self.error_correction(payload=unique_payload)
        binary = self.unique_payload_to_binary(payload=unique_payload)
        return binary

    def unique_payload_to_binary(self, payload):
        binary = []
        for z in payload:
            try:
                binary.append(self.z_to_binary[z])
            except KeyError:
                return ''
        binary = ["".join(map(str, tup)) for tup in binary]
        return "".join(binary)

    def shrink_payload(self, payload_accumulation):
        """ Note that missing k-mers will be removed from the oligo_accumulation """
        if self.k_mer == 1:
            return payload_accumulation
        k_mer_accumulation = []
        for payload in payload_accumulation:
            k_mer_list = []
            oligo_valid = True
            for k_letters in [payload[i:i + self.k_mer] for i in range(0, self.oligo_length, self.k_mer)]:
                try:
                    k_mer_list.append(self.shrink_dict[k_letters])
                except KeyError:
                    oligo_valid = False
                    break
            if oligo_valid:
                k_mer_accumulation.append(k_mer_list)
        return k_mer_accumulation

    def payload_histogram(self, payload):
        hist = []
        for col_idx in range(int(self.oligo_length / self.k_mer)):
            col = [letter[col_idx] for letter in payload]
            letter_counts = Counter(col)
            hist.append(letter_counts)
        return hist

    def error_correction(self, payload):
        return payload

    def payload_histogram_to_payload(self, payload_histogram):
        result_payload = []
        for counter in payload_histogram:
            reps = counter.most_common(self.algorithm.subset_size)
            if len(reps) != self.algorithm.subset_size:
                return []
            k_mer_rep = tuple(self.sorted_human([rep[0] for rep in reps]))
            result_payload.append(self.k_mer_representative_to_z[k_mer_rep])
        return result_payload

    def save_binary(self, binary, barcode_prev):
        with open(self.results_file, 'a+', encoding='utf-8') as f:
            f.write(barcode_prev + binary + '\n')

    @staticmethod
    def sorted_human(iterable):
        """ Sort the given iterable in the way that humans expect."""
        convert = lambda text: int(text) if text.isdigit() else text
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(iterable, key=alphanum_key)
