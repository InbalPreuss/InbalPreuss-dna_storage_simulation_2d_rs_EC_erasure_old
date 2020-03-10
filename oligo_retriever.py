from collections import Counter
from typing import Union, Dict
from pathlib import Path

from compoiste_algorithm import CompositeAlgorithm
from k_mer_algorithm import KMerAlgorithm

#################################################################
# @ Class: OligoRetriever
# @ Description: Retrieve the oligo to the oligo that was written
#                in originally
#################################################################
class OligoRetriever:
    def __init__(self, number_of_barcode_letters: int,
                 oligo_length: int,
                 oligo_sorted_file_name: str,
                 algorithm: Union[CompositeAlgorithm, KMerAlgorithm],
                 shrink_dict: Dict,
                 k_mer: int,
                 k_mer_representative_to_z: Dict,
                 z_to_binary: Dict,
                 unique_oligo_results_file: Union[Path, str],
                 ):
        self.file_name = oligo_sorted_file_name
        self.number_of_barcode_letters = number_of_barcode_letters
        self.oligo_length = oligo_length
        self.algorithm = algorithm
        self.shrink_dict = shrink_dict
        self.k_mer = k_mer
        self.k_mer_representative_to_z = k_mer_representative_to_z
        self.z_to_binary = z_to_binary
        self.binary_results_file = open(unique_oligo_results_file, 'w+')

    def __del__(self):
        self.binary_results_file.close()

    def run(self):
        barcode_prev = ''
        payload_accumulation = []
        with open(self.file_name, 'r') as file:
            for line in file:
                barcode_and_payload = line.split(sep=' ')[0].rstrip()
                barcode, payload = barcode_and_payload[:self.number_of_barcode_letters], barcode_and_payload[self.number_of_barcode_letters:]
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
            for key, val in self.z_to_binary.items():
                if key == z:
                    binary.append(val)
                    break
        binary = ["".join(map(str, tup)) for tup in binary]
        return "".join(binary)

    def get_unique_payload(self, payload_accumulation):
        payload_accumulation_histogram = self.payload_histogram(oligo_accumulation=payload_accumulation)
        unique_payload_accumulation = self.payload_histogram_to_payload(hist_accumulation=payload_accumulation_histogram)
        unique_payload_accumulation = self.error_correction(unique_payload_accumulation)
        return unique_payload_accumulation

    def get_shrunk_unique_payload(self, unique_payload_accumulation):
        shrunk_unique_payload = self.shrink_payload(payload_accumulation=unique_payload_accumulation)
        # barcode_histogram = self.barcode_histogram(shrinked_oligo=shrinked_oligo)
        return shrunk_unique_payload

    def shrink_payload(self, payload_accumulation):
        """ Note that missing k-mers will be removed from the oligo_accumulation """
        if self.k_mer == 1:
            return payload_accumulation
        k_mer_accumulation = []
        for payload in payload_accumulation:
            k_mer_list = []
            oligo_valid = True
            for k_letters in [payload[i:i+self.k_mer] for i in range(0, self.oligo_length, self.k_mer)]:
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
        for col_idx in range(int(self.oligo_length/self.k_mer)):
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
            k_mer_rep = set([rep[0] for rep in reps])
            for key, val in self.k_mer_representative_to_z.items():
                if set(key) == k_mer_rep:
                    result_payload.append(val)
                    break
        return result_payload

    def save_binary(self, binary, barcode_prev):
        self.binary_results_file.write(barcode_prev + binary + '\n')

