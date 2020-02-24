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
                 unique_oligo_results_file: Union[Path, str]
                 ):
        self.file_name = oligo_sorted_file_name
        self.number_of_barcode_letters = number_of_barcode_letters
        self.oligo_length = oligo_length
        self.algorithm = algorithm
        self.shrink_dict = shrink_dict
        self.k_mer = k_mer
        self.unique_oligo_results_file = open(unique_oligo_results_file, 'w+')

    def __del__(self):
        self.unique_oligo_results_file.close()

    def run(self):
        barcode_prev = ''
        oligo_accumulation = []
        with open(self.file_name, 'r') as file:
            for line in file:
                barcode_and_oligo = line.split(sep=' ')[0].rstrip()
                barcode, oligo = barcode_and_oligo[:self.number_of_barcode_letters], barcode_and_oligo[self.number_of_barcode_letters:]
                if barcode != barcode_prev:
                    if len(oligo_accumulation) != 0:
                        unique_oligo_accumulation = self.retrieve_unique_oligo(oligo_accumulation=oligo_accumulation)
                        for unique_oligo in unique_oligo_accumulation:
                            self.save_unique_oligo(unique_oligo, barcode_prev)
                    oligo_accumulation = [oligo]
                    barcode_prev = barcode
                else:
                    oligo_accumulation.append(oligo)
            unique_oligo = self.retrieve_unique_oligo(oligo_accumulation=oligo_accumulation)
            for unique_oligo in unique_oligo_accumulation:
                self.save_unique_oligo(unique_oligo, barcode_prev)

    def retrieve_unique_oligo(self, oligo_accumulation):
        unique_oligo_accumulation = self.get_unique_oligo(oligo_accumulation=oligo_accumulation)
        shrinked_unique_oligos = self.get_shrinked_unique_oligo(unique_oligo_accumulation=unique_oligo_accumulation)
        binary_data = self.algorithm.decode(shrinked_oligos=shrinked_unique_oligos)
        return unique_oligo_accumulation

    def get_unique_oligo(self, oligo_accumulation):
        oligo_accumulation_histogram = self.oligo_histogram(oligo_accumulation=oligo_accumulation)
        unique_oligo_accumulation = self.hist_accumulation_to_unique_oligo_accumulation(hist_accumulation=oligo_accumulation_histogram)
        unique_oligo_accumulation = self.error_correction(unique_oligo_accumulation)
        return unique_oligo_accumulation

    def get_shrinked_unique_oligo(self, unique_oligo_accumulation):
        shrinked_unique_oligos = self.shrink_oligo(oligo_accumulation=unique_oligo_accumulation)
        # barcode_histogram = self.barcode_histogram(shrinked_oligo=shrinked_oligo)
        return shrinked_unique_oligos

    def shrink_oligo(self, oligo_accumulation):
        """ Note that missing k-mers will be removed from the oligo_accumulation """
        if self.k_mer == 1:
            return oligo_accumulation
        k_mer_accumulation = []
        for oligo in oligo_accumulation:
            k_mer_list = []
            oligo_valid = True
            for k_letters in [oligo[i:i+self.k_mer] for i in range(0, self.oligo_length, self.k_mer)]:
                try:
                    k_mer_list.append(self.shrink_dict[k_letters])
                except KeyError:
                    oligo_valid = False
                    break
            if oligo_valid:
                k_mer_accumulation.append(k_mer_list)
        return k_mer_accumulation

    def barcode_histogram(self, shrinked_oligo):
        hist = []
        for col_idx in range(int(self.oligo_length/self.k_mer)):
            col = [letter[col_idx] for letter in shrinked_oligo]
            letter_counts = Counter(col)
            hist.append(letter_counts)
        return hist

    def oligo_histogram(self, oligo_accumulation):
        return Counter(tuple(item) for item in oligo_accumulation)

    def error_correction(self, oligo):
        return oligo

    def hist_accumulation_to_unique_oligo_accumulation(self, hist_accumulation):
        list_of_hists = hist_accumulation.most_common(self.algorithm.binary_bits_on)
        return [''.join(tup[0]) for tup in list_of_hists]

    def hist_to_unique_oligo(self, barcode_histogram):
        unique_oligo = self.algorithm.decode(barcode_histogram=barcode_histogram)
        return unique_oligo

    def save_unique_oligo(self, unique_oligo, barcode_prev):
        self.unique_oligo_results_file.write(barcode_prev + unique_oligo + '\n')

