import itertools
import pathlib
from compoiste_algorithm import CompositeAlgorithm
from k_mer_algorithm import KMerAlgorithm

shrink_dict_3_mer = {'AAT': 'X1',
                     'ACA': 'X2',
                     'ATG': 'X3',
                     'AGC': 'X4',
                     'TAA': 'X5',
                     'TCT': 'X6',
                     'TTC': 'X7',
                     'TGG': 'X8',
                     'GAG': 'X9',
                     'GCC': 'X10',
                     'GTT': 'X11',
                     'GGA': 'X12',
                     'CAC': 'X13',
                     'CCG': 'X14',
                     'CTA': 'X15',
                     'CGT': 'X16'}

k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, 16 + 1)], 2)
r = [set(k) for k in k_mer_representative]
k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, 16 + 1)], 2)
all_binary_combinations = itertools.product([0, 1], repeat=12)
z = itertools.combinations(['Z' + str(i) for i in range(1, len(r) + 1)], 1)
z = [i[0] for i in z]
z_to_binary = dict(zip(z, all_binary_combinations))

k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, 16 + 1)], 2)
k_mer_representative_to_z = dict(zip(k_mer_representative, z))
# binary_to_k_mer_representative = dict(zip(all_binary_combinations, k_mer_representative))

config = {
    # 'NUMBER_OF_BARCODE_LETTERS': 16,
    # 'OLIGO_LENGTH': 151,
    'NUMBER_OF_BARCODE_LETTERS': 3,
    'OLIGO_LENGTH': 9,
    'K_MER': 3,
    'shrink_dict': shrink_dict_3_mer,
    'OLIGO_FILE_NAME': 'Oligo_Input',
    'FASTQ_FILE_NAME': 'Bible4_sample',
    'file_name_sorted': pathlib.Path(r'data/testing/small_data_3_barcode_9_oligo.dna'),
    'unique_oligo_results_file': pathlib.Path(
        r'data/testing/small_data_3_barcode_9_oligo.unique_oligo_results_file.dna'),
    'do_fastq_handling': False,
    'do_retrieved_oligo': True,
    'do_oligo_handling': False,
    'algorithm': KMerAlgorithm,
    'algorithm_config': {'subset_size': 2,
                         'k_mer_representative_to_z': k_mer_representative_to_z,
                         'z_to_binary': z_to_binary}

}
