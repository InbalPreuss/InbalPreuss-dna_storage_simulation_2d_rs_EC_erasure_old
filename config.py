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

subset_size = 5
shrink_dict_size = len(shrink_dict_3_mer)
bits_per_z = 12

k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, shrink_dict_size + 1)], subset_size)
r = [set(k) for k in k_mer_representative]
k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, shrink_dict_size + 1)], subset_size)
all_binary_combinations = itertools.product([0, 1], repeat=bits_per_z)
z = itertools.combinations(['Z' + str(i) for i in range(1, len(r) + 1)], 1)
z = [i[0] for i in z]
z_to_binary = dict(zip(z, all_binary_combinations))
all_binary_combinations = itertools.product([0, 1], repeat=bits_per_z)
binary_to_z = dict(zip(all_binary_combinations, z))

k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, shrink_dict_size + 1)], subset_size)
k_mer_representative_to_z = dict(zip(k_mer_representative, z))
k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, shrink_dict_size + 1)], subset_size)
z_to_k_mer_representative = dict(zip(z, k_mer_representative))

k_mer_to_dna = {v: k for k, v in shrink_dict_3_mer.items()}

config = {
    # 'NUMBER_OF_BARCODE_LETTERS': 16,
    # 'OLIGO_LENGTH': 150,
    'NUMBER_OF_BARCODE_LETTERS': 3,
    'OLIGO_LENGTH': 9,
    'K_MER': 3,
    'shrink_dict': shrink_dict_3_mer,
    'OLIGO_FILE_NAME': 'Oligo_Input',
    'FASTQ_FILE_NAME': 'Bible4_sample',
    'file_name_sorted': pathlib.Path(r'data/testing/small_data_3_barcode_9_oligo.dna'),
    'input_text_file': pathlib.Path(r'data/testing/small_data.input_text.dna'),
    'binary_file_name': pathlib.Path(r'data/testing/small_data.binary.dna'),
    'encoder_results_file': pathlib.Path(
            r'data/testing/small_data_binary.encoder_results_file.dna'),
    'synthesis_results_file': pathlib.Path(
                r'data/testing/small_data_binary.synthesis_results_file.dna'),
    'decoder_results_file': pathlib.Path(
            r'data/testing/small_data_binary.decoder_results_file.dna'),
    'binary_results_file': pathlib.Path(r'data/testing/small_data.binary_results_file.dna'),
    'text_results_file':  pathlib.Path(r'data/testing/small_data.text_results_file.dna'),
    'do_oligo_handling': False,
    'write_text_to_binary': True,
    'do_encode': True,
    'do_synthesize': True,
    'do_fastq_handling': False,
    'do_decode': True,
    'decoder_results_to_binary': True,
    'binary_results_to_text': True,
    'algorithm': KMerAlgorithm,
    'algorithm_config': {'subset_size': subset_size,
                         'bits_per_z': bits_per_z,
                         'shrink_dict_size': shrink_dict_size,
                         'k_mer_representative_to_z': k_mer_representative_to_z,
                         'z_to_k_mer_representative': z_to_k_mer_representative,
                         'z_to_binary': z_to_binary,
                         'binary_to_z': binary_to_z,
                         'k_mer_to_dna': k_mer_to_dna},
    'synthesis': {'number_of_oligos_per_barcode': 20,
                  'letter_replace_error_ratio': 0,
                  'letter_remove_error_ratio': 0,
                  'letter_add_error_ratio': 0,
                  'seed': 0
                  }

}
