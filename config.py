import itertools
import pathlib
from typing import Union
from compoiste_algorithm import CompositeAlgorithm
from k_mer_algorithm import KMerAlgorithm

PathLike = Union[str, pathlib.Path]


def build_config(subset_size: int = 5,
                 bits_per_z: int = 12,
                 letter_replace_error_ratio: int = 0,
                 letter_remove_error_ratio: int = 0,
                 letter_add_error_ratio: int = 0,
                 number_of_oligos_per_barcode: int = 20,
                 input_text_file: PathLike = pathlib.Path(r'data/testing/input_text.dna'),
                 output_dir: PathLike = pathlib.Path(r'data/testing/output')):

    pathlib.Path('data/testing').mkdir(parents=True, exist_ok=True)
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

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

    shrink_dict_size = len(shrink_dict_3_mer)

    k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, shrink_dict_size + 1)], subset_size)
    r = [set(k) for k in k_mer_representative]
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
        # 'mode': 'prod',
        'mode': 'test',
        'k_mer': 3,
        'oligo_len_binary': None,
        'shrink_dict': shrink_dict_3_mer,
        'oligo_file_name': 'Oligo_Input',
        'fastq_file_name': 'Bible4_sample',
        'file_name_sorted': pathlib.Path(r'data/testing/small_data_3_barcode_9_oligo.dna'),
        'input_text_file': input_text_file,
        'binary_file_name': pathlib.Path(r'data/testing/small_data.binary.dna'),
        'encoder_results_file': pathlib.Path(
            r'data/testing/small_data_binary.encoder_results_file.dna'),
        'synthesis_results_file': pathlib.Path(
            r'data/testing/small_data_binary.synthesis_results_file.dna'),
        'decoder_results_file': pathlib.Path(
            r'data/testing/small_data_binary.decoder_results_file.dna'),
        'binary_results_file': pathlib.Path(r'data/testing/small_data.binary_results_file.dna'),
        'text_results_file': pathlib.Path(r'data/testing/small_data.text_results_file.dna'),
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
        'synthesis': {'number_of_oligos_per_barcode': number_of_oligos_per_barcode,
                      'letter_replace_error_ratio': letter_replace_error_ratio,
                      'letter_remove_error_ratio': letter_remove_error_ratio,
                      'letter_add_error_ratio': letter_add_error_ratio,
                      'seed': 0
                      }

    }

    if config['mode'] == 'prod':
        config['barcode_len'] = 12
        config['barcode_rs_len'] = 4
        config['payload_len'] = 120
        config['payload_rs_len'] = 14
    elif config['mode'] == 'test':
        config['barcode_len'] = 12
        config['barcode_rs_len'] = 4
        config['payload_len'] = 9
        config['payload_rs_len'] = 3*3
    
    config['barcode_total_len'] = config['barcode_len'] + config['barcode_rs_len']
    config['payload_total_len'] = config['payload_len'] + config['payload_rs_len']

    config['oligo_len_binary'] = int(
        config['payload_len'] / config['k_mer'] * config['algorithm_config']['bits_per_z'])

    return config


config = build_config()
