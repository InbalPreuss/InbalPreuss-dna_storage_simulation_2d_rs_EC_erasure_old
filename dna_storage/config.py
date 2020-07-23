import itertools
import pathlib
from typing import Union

from dna_storage.rs_adapter import RSBarcodeAdapter, RSPayloadAdapter, RSWideAdapter

PathLike = Union[str, pathlib.Path]


def build_config(subset_size: int = 5,
                 bits_per_z: int = 12,
                 letter_replace_error_ratio: int = 0,
                 letter_remove_error_ratio: int = 0,
                 letter_add_error_ratio: int = 0,
                 number_of_oligos_per_barcode: int = 20,
                 number_of_sampled_oligos_from_file: int = 10000,
                 input_text_file: PathLike = pathlib.Path(r'data/testing/input_text.dna'),
                 output_dir: PathLike = pathlib.Path(r'data/testing/output')):

    pathlib.Path('data/testing').mkdir(parents=True, exist_ok=True)
    pathlib.Path('data/testing/temp').mkdir(parents=True, exist_ok=True)
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
    x_combinations = [set(k) for k in k_mer_representative]
    all_binary_combinations = itertools.product([0, 1], repeat=bits_per_z)
    z = itertools.combinations(['Z' + str(i) for i in range(1, len(x_combinations) + 1)], 1)
    z = [i[0] for i in z]
    z_to_binary = dict(zip(z, all_binary_combinations))
    all_binary_combinations = itertools.product([0, 1], repeat=bits_per_z)
    binary_to_z = dict(zip(all_binary_combinations, z))

    z = z[:2**bits_per_z]
    k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, shrink_dict_size + 1)], subset_size)
    k_mer_representative = list(k_mer_representative)[:2**bits_per_z]
    k_mer_representative_to_z = dict(zip(k_mer_representative, z))
    z_to_k_mer_representative = dict(zip(z, k_mer_representative))

    k_mer_to_dna = {v: k for k, v in shrink_dict_3_mer.items()}

    config = {
        'mode': 'prod',
        # 'mode': 'test',
        'k_mer': 3,
        'number_of_sampled_oligos_from_file': number_of_sampled_oligos_from_file,
        'shrink_dict': shrink_dict_3_mer,
        'fastq_file_name': 'Bible4_sample',
        'file_name_sorted': pathlib.Path(r'data/testing/small_data_3_barcode_9_oligo.dna'),
        'input_text_file': input_text_file,
        'binary_file_name': pathlib.Path(r'data/testing/simulation_data.1.binary.dna'),
        'encoder_results_file': pathlib.Path(
            r'data/testing/simulation_data.2.encoder_results_file.dna'),
        'synthesis_results_file': pathlib.Path(
            r'data/testing/simulation_data.3.synthesis_results_file.dna'),
        'shuffle_db_file': pathlib.Path(r'data/testing/temp_shuffle_db'),
        'shuffle_results_file': pathlib.Path(
            r'data/testing/simulation_data.4.shuffle_results_file.dna'),
        'sample_oligos_results_file': pathlib.Path(
            r'data/testing/simulation_data.5.sample_oligos_results_file.dna'),
        'sort_oligo_db_file': pathlib.Path(r'data/testing/temp_sort_oligo_db'),
        'sort_oligo_results_file': pathlib.Path(
            r'data/testing/simulation_data.6.sort_oligo_results_file.dna'),
        'decoder_results_file': pathlib.Path(
            r'data/testing/simulation_data.7.decoder_results_file.dna'),
        'binary_results_file': pathlib.Path(r'data/testing/simulation_data.8.binary_results_file.dna'),
        'text_results_file': pathlib.Path(r'data/testing/simulation_data.9.text_results_file.dna'),
        'write_text_to_binary': True,
        'do_encode': True,
        'do_synthesize': True,
        'do_shuffle': True,
        'do_sample_oligos_from_file': True,
        'do_sort_oligo_file': True,
        'do_fastq_handling': False,
        'do_decode': True,
        'decoder_results_to_binary': True,
        'binary_results_to_text': True,
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

    wide_n_k = {3: {'block_len': 430, 'block_rs_len': 82},
                5: {'block_len': 3500, 'block_rs_len': 596},
                7: {'block_len': 6900, 'block_rs_len': 1292}}

    if config['mode'] == 'prod':
        config['barcode_len'] = 12  # in ACGT
        config['barcode_rs_len'] = 4  # in ACGT
        config['payload_len'] = 120  # in Z
        config['payload_rs_len'] = 14  # in Z
        config['oligos_per_block_len'] = wide_n_k[subset_size]['block_len']
        config['oligos_per_block_rs_len'] = wide_n_k[subset_size]['block_rs_len'] - 1
        config['number_of_sampled_oligos_from_file'] = number_of_sampled_oligos_from_file * (wide_n_k[subset_size]['block_len'] + (wide_n_k[subset_size]['block_rs_len'] - 1))
    elif config['mode'] == 'test':
        config['barcode_len'] = 12  # in ACGT
        config['barcode_rs_len'] = 4  # in ACGT
        config['payload_len'] = 120  # in Z
        config['payload_rs_len'] = 14  # in Z
        config['oligos_per_block_len'] = 12
        config['oligos_per_block_rs_len'] = 4
        config['number_of_sampled_oligos_from_file'] = config['oligos_per_block_len'] + config['oligos_per_block_rs_len']
    
    config['barcode_total_len'] = config['barcode_len'] + config['barcode_rs_len']  # in ACGT
    config['payload_total_len'] = config['payload_len'] + config['payload_rs_len']  # in Z

    config['barcode_coder'] = RSBarcodeAdapter(bits_per_z=bits_per_z, barcode_len=config['barcode_len'],
                                               barcode_rs_len=config['barcode_rs_len'])
    config['payload_coder'] = RSPayloadAdapter(bits_per_z=bits_per_z, payload_len=config['payload_len'],
                                               payload_rs_len=config['payload_rs_len'])
    config['wide_coder'] = RSWideAdapter(bits_per_z=bits_per_z, payload_len=config['oligos_per_block_len'],
                                         payload_rs_len=config['oligos_per_block_rs_len'])

    return config


config = build_config()
