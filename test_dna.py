import copy
import os
import pytest
import pathlib

from main import main


# def test_small_data_barcode_4_oligo_8():
#     from config import config
#     config = copy.deepcopy(config)
#     config['NUMBER_OF_BARCODE_LETTERS'] = 4
#     config['OLIGO_LENGTH'] = 8
#     config['K_MER'] = 1
#     config['file_name_sorted'] = pathlib.Path(r'data/testing/small_data_4_barcode_8_oligo.dna')
#     config['unique_oligo_results_file'] = pathlib.Path(
#         r'data/testing/small_data_4_barcode_8_oligo.unique_oligo_results_file.dna')
#     try:
#         os.remove('data/testing/small_data_4_barcode_8_oligo.unique_oligo_results_file.dna')
#     except OSError:
#         pass
#     main(config)
#     with open('data/testing/small_data_4_barcode_8_oligo.unique_oligo_results_file.dna', 'r') as file:
#         data = file.read()
#     assert data == """AAAAACGTACGT
# ACGTAAAAAAAA
# ACGGAAAAAAAA
# ACACCCCCCCCC
# """


# def test_small_data_barcode_3_oligo_9():
#     from config import config
#     config = copy.deepcopy(config)
#     try:
#         os.remove('data/testing/small_data_3_barcode_9_oligo.unique_oligo_results_file.dna')
#     except OSError:
#         pass
#     main(config)
#     with open('data/testing/small_data_3_barcode_9_oligo.unique_oligo_results_file.dna', 'r') as file:
#         data = file.read()
#     assert data == """AAAX1X2X1
# AGCX7X7X7
# AGTX7X7X7
# ACAX14X14X14
# """

def test_full_flow():
    from config import config
    config = copy.deepcopy(config)
    config['synthesis']['letter_replace_error_ratio'] = 0
    config['synthesis']['letter_remove_error_ratio'] = 0
    config['synthesis']['letter_add_error_ratio'] = 0
    input_data = ''
    with open('data/testing/small_data_binary', 'r') as input_file:
        input_data = input_file.read()
    input_data = input_data.rsplit()
    main(config)
    with open('data/testing/small_data_binary.decoder_results_file.dna', 'r') as file:
        data = file.read()
    data = data.rsplit()
    data = [d[3:] for d in data]
    bpz = config['algorithm_config']['bits_per_z']
    result = []
    for d in data:
        res = [d[i:i+bpz] for i in range(0, len(d), bpz)]
        result.extend(res)
    assert input_data == result


if __name__ == '__main__':
    test_full_flow()
