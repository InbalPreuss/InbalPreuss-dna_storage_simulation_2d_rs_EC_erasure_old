import copy
import os
import pytest
import pathlib

from main import main


def test_full_flow():
    from config import config
    config = copy.deepcopy(config)
    config['synthesis']['letter_replace_error_ratio'] = 0
    config['synthesis']['letter_remove_error_ratio'] = 0
    config['synthesis']['letter_add_error_ratio'] = 0
    with open('data/testing/input_text.dna', 'r', encoding='utf-8') as input_file:
        input_data = input_file.read()
    input_data = input_data.rsplit()
    main(config)
    with open('data/testing/small_data.text_results_file.dna', 'r', encoding='utf-8') as file:
        data = file.read()
    data = data.rsplit()
    assert input_data == data


if __name__ == '__main__':
    test_full_flow()
