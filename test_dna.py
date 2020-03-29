import copy
import os
import pytest
import pathlib
import itertools
import pickle

import numpy as np
import matplotlib.pyplot as plt
import Levenshtein as levenshtein

from config import build_config
from main import main


def test_with_subset_size():
    errors = [0.01, 0.001, 0.0001, 0]
    results = {}
    sizes_and_bit_sizes = [(3, 9), (5, 12), (7, 13)]
    res_file = pathlib.Path(r'data/testing/output/error_results.pk')
    if not res_file.is_file():
        for size, bits_per_z in sizes_and_bit_sizes:
            products = itertools.product(errors, repeat=3)
            for prod in products:
                config = build_config(subset_size=size,
                                      bits_per_z=bits_per_z,
                                      letter_replace_error_ratio=prod[0],
                                      letter_remove_error_ratio=prod[1],
                                      letter_add_error_ratio=prod[2])
                dist, input_data, output_data = test_config(config)
                pos = size, *prod
                print(pos)
                results[pos] = {'dist': dist, 'input_data': input_data, 'output_data': output_data}

        with open(res_file, 'wb') as f:
            pickle.dump(results, f)

    else:
        with open(res_file, 'rb') as f:
            results = pickle.load(f)

    triples = [(errors, [0], [0]), ([0], errors, [0]), ([0], [0], errors)]
    for size, _ in sizes_and_bit_sizes:
        for idx, triple in enumerate(triples, 1):
            prod0 = [{key: val} for key, val in results.items() if
                     key[0] == size and
                     any(key[1] == v for v in triple[0]) and
                     any(key[2] == v for v in triple[1]) and
                     any(key[3] == v for v in triple[2])]
            x = [list(r.keys())[0][idx] for r in prod0]
            y = [list(r.values())[0]['dist'] for r in prod0]
            x, y = zip(*sorted(zip(x, y), reverse=True))

            fig, ax = plt.subplots()
            ax.plot(x, y)
            ax.set_xlabel('error rate')
            ax.set_ylabel('Levenshtein distance')
            ax.set_xscale('log')
            name = f'subset size {size} ' \
                   f'replace {1 if len(triple[0]) > 1 else 0} ' \
                   f'remove {1 if len(triple[1]) > 1 else 0} ' \
                   f'add {1 if len(triple[2]) > 1 else 0}'
            ax.title(name)
            plt.savefig(f'data/testing/output/{name}.png')
            fig = plt.gcf()
            plt.close(fig)


def test_config(config):
    with open('data/testing/input_text.dna', 'r', encoding='utf-8') as input_file:
        input_data = input_file.read()
    input_data = input_data.rsplit()
    main(config)
    with open('data/testing/small_data.text_results_file.dna', 'r', encoding='utf-8') as file:
        output_data = file.read()
    output_data = output_data.rsplit()
    input_data = ''.join(input_data)
    output_data = ''.join(output_data)
    dist = levenshtein.distance(input_data, output_data)
    return dist, input_data, output_data


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
    test_with_subset_size()
    # test_full_flow()
