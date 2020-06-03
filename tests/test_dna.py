import copy
from typing import Union
import pathlib
import itertools
import pickle

import matplotlib.pyplot as plt
import Levenshtein as levenshtein

from dna_storage.config import build_config
from dna_storage.main import main


def test_number_of_oligos_per_barcode():
    for number_of_oligos_per_barcode in [20, 100, 1000, 10000]:
        for number_of_sampled_oligos_from_file in [20, 50, 100, None]:
            subset_size_and_error_plot(number_of_oligos_per_barcode=number_of_oligos_per_barcode,
                                       number_of_sampled_oligos_from_file=number_of_sampled_oligos_from_file)


def subset_size_and_error_plot(number_of_oligos_per_barcode: int = 20,
                               number_of_sampled_oligos_from_file: Union[int, None] = 10000):
    plt.ion()
    errors = [0.01, 0.001, 0.0001, 0]
    results = {}
    sizes_and_bit_sizes = [(3, 9), (5, 12), (7, 13)]
    res_file = pathlib.Path(f'data/testing/output/error_results_n_oligos_{number_of_oligos_per_barcode}_n_sampled_{number_of_sampled_oligos_from_file}.pk')
    if not res_file.is_file():
        for size, bits_per_z in sizes_and_bit_sizes:
            # products = itertools.product(errors, repeat=3)
            prods = list(itertools.product(errors, repeat=3))
            products = [i for i in prods if (i[0] == 0 and i[1] == 0) or (i[0] == 0 and i[2] == 0) or (i[1] == 0 and i[2] == 0)]
            for prod in products:
                config = build_config(subset_size=size,
                                      bits_per_z=bits_per_z,
                                      letter_replace_error_ratio=prod[0],
                                      letter_remove_error_ratio=prod[1],
                                      letter_add_error_ratio=prod[2],
                                      number_of_oligos_per_barcode=number_of_oligos_per_barcode,
                                      number_of_sampled_oligos_from_file=number_of_sampled_oligos_from_file)
                dist, input_data, output_data = run_pipe_with_config(config)
                pos = size, *prod
                if number_of_sampled_oligos_from_file is None:
                    n_samples = 'all'
                else:
                    n_samples = number_of_sampled_oligos_from_file
                print(f"[ size: {size:>1} ]",
                      f"[ errors: (replace: {prod[0]:<6}), (remove: {prod[1]:<6}), (add: {prod[2]:<6}) ]",
                      f"[ dist: {dist:>5}]",
                      f"[ n_oligos_per_barcode: {number_of_oligos_per_barcode:>6} ]",
                      f"[ m_samples_from_synthesis_file: {n_samples:>6}]")
                print(f"[ input data: {input_data} ]",
                      f"[ output data: {output_data}]")
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
            ax.plot(x, y, marker='o', linestyle='-')
            ax.set_xlabel('error rate')
            ax.set_ylabel('Levenshtein distance')
            ax.set_xscale('symlog', linthreshx=errors[-2])
            ax.set_xticks(errors)
            ax.set_xlim([0, errors[0]])

            name = f'[ subset size {size} ]' \
                   f'[ number of oligos per barcode {number_of_oligos_per_barcode} ]\n' \
                   f'[ number of oligos sampled after synthesis {number_of_sampled_oligos_from_file} ]\n' \
                   f'[ replace {1 if len(triple[0]) > 1 else 0} ]' \
                   f'[ remove {1 if len(triple[1]) > 1 else 0} ]' \
                   f'[ add {1 if len(triple[2]) > 1 else 0} ]'
            ax.set_title(name)
            fig.tight_layout(rect=[0, 0, 1, 1])
            name = name.replace("\n", "")
            fig.savefig(f'data/testing/output/{name}.png')
            plt.close(fig)


def timing():
    import time
    from dna_storage.config import build_config
    from dna_storage.text_handling import generate_random_text_file
    config = build_config(input_text_file=pathlib.Path(r'data/testing/random_file_1_KiB.txt'))
    print('Generating 1 KiB file')
    generate_random_text_file(size_kb=1, file='./data/testing/random_file_1_KiB.txt')
    t = time.time()
    print('Running 1 KiB file')
    main(config)
    kb_time = time.time() - t
    print(f"1 KiB took {kb_time} seconds")

    config = build_config(input_text_file=pathlib.Path(r'data/testing/random_file_1_MiB.txt'))
    print('\nGenerating 1 MiB file')
    generate_random_text_file(size_kb=16, file='./data/testing/random_file_1_MiB.txt')
    t = time.time()
    print('Running 1 MiB file')
    main(config)
    mb_time = time.time() - t
    print(f"1 MiB took {mb_time} seconds")

    with open('./data/testing/timings.txt', 'w') as f:
        f.write(f"Timings:\n========\n1 KiB time: {kb_time}\n1 MiB time: {mb_time}\n")


def code_profiling(size_kb: int = 1):
    from dna_storage.text_handling import generate_random_text_file
    file = pathlib.Path(f'data/testing/random_file_{size_kb}_KiB.txt')
    generate_random_text_file(size_kb=size_kb, file=file)
    config = build_config(input_text_file=file)
    main(config)


def run_pipe_with_config(config):
    with open('./data/testing/input_text.dna', 'r', encoding='utf-8') as input_file:
        input_data = input_file.read()

    main(config)
    with open('./data/testing/simulation_data.9.text_results_file.dna', 'r', encoding='utf-8') as file:
        output_data = file.read()

    dist = levenshtein.distance(input_data, output_data)
    return dist, input_data, output_data


def test_full_flow():
    from dna_storage.config import build_config
    config = build_config(bits_per_z=13, subset_size=7)
    with open('./data/testing/input_text.dna', 'r', encoding='utf-8') as input_file:
        input_data = input_file.read()
    main(config)
    with open('./data/testing/simulation_data.9.text_results_file.dna', 'r', encoding='utf-8') as file:
        data = file.read()
    assert input_data == data


if __name__ == '__main__':
    test_number_of_oligos_per_barcode()
    # code_profiling(size_kb=1024)
    # timing()
    # subset_size_and_error_plot()
    # test_full_flow()
