import datetime
import os
import itertools
import pathlib
import gzip
import os
from typing import Dict
import json
from pathlib import Path

import Levenshtein as levenshtein

from dna_storage.config import build_config
from dna_storage.main import main
from dna_storage.text_handling import generate_random_text_file


def build_runs():
    number_of_oligos_per_barcode = [1000]
    # number_of_oligos_per_barcode = [20, 100, 1000, 10000]`
    number_of_sampled_oligos_from_file = [-1, 20, 50, 100, 1000]
    # number_of_sampled_oligos_from_file = [20, 50, 100, 1000, float('inf')]
    oligos_and_samples = list(itertools.product(number_of_oligos_per_barcode, number_of_sampled_oligos_from_file))
    oligos_and_samples = [s for s in oligos_and_samples if s[0] >= s[1]]

    errors = [0.01, 0.001, 0.0001, 0]
    sizes_and_bit_sizes = [(3, 9), (5, 12), (7, 13)]
    variable_number_of_sampled_oligos_from_file = {3: 5, 5: 10, 7: 15}

    runs = []
    was_variable = False
    for number_of_oligos_per_barcode, number_of_sampled_oligos_from_file in oligos_and_samples:
        for size, bits_per_z in sizes_and_bit_sizes:
            prods = list(itertools.product(errors, repeat=3))
            products = [i for i in prods if
                        (i[0] == 0 and i[1] == 0) or (i[0] == 0 and i[2] == 0) or (i[1] == 0 and i[2] == 0)]
            for prod in products:
                if number_of_sampled_oligos_from_file == -1:
                    was_variable = True
                    number_of_sampled_oligos_from_file = variable_number_of_sampled_oligos_from_file[size]
                name = f'[ subset size {size}, bits per z {bits_per_z:>2} ]' \
                       f'[ number of oligos per barcode {number_of_oligos_per_barcode:>6} ]\n' \
                       f'[ number of oligos sampled after synthesis {number_of_sampled_oligos_from_file:>6} ]\n' \
                       f'[ errors, substitution {prod[0]:<6}, deletion {prod[1]:<6}, insertion {prod[2]:<6} ]\n'
                output_dir = os.path.join("data/testing/", name.replace("\n", ""))
                runs.append({
                    "number_of_oligos_per_barcode": number_of_oligos_per_barcode,
                    "number_of_sampled_oligos_from_file": number_of_sampled_oligos_from_file,
                    "size": size,
                    "bits_per_z": bits_per_z,
                    "substitution_error": prod[0],
                    "deletion_error": prod[1],
                    "insertion_error": prod[2],
                    "output_dir": output_dir
                })
                if was_variable:
                    number_of_sampled_oligos_from_file = -1
                    was_variable = False
    return runs


def run_config_n_times(config_for_run: Dict, n: int = 10):
    for run_number in range(n):
        run_config(config_for_run=config_for_run, run_number=run_number)


def run_config(config_for_run: Dict, run_number):
    output_dir = f"{config_for_run['output_dir']} trial {run_number:>2}"
    input_text = output_dir + '/random_file_10_KiB.txt'
    config = build_config(
        subset_size=config_for_run["size"],
        bits_per_z=config_for_run["bits_per_z"],
        letter_substitution_error_ratio=config_for_run["substitution_error"],
        letter_deletion_error_ratio=config_for_run["deletion_error"],
        letter_insertion_error_ratio=config_for_run["insertion_error"],
        number_of_oligos_per_barcode=config_for_run["number_of_oligos_per_barcode"],
        number_of_sampled_oligos_from_file=config_for_run["number_of_sampled_oligos_from_file"],
        output_dir=output_dir,
        input_text_file=input_text
    )

    generate_random_text_file(size_kb=0.1, file=input_text)
    print(f"$$$$$$$$ Running {output_dir} $$$$$$$$")
    main(config)
    with open(input_text, 'r', encoding='utf-8') as input_file:
        input_data = input_file.read()
    with open(Path(output_dir) / 'simulation_data.9.text_results_file.dna', 'r',
              encoding='utf-8') as file:
        output_data = file.read()

    # gzip and delete files
    files_in_dir = list(Path(output_dir).iterdir())
    exclude = ["temp_shuffle_db", "temp_sort_oligo_db"]
    files = [f for f in files_in_dir if f.name not in exclude]
    files_delete = [f for f in files_in_dir if f.name in exclude]
    for file in files:
        with open(file, "rb") as f_in, gzip.open(file.with_suffix(file.suffix + ".gz"), "wb") as f_out:
            f_out.writelines(f_in)
        os.remove(file)
    [os.remove(f) for f in files_delete]

    # write a json results file
    dist = levenshtein.distance(input_data, output_data)
    res_file = Path(output_dir) / f"config_and_levenshtein_distance_{dist}.json"
    res = {**config_for_run, "output_dir": output_dir, "levenshtein_distance": dist}
    with open(res_file, 'w') as f:
        json.dump(res, f, indent=4)
    print(f"@@@@@@@@ Finished {output_dir} @@@@@@@@")
    finished_dir = Path('data/finished')
    finished_dir.mkdir(parents=True, exist_ok=True)
    with open(finished_dir / Path(output_dir).name, "w") as f:
        f.write("1")
    return dist, input_data, output_data


def main_fn():
    from multiprocessing import Pool, cpu_count
    configs_for_run = build_runs()
    with Pool(cpu_count()) as p:
        p.map(run_config_n_times, configs_for_run)


if __name__ == '__main__':
    main_fn()
