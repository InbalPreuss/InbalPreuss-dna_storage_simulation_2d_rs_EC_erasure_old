import os
import itertools
from typing import Dict
import json
from pathlib import Path

import Levenshtein as levenshtein

from dna_storage.config import build_config
from dna_storage.main import main


def build_runs():
    number_of_oligos_per_barcode = [20, 100, 1000, 10000]
    number_of_sampled_oligos_from_file = [20, 50, 100, 1000, float('inf')]
    oligos_and_samples = list(itertools.product(number_of_oligos_per_barcode, number_of_sampled_oligos_from_file))
    oligos_and_samples = [s for s in oligos_and_samples if s[0] <= s[1]]

    errors = [0.01, 0.001, 0.0001, 0]
    sizes_and_bit_sizes = [(3, 9), (5, 12), (7, 13)]

    runs = []
    for number_of_oligos_per_barcode, number_of_sampled_oligos_from_file in oligos_and_samples:
        for size, bits_per_z in sizes_and_bit_sizes:
            prods = list(itertools.product(errors, repeat=3))
            products = [i for i in prods if
                        (i[0] == 0 and i[1] == 0) or (i[0] == 0 and i[2] == 0) or (i[1] == 0 and i[2] == 0)]
            for prod in products:
                name = f'[ subset size {size}, bits per z {bits_per_z:>2} ]' \
                       f'[ number of oligos per barcode {number_of_oligos_per_barcode:>6} ]\n' \
                       f'[ number of oligos sampled after synthesis {number_of_sampled_oligos_from_file:>6} ]\n' \
                       f'[ errors, replace {prod[0]:<6}, remove {prod[1]:<6}, add {prod[2]:<6} ]\n'
                output_dir = os.path.join("data/testing/", name.replace("\n", ""))
                runs.append({
                    "number_of_oligos_per_barcode": number_of_oligos_per_barcode,
                    "number_of_sampled_oligos_from_file": number_of_sampled_oligos_from_file,
                    "size": size,
                    "bits_per_z": bits_per_z,
                    "replace_error": prod[0],
                    "remove_error": prod[1],
                    "add_error": prod[2],
                    "output_dir": output_dir
                })
    return runs


def run_config(config_for_run: Dict):
    config = build_config(
        subset_size=config_for_run["size"],
        bits_per_z=config_for_run["bits_per_z"],
        letter_replace_error_ratio=config_for_run["replace_error"],
        letter_remove_error_ratio=config_for_run["remove_error"],
        letter_add_error_ratio=config_for_run["add_error"],
        number_of_oligos_per_barcode=config_for_run["number_of_oligos_per_barcode"],
        number_of_sampled_oligos_from_file=config_for_run["number_of_sampled_oligos_from_file"],
        output_dir=config_for_run["output_dir"]
    )

    print(f"$$$$$$$$ Running {config_for_run['output_dir']} $$$$$$$$")
    main(config)
    with open('./data/testing/input_text.dna', 'r', encoding='utf-8') as input_file:
        input_data = input_file.read()
    with open(Path(config_for_run["output_dir"]) / 'simulation_data.9.text_results_file.dna', 'r', encoding='utf-8') as file:
        output_data = file.read()

    dist = levenshtein.distance(input_data, output_data)
    res_file = Path(config_for_run["output_dir"]) / f"config_and_levenshtein_distance_{dist}.json"
    res = {**config_for_run, "levenshtein_distance": dist}
    with open(res_file, 'w') as f:
        json.dump(res, f, indent=4)
    print(f"@@@@@@@@ Finished {config_for_run['output_dir']} @@@@@@@@")
    finished_dir = Path('data/finished')
    finished_dir.mkdir(parents=True, exist_ok=True)
    with open(finished_dir / Path(config_for_run['output_dir']).name, "w") as f:
        f.write("1")
    return dist, input_data, output_data


def main_fn():
    from multiprocessing import Pool, cpu_count
    configs_for_run = build_runs()
    with Pool(cpu_count()) as p:
        p.map(run_config, configs_for_run)


if __name__ == '__main__':
    main_fn()
