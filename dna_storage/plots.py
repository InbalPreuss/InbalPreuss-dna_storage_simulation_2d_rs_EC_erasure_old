import json
from pathlib import Path
from textwrap import wrap

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


output_dir = Path("data/testing")

run_dirs = list(f for f in output_dir.iterdir() if f.name.startswith("["))

json_files = []

for run_dir in run_dirs:
    try:
        json_file = [f for f in run_dir.iterdir() if f.suffix == ".json"][0]
    except IndexError:
        continue
    json_files.append(json_file)

info_list = []

for json_file in json_files:
    with json_file.open("rb") as f:
        info = json.load(f)
    info_list.append(info)

df = pd.DataFrame(info_list)
df["output_dir"] = df["output_dir"].apply(lambda x: x.split("data/testing/")[1].split(" trial")[0])

trials_group = df.groupby([
    'number_of_oligos_per_barcode',
    'number_of_sampled_oligos_from_file',
    'size',
    'bits_per_z'
])

for idx, trial_group in trials_group:
    fig, ax = plt.subplots()
    fig.suptitle("\n".join(wrap(trial_group["output_dir"].iloc[0].split("[ errors")[0], 71)))
    fig.subplots_adjust(top=0.8)
    ax = sns.boxplot(x="substitution_error", y="levenshtein_distance", data=trial_group, ax=ax)
    ax = sns.swarmplot(x="substitution_error", y="levenshtein_distance", data=trial_group, color=".25", ax=ax)


input("?")
