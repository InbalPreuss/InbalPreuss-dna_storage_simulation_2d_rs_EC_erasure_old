import json
from pathlib import Path
from textwrap import wrap

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def load_data_to_df() -> pd.DataFrame:
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
    return df


def draw_boxplots(df: pd.DataFrame, percentage: bool = False):
    trials_group = df.groupby([
        'number_of_oligos_per_barcode',
        'number_of_sampled_oligos_from_file',
        'size',
        'bits_per_z'
    ])

    errors = ["substitution_error", "deletion_error", "insertion_error"]

    for idx, trial_group in trials_group:
        fig, axes = plt.subplots(nrows=3)
        fig.suptitle("\n".join(wrap(trial_group["output_dir"].iloc[0].split("[ errors")[0], 71)))
        fig.subplots_adjust(top=0.9)
        for ax, error in zip(axes, errors):
            zero_cols = [e for e in errors if e != error]
            df_for_err = trial_group
            for col in zero_cols:
                df_for_err = df_for_err[df_for_err[col] == 0]
            if percentage:
                sns.barplot(x=error, y="levenshtein_distance", data=df_for_err, estimator=estimate, ax=ax)
                ax.set(ylabel="D [levenshtein]")
            else:
                ax = sns.boxplot(x=error, y="levenshtein_distance", data=df_for_err, ax=ax)
                ax = sns.swarmplot(x=error, y="levenshtein_distance", data=df_for_err, color=".25", ax=ax)
                ax.set(ylabel="D [levenshtein]")


def draw_sampled_vs_error(df: pd.DataFrame):
    fig, ax = plt.subplots()
    ax = sns.boxplot(x="number_of_sampled_oligos_from_file", y="levenshtein_distance", data=df, ax=ax)
    ax = sns.swarmplot(x="number_of_sampled_oligos_from_file", y="levenshtein_distance", data=df, color=".25", ax=ax)
    ax.set(xlabel="number_of_sampled_oligos_from_file".replace("_", " "), ylabel="D [levenshtein]")


def estimate(x):
    res = sum(x == 0) * 100.0 / len(x)
    return res


def draw_zero_error_percentage(df: pd.DataFrame):
    fig, ax = plt.subplots()
    sns.barplot(x="number_of_sampled_oligos_from_file", y="levenshtein_distance", data=df, estimator=estimate, ax=ax)
    ax.set(xlabel="number_of_sampled_oligos_from_file".replace("_", " "), ylabel="D [levenshtein]")


def main():
    plt.ion()
    df = load_data_to_df()
    draw_zero_error_percentage(df=df)
    draw_boxplots(df=df)
    draw_boxplots(df=df, percentage=True)
    draw_sampled_vs_error(df=df)
    input("Hit enter to terminate")


if __name__ == '__main__':
    main()
