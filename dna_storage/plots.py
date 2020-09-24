import json
import uuid
from pathlib import Path
from textwrap import wrap

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px


def load_sorted_oligos_to_df() -> pd.DataFrame:
    output_dir = Path("data/testing")

    run_dirs = list(f for f in output_dir.iterdir() if f.name.startswith("["))

    sorted_files = []

    for run_dir in run_dirs:
        try:
            sorted_file = [f for f in run_dir.iterdir() if "sort_oligo_results_file" in f.name][0]
        except IndexError:
            continue
        sorted_files.append(sorted_file)

    info_list = []

    for sorted_file in sorted_files:
        uid = uuid.uuid4()
        with sorted_file.open("rb") as f:
            for line in f:
                line = line.strip()
                info_list.append({
                    "uid": uid,
                    "barcode": line[:12].decode(),
                    "read_len": len(line[12:]),
                })

    df_reads = pd.DataFrame(info_list)
    return df_reads


def draw_reads_histograms(df: pd.DataFrame):
    plt.subplots()
    fig = sns.distplot(df["read_len"], kde=False)
    fig.set_ylabel("count")
    fig = plt.gcf()
    fig.savefig("read_len_count")

    groups = df.groupby(["uid", "barcode"])
    reads_per_barcode = []
    for _, group in groups:
        reads_per_barcode.append({"reads": len(group)})

    plt.subplots()
    df_reads_per_barcode = pd.DataFrame(reads_per_barcode)
    fig = sns.distplot(df_reads_per_barcode["reads"], kde=False)
    fig.set_ylabel("count")
    fig = plt.gcf()
    fig.savefig("reads per barcode")


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
        fig.subplots_adjust(top=0.85, hspace=0.5)
        for ax, error in zip(axes, errors):
            zero_cols = [e for e in errors if e != error]
            df_for_err = trial_group
            for col in zero_cols:
                df_for_err = df_for_err[df_for_err[col] == 0]
            if percentage:
                sns.barplot(x=error, y="levenshtein_distance", data=df_for_err, estimator=estimate, ax=ax, palette="Blues")
                ax.set(ylabel="D% [levenshtein]")
                fig.savefig(" ".join(wrap(trial_group["output_dir"].iloc[0].split("[ errors")[0].replace("[", " ").replace("]", " ") + "success", 71)))
                plt.close(fig)
            else:
                ax = sns.boxplot(x=error, y="levenshtein_distance", data=df_for_err, ax=ax, palette="Blues")
                ax = sns.swarmplot(x=error, y="levenshtein_distance", data=df_for_err, color=".25", ax=ax)
                ax.set(ylabel="D [levenshtein]")
                fig.savefig("".join(wrap(trial_group["output_dir"].iloc[0].split("[ errors")[0].replace("[", " ").replace("]", " "), 71)))
                plt.close(fig)


def draw_sampled_vs_error(df: pd.DataFrame):
    fig, ax = plt.subplots()
    ax = sns.boxplot(x="number_of_sampled_oligos_from_file", y="levenshtein_distance", data=df, ax=ax, palette="Blues")
    ax = sns.swarmplot(x="number_of_sampled_oligos_from_file", y="levenshtein_distance", data=df, color=".25", ax=ax)
    ax.set(xlabel="number_of_sampled_oligos_from_file".replace("_", " "), ylabel="D [levenshtein]")
    fig.savefig("number of sampled oligos and levenshtein distance")
    plt.close(fig)


def estimate(x):
    res = sum(x == 0) * 100.0 / len(x)
    return res


def draw_zero_error_percentage(df: pd.DataFrame):
    fig, ax = plt.subplots()
    sns.barplot(x="number_of_sampled_oligos_from_file", y="levenshtein_distance", data=df, estimator=estimate, ax=ax, palette="Blues")
    ax.set(xlabel="number_of_sampled_oligos_from_file".replace("_", " "), ylabel="D% [levenshtein]")
    fig.savefig("number of sampled oligos and levenshtein distance success")
    plt.close(fig)


def draw_error_per_number_of_sampled_oligos(df: pd.DataFrame):
    for idx, row in df.iterrows():
        if row["substitution_error"] != 0:
            df.loc[idx, "error_type"] = "substitution_error"
            df.loc[idx, "error"] = row["substitution_error"]
        elif row["deletion_error"] != 0:
            df.loc[idx, "error_type"] = "deletion_error"
            df.loc[idx, "error"] = row["deletion_error"]
        elif row["insertion_error"] != 0:
            df.loc[idx, "error_type"] = "insertion_error"
            df.loc[idx, "error"] = row["insertion_error"]
        else:
            df.loc[idx, "error_type"] = "no_error"
            df.loc[idx, "error"] = 0

    ax = sns.catplot(
        data=df,
        x="number_of_sampled_oligos_from_file",
        y="levenshtein_distance",
        row="error_type",
        col="error",
        palette="Blues",
        kind="bar",
    )

    fig = px.bar(
        df, x="number_of_sampled_oligos_from_file", y="levenshtein_distance",
        facet_row="error_type", facet_col="error",
        barmode="group",
    )
    fig.show()

def main():
    plt.ion()
    df_reads = load_sorted_oligos_to_df()
    draw_reads_histograms(df=df_reads)
    
    df = load_data_to_df()
    draw_zero_error_percentage(df=df)
    draw_boxplots(df=df)
    draw_boxplots(df=df, percentage=True)
    draw_sampled_vs_error(df=df)
    draw_error_per_number_of_sampled_oligos(df=df)
    input("Hit enter to terminate")


if __name__ == '__main__':
    main()
