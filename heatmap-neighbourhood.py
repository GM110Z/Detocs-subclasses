#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def load_table(filepath):
    ext = os.path.splitext(filepath)[1].lower()

    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(filepath)
    if ext == ".csv":
        return pd.read_csv(filepath)
    if ext in [".tsv", ".txt"]:
        return pd.read_csv(filepath, sep="\t")

    raise ValueError(f"Unsupported input format: {ext}")


def parse_combined_column(series, label):
    series = series.dropna().astype(str).str.strip()

    parsed = series.str.extract(r"^(.*?)\s+(\d+)$")
    parsed.columns = ["defence", "count"]

    bad = series[parsed["defence"].isna()]
    if len(bad) > 0:
        raise ValueError(
            f"Parsing error in {label}. Expected 'Name 10'.\n"
            + "\n".join(bad.tolist())
        )

    parsed["count"] = parsed["count"].astype(int)
    parsed["focal_system"] = label
    return parsed


def normalise_input(df, left_name="Detocs", right_name="M60"):
    cols = list(df.columns)

    # 4-column format
    if {left_name, f"{left_name}_count", right_name, f"{right_name}_count"}.issubset(cols):
        left = df[[left_name, f"{left_name}_count"]]
        left.columns = ["defence", "count"]
        left = left.dropna(subset=["defence"])
        left["focal_system"] = left_name

        right = df[[right_name, f"{right_name}_count"]]
        right.columns = ["defence", "count"]
        right = right.dropna(subset=["defence"])
        right["focal_system"] = right_name

        return pd.concat([left, right], ignore_index=True)

    # 2-column format
    if {left_name, right_name}.issubset(cols):
        left = parse_combined_column(df[left_name], left_name)
        right = parse_combined_column(df[right_name], right_name)
        return pd.concat([left, right], ignore_index=True)

    raise ValueError("Input format not recognised")


def make_matrix(long_df, left_name="Detocs", right_name="M60"):
    mat = long_df.pivot_table(
        index="defence",
        columns="focal_system",
        values="count",
        aggfunc="sum",
        fill_value=0
    )

    for col in [left_name, right_name]:
        if col not in mat.columns:
            mat[col] = 0

    return mat[[left_name, right_name]]


def to_percentage(mat, n_left=None, n_right=None):
    """
    Convert to % of loci containing that system
    """
    if n_left:
        mat.iloc[:, 0] = mat.iloc[:, 0] / n_left * 100
    if n_right:
        mat.iloc[:, 1] = mat.iloc[:, 1] / n_right * 100
    return mat


def zscore_rows(mat):
    """
    Row-wise z-score → highlights relative enrichment
    """
    return mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, 1), axis=0)


def plot_heatmap(mat, out_png, title=None):
    plt.figure(figsize=(4, max(4, len(mat) * 0.4)))

    vmax = abs(mat.values).max()

    im = plt.imshow(
        mat,
        aspect="auto",
        cmap="Pastel2",   # blue = low, red = high
        vmin=-vmax,
        vmax=vmax
    )

    plt.xticks(range(len(mat.columns)), mat.columns)
    plt.yticks(range(len(mat.index)), mat.index)

    plt.colorbar(label="Z-score (relative enrichment)")

    if title:
        plt.title(title)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Heatmap of neighbourhood composition")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--outprefix", default="heatmap")
    parser.add_argument("--left-name", default="Detocs")
    parser.add_argument("--right-name", default="M60")

    # optional normalisation
    parser.add_argument("--n-left", type=float, help="Total number of Detocs loci")
    parser.add_argument("--n-right", type=float, help="Total number of M60 loci")

    parser.add_argument(
        "--mode",
        choices=["raw", "percent", "zscore"],
        default="zscore",
        help="How to scale values"
    )

    parser.add_argument("--title", default=None)

    args = parser.parse_args()

    raw = load_table(args.input)
    long_df = normalise_input(raw, args.left_name, args.right_name)
    mat = make_matrix(long_df, args.left_name, args.right_name)

    # apply scaling
    if args.mode == "percent":
        mat = to_percentage(mat, args.n_left, args.n_right)

    elif args.mode == "zscore":
        if args.n_left and args.n_right:
            mat = to_percentage(mat, args.n_left, args.n_right)
        mat = zscore_rows(mat)

    # save table for Prism
    mat_out = f"{args.outprefix}.heatmap_table.csv"
    mat.to_csv(mat_out)

    # plot
    png_out = f"{args.outprefix}.heatmap.png"
    plot_heatmap(mat, png_out, args.title)

    print(f"Saved table: {mat_out}")
    print(f"Saved heatmap: {png_out}")


if __name__ == "__main__":
    main()
