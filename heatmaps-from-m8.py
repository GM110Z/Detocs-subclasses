#!/usr/bin/env python3
"""
Dereplicate two sequence families and generate three pairwise-similarity heatmaps
from an all-vs-all m8-like table.

Expected input columns:
    1 = query
    2 = target
    3 = identity (0-1 or 0-100)
    4 = alignment length
    5 = e-value

Outputs:
- representatives for family A and family B
- 3 TSV matrices
- 3 PNG heatmaps:
    * A_vs_A
    * B_vs_B
    * A_vs_B

Example:
python family_heatmaps_from_m8.py \
    --m8 all_vs_all.m8 \
    --prefix-a Cap17_ \
    --prefix-b M60C_ \
    --outdir heatmaps_cap17_m60c \
    --derep-id 0.95 \
    --min-id 0.30 \
    --min-aln 50 \
    --max-evalue 1e-3
"""

from __future__ import annotations

import argparse
import os
from collections import deque
from typing import Dict, List, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Dereplicate two families and generate three pairwise heatmaps from an all-vs-all m8 file."
    )
    parser.add_argument("--m8", required=True, help="Input all-vs-all m8 file")
    parser.add_argument("--prefix-a", required=True, help="Prefix identifying family A, e.g. Cap17_")
    parser.add_argument("--prefix-b", required=True, help="Prefix identifying family B, e.g. M60C_")
    parser.add_argument(
        "--label-a",
        default=None,
        help="Optional display label for family A. If omitted, derived from prefix-a.",
    )
    parser.add_argument(
        "--label-b",
        default=None,
        help="Optional display label for family B. If omitted, derived from prefix-b.",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--derep-id",
        type=float,
        default=0.95,
        help="Dereplication identity threshold. Accepts 0-1 or 0-100.",
    )
    parser.add_argument(
        "--min-id",
        type=float,
        default=0.0,
        help="Minimum identity for heatmap matrix entries. Accepts 0-1 or 0-100.",
    )
    parser.add_argument("--min-aln", type=int, default=0, help="Minimum alignment length")
    parser.add_argument("--max-evalue", type=float, default=1.0, help="Maximum e-value")
    parser.add_argument(
        "--max-reps-per-family",
        type=int,
        default=0,
        help="Optional cap on number of dereplicated representatives per family. 0 = no cap.",
    )
    parser.add_argument(
        "--cross-fill",
        choices=["zero", "nan"],
        default="nan",
        help="How to fill missing A-vs-B values in the rectangular matrix",
    )
    return parser.parse_args()


def normalize_identity_value(x: float) -> float:
    if x > 1.0:
        return x / 100.0
    return x


def clean_label(prefix: str) -> str:
    return prefix.rstrip("_")


def read_m8(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", header=None, comment="#", engine="python")
    if df.shape[1] < 5:
        raise ValueError("Input file must have at least 5 columns: query target identity aln_len evalue")

    df = df.iloc[:, :5].copy()
    df.columns = ["query", "target", "identity", "aln_len", "evalue"]

    df["identity"] = pd.to_numeric(df["identity"], errors="coerce")
    df["aln_len"] = pd.to_numeric(df["aln_len"], errors="coerce")
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df = df.dropna(subset=["query", "target", "identity", "aln_len", "evalue"]).copy()

    df["identity"] = df["identity"].map(normalize_identity_value)
    return df


def filter_hits(df: pd.DataFrame, min_id: float, min_aln: int, max_evalue: float) -> pd.DataFrame:
    min_id = normalize_identity_value(min_id)
    keep = (
        (df["query"] != df["target"]) &
        (df["identity"] >= min_id) &
        (df["aln_len"] >= min_aln) &
        (df["evalue"] <= max_evalue)
    )
    return df.loc[keep].copy()


def collect_family_ids(df: pd.DataFrame, prefix: str) -> Set[str]:
    ids = set(df.loc[df["query"].str.startswith(prefix), "query"])
    ids |= set(df.loc[df["target"].str.startswith(prefix), "target"])
    return ids


def build_best_identity_map(df: pd.DataFrame) -> Dict[Tuple[str, str], float]:
    best: Dict[Tuple[str, str], float] = {}
    for row in df.itertuples(index=False):
        a = row.query
        b = row.target
        if a == b:
            continue
        key = tuple(sorted((a, b)))
        val = float(row.identity)
        if key not in best or val > best[key]:
            best[key] = val
    return best


def build_family_graph(
    ids: Set[str],
    best_map: Dict[Tuple[str, str], float],
    derep_id: float,
) -> Dict[str, Set[str]]:
    derep_id = normalize_identity_value(derep_id)
    graph: Dict[str, Set[str]] = {x: set() for x in ids}
    for (a, b), ident in best_map.items():
        if a in ids and b in ids and ident >= derep_id:
            graph[a].add(b)
            graph[b].add(a)
    return graph


def connected_components(graph: Dict[str, Set[str]]) -> List[List[str]]:
    seen: Set[str] = set()
    comps: List[List[str]] = []

    for node in sorted(graph):
        if node in seen:
            continue
        comp = []
        dq = deque([node])
        seen.add(node)
        while dq:
            cur = dq.popleft()
            comp.append(cur)
            for nbr in graph[cur]:
                if nbr not in seen:
                    seen.add(nbr)
                    dq.append(nbr)
        comps.append(sorted(comp))
    return comps


def choose_representatives(components: List[List[str]], max_reps: int = 0) -> List[str]:
    reps = [sorted(comp)[0] for comp in components]
    reps = sorted(reps)
    if max_reps > 0:
        reps = reps[:max_reps]
    return reps


def save_representatives(path: str, reps: List[str], label: str, components: List[List[str]]) -> None:
    rep_set = set(reps)
    rows = []
    cluster_num = 1
    for comp in components:
        rep = sorted(comp)[0]
        if rep not in rep_set:
            continue
        for member in comp:
            rows.append({
                "family": label,
                "cluster": f"{label}_cluster_{cluster_num}",
                "representative": rep,
                "member": member,
                "cluster_size": len(comp),
            })
        cluster_num += 1
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def make_matrix(
    rows: List[str],
    cols: List[str],
    best_map: Dict[Tuple[str, str], float],
    fill_mode: str = "nan",
) -> pd.DataFrame:
    fill_value = np.nan if fill_mode == "nan" else 0.0
    mat = pd.DataFrame(fill_value, index=rows, columns=cols, dtype=float)

    for r in rows:
        for c in cols:
            if r == c:
                mat.loc[r, c] = 1.0
                continue
            key = tuple(sorted((r, c)))
            if key in best_map:
                mat.loc[r, c] = best_map[key]

    return mat


def sort_by_mean_similarity_square(mat: pd.DataFrame) -> pd.DataFrame:
    means = mat.replace(0, np.nan).mean(axis=1, skipna=True).fillna(0)
    order = means.sort_values(ascending=False).index.tolist()
    return mat.loc[order, order]


def sort_cross_matrix(mat: pd.DataFrame) -> pd.DataFrame:
    row_means = mat.mean(axis=1, skipna=True).fillna(0)
    col_means = mat.mean(axis=0, skipna=True).fillna(0)
    row_order = row_means.sort_values(ascending=False).index.tolist()
    col_order = col_means.sort_values(ascending=False).index.tolist()
    return mat.loc[row_order, col_order]


def plot_heatmap(
    mat: pd.DataFrame,
    title: str,
    outfile: str,
    cmap: str = "inferno",
    vmin: float = 0.0,
    vmax: float = 1.0,
) -> None:
    arr = mat.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(arr)

    nrows, ncols = arr.shape
    fig_w = min(max(6, ncols * 0.18), 20)
    fig_h = min(max(6, nrows * 0.18), 20)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(masked, aspect="auto", interpolation="nearest", cmap=cmap, vmin=vmin, vmax=vmax)

    ax.set_title(title, fontsize=12)
    ax.set_xlabel("Target")
    ax.set_ylabel("Query")

    if ncols <= 80:
        ax.set_xticks(np.arange(ncols))
        ax.set_xticklabels(mat.columns, rotation=90, fontsize=6)
    else:
        ax.set_xticks([])

    if nrows <= 80:
        ax.set_yticks(np.arange(nrows))
        ax.set_yticklabels(mat.index, fontsize=6)
    else:
        ax.set_yticks([])

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Sequence identity")

    fig.tight_layout()
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    label_a = args.label_a if args.label_a else clean_label(args.prefix_a)
    label_b = args.label_b if args.label_b else clean_label(args.prefix_b)

    df = read_m8(args.m8)

    ids_a = collect_family_ids(df, args.prefix_a)
    ids_b = collect_family_ids(df, args.prefix_b)

    if not ids_a:
        raise ValueError(f"No IDs found with prefix-a: {args.prefix_a}")
    if not ids_b:
        raise ValueError(f"No IDs found with prefix-b: {args.prefix_b}")

    derep_df = filter_hits(df, min_id=args.derep_id, min_aln=args.min_aln, max_evalue=args.max_evalue)
    derep_best = build_best_identity_map(derep_df)

    graph_a = build_family_graph(ids_a, derep_best, args.derep_id)
    graph_b = build_family_graph(ids_b, derep_best, args.derep_id)

    comps_a = connected_components(graph_a)
    comps_b = connected_components(graph_b)

    reps_a = choose_representatives(comps_a, max_reps=args.max_reps_per_family)
    reps_b = choose_representatives(comps_b, max_reps=args.max_reps_per_family)

    save_representatives(
        os.path.join(args.outdir, f"{label_a}_representatives.tsv"),
        reps_a,
        label_a,
        comps_a,
    )
    save_representatives(
        os.path.join(args.outdir, f"{label_b}_representatives.tsv"),
        reps_b,
        label_b,
        comps_b,
    )

    plot_df = filter_hits(df, min_id=args.min_id, min_aln=args.min_aln, max_evalue=args.max_evalue)
    plot_best = build_best_identity_map(plot_df)

    mat_aa = make_matrix(reps_a, reps_a, plot_best, fill_mode="nan")
    mat_bb = make_matrix(reps_b, reps_b, plot_best, fill_mode="nan")
    mat_ab = make_matrix(reps_a, reps_b, plot_best, fill_mode=args.cross_fill)

    mat_aa = sort_by_mean_similarity_square(mat_aa)
    mat_bb = sort_by_mean_similarity_square(mat_bb)
    mat_ab = sort_cross_matrix(mat_ab)

    mat_aa.to_csv(os.path.join(args.outdir, f"{label_a}_vs_{label_a}_matrix.tsv"), sep="\t")
    mat_bb.to_csv(os.path.join(args.outdir, f"{label_b}_vs_{label_b}_matrix.tsv"), sep="\t")
    mat_ab.to_csv(os.path.join(args.outdir, f"{label_a}_vs_{label_b}_matrix.tsv"), sep="\t")

    plot_heatmap(
        mat_aa,
        title=f"{label_a} vs {label_a} (dereplicated)",
        outfile=os.path.join(args.outdir, f"{label_a}_vs_{label_a}_heatmap.png"),
    )
    plot_heatmap(
        mat_bb,
        title=f"{label_b} vs {label_b} (dereplicated)",
        outfile=os.path.join(args.outdir, f"{label_b}_vs_{label_b}_heatmap.png"),
    )
    plot_heatmap(
        mat_ab,
        title=f"{label_a} vs {label_b} (dereplicated)",
        outfile=os.path.join(args.outdir, f"{label_a}_vs_{label_b}_heatmap.png"),
    )

    with open(os.path.join(args.outdir, "run_summary.txt"), "w") as fh:
        fh.write(f"Input file: {args.m8}\n")
        fh.write(f"Prefix A: {args.prefix_a}\n")
        fh.write(f"Prefix B: {args.prefix_b}\n")
        fh.write(f"Label A: {label_a}\n")
        fh.write(f"Label B: {label_b}\n")
        fh.write(f"Dereplication identity: {normalize_identity_value(args.derep_id):.3f}\n")
        fh.write(f"Plot minimum identity: {normalize_identity_value(args.min_id):.3f}\n")
        fh.write(f"Minimum alignment length: {args.min_aln}\n")
        fh.write(f"Maximum e-value: {args.max_evalue}\n")
        fh.write(f"{label_a} IDs before dereplication: {len(ids_a)}\n")
        fh.write(f"{label_b} IDs before dereplication: {len(ids_b)}\n")
        fh.write(f"{label_a} representatives after dereplication: {len(reps_a)}\n")
        fh.write(f"{label_b} representatives after dereplication: {len(reps_b)}\n")

    print("Done.")
    print(f"{label_a} IDs before dereplication: {len(ids_a)}")
    print(f"{label_b} IDs before dereplication: {len(ids_b)}")
    print(f"{label_a} representatives after dereplication: {len(reps_a)}")
    print(f"{label_b} representatives after dereplication: {len(reps_b)}")
    print(f"Outputs written to: {args.outdir}")


if __name__ == "__main__":
    main()
