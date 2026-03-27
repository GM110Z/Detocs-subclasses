#!/usr/bin/env python3

import argparse
import csv
import requests
import time
import subprocess
from collections import defaultdict


NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def run(cmd):
    print("[RUN]", cmd)
    subprocess.run(cmd, shell=True, check=True)


def fetch_ipg(pid):
    params = {
        "db": "protein",
        "id": pid,
        "rettype": "ipg",
        "retmode": "text"
    }
    r = requests.get(NCBI, params=params)
    r.raise_for_status()
    lines = [x for x in r.text.splitlines() if x.strip()]
    if len(lines) < 2:
        return []
    return list(csv.DictReader(lines, delimiter="\t"))


def is_refseq(row):
    return row.get("Source", "").lower() == "refseq"


def pick_rep(rows):
    def score(acc):
        if acc.startswith("WP_"): return 0
        if acc.startswith("NP_"): return 1
        if acc.startswith("YP_"): return 2
        if acc.startswith("XP_"): return 3
        return 9

    key = None
    for k in rows[0]:
        if "protein" in k.lower():
            key = k
            break

    if key is None:
        return rows[0]

    return sorted(rows, key=lambda r: score(r.get(key, "")))[0]


def merge_intervals(intervals):
    intervals = sorted(intervals)
    merged = []

    for s, e in intervals:
        if not merged:
            merged.append([s, e])
            continue

        last = merged[-1]
        if s <= last[1]:
            last[1] = max(last[1], e)
        else:
            merged.append([s, e])

    return merged


def detect_cols(row):
    keys = row.keys()

    def find(opts):
        for o in opts:
            for k in keys:
                if o in k.lower():
                    return k
        return None

    return (
        find(["nucleotide", "genomic", "accession"]),
        find(["start"]),
        find(["stop", "end"])
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input", required=True, help="protein ID list")
    p.add_argument("-k", "--kb", type=int, default=25)
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("--delay", type=float, default=0.3)
    p.add_argument("--largest-only", action="store_true")
    args = p.parse_args()

    flank = args.kb * 1000

    proteins = [x.strip() for x in open(args.input) if x.strip()]

    regions = defaultdict(list)

    print("[INFO] Fetching IPG + selecting representatives")

    for pid in proteins:
        try:
            rows = fetch_ipg(pid)
        except Exception as e:
            print("[WARN]", pid, e)
            continue

        ref = [r for r in rows if is_refseq(r)]
        if not ref:
            continue

        rep = pick_rep(ref)

        acc_col, start_col, stop_col = detect_cols(rep)
        if not acc_col or not start_col or not stop_col:
            continue

        try:
            acc = rep[acc_col]
            s = int(rep[start_col])
            e = int(rep[stop_col])
        except:
            continue

        s2 = max(1, min(s, e) - flank)
        e2 = max(s, e) + flank

        regions[acc].append((s2, e2))

        time.sleep(args.delay)

    print(f"[INFO] Regions collected: {sum(len(v) for v in regions.values())}")

    final = []

    for acc, ivals in regions.items():
        merged = merge_intervals(ivals)

        if args.largest_only:
            largest = max(merged, key=lambda x: x[1]-x[0])
            final.append((acc, largest[0], largest[1]))
        else:
            for s, e in merged:
                final.append((acc, s, e))

    print(f"[INFO] Final regions: {len(final)}")

    with open("regions.tsv", "w") as out:
        for r in final:
            out.write(f"{r[0]}\t{r[1]}\t{r[2]}\n")

    print("[INFO] Downloading regions")

    run(
        f"""parallel -a regions.tsv -C '\\t' -j{args.threads} --delay 0.4 """
        f""""efetch -db nuccore -id {{1}} -seq_start {{2}} -seq_stop {{3}} -format gbwithparts > {{1}}_{{2}}_{{3}}.gb" """
    )

    print("[DONE]")


if __name__ == "__main__":
    main()
