#!/usr/bin/env python3
"""
04_antiparallel_consensus -- Extract per-gene Wilcoxon stats limited to overlapped genes

Filters an `extract_calls.py` output TSV down to the antiparallel consensus
gene set. The resulting per-group tables are then combined by `merge_g1_g2.py`.

Usage:
    python extract_overlapped_stats.py -i calls_g1.tsv -g overlapped_genes.txt \
                                       -o calls_g1_overlapped.tsv
"""
import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input",  required=True)
    ap.add_argument("-g", "--genes",  required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    with open(args.genes, encoding="utf-8") as f:
        gene_set = {line.strip().split()[0] for line in f if line.strip()}

    with open(args.input, encoding="utf-8") as fin, \
         open(args.output, "w", encoding="utf-8") as fout:
        for i, line in enumerate(fin):
            if i == 0:
                fout.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            if cols and cols[0] in gene_set:
                fout.write(line)
    print(f"[OK] Filtered stats to {len(gene_set)} genes -> {args.output}")


if __name__ == "__main__":
    main()
