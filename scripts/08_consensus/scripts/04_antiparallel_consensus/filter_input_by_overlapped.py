#!/usr/bin/env python3
"""
04_antiparallel_consensus -- Filter the input log2FC matrix to overlapped genes

Produces a smaller log2FC matrix containing only the antiparallel consensus
genes, suitable as input to `heatmap_cluster.R`.

Usage:
    python filter_input_by_overlapped.py -i input_matrix.tsv \
                                         -g overlapped_genes.txt \
                                         -o filtered_matrix.tsv

Input matrix layout: row 1 group label, row 2 sample ID, row 3+ gene + log2FC.
"""
import argparse


def load_genes(path: str) -> set:
    s = set()
    with open(path, encoding="utf-8") as f:
        for line in f:
            tok = line.strip().split()
            if tok:
                s.add(tok[0])
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input",  required=True, help="Input log2FC matrix")
    ap.add_argument("-g", "--genes",  required=True, help="Overlapped gene list (one per line)")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    genes = load_genes(args.genes)

    with open(args.input, encoding="utf-8") as fin, \
         open(args.output, "w", encoding="utf-8") as fout:
        for i, line in enumerate(fin):
            if i < 2:           # preserve the two header rows
                fout.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            if cols and cols[0] in genes:
                fout.write(line)
    print(f"[OK] Filtered to {len(genes)} overlapped genes -> {args.output}")


if __name__ == "__main__":
    main()
