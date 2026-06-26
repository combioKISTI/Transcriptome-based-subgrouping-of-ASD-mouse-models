#!/usr/bin/env python3
"""
03_wilcoxon -- Extract the Wilcoxon call table from the labeled sorted TSV

`wilcoxon_one_sided_BH.R` writes a single TSV containing
  (a) the sorted log2FC matrix,
  (b) one blank separator line, and
  (c) a per-gene table with columns Gene / MeanFC / P_* / Padj_* / WilcoxonCall.

This script extracts section (c) so the next step (antiparallel consensus)
gets a clean tabular input.

Usage:
    python extract_calls.py -i <labeled_sorted.tsv> -o <calls.tsv>
"""
import argparse


def extract_from_gene(input_file: str, output_file: str, delimiter: str = "\t") -> None:
    write_flag = False
    with open(input_file, encoding="utf-8") as fin, \
         open(output_file, "w", encoding="utf-8") as fout:
        for line in fin:
            cols = line.rstrip("\n").split(delimiter)
            if not write_flag and cols and cols[0] == "Gene":
                write_flag = True
            if write_flag:
                fout.write(line)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Extract the Wilcoxon call table from the wilcoxon_one_sided_BH.R output"
    )
    ap.add_argument("-i", "--input",  required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()
    extract_from_gene(args.input, args.output)
