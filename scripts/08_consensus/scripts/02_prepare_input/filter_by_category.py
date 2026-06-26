#!/usr/bin/env python3
"""
02_prepare_input -- Filter rows by functional category gene list

Curation sources used in the manuscript:
  - synaptic       : SynGO release 2023-12-01
  - chromatin      : MSigDB v2023.2.Hs GOBP "chromatin organization" + descendants
  - mRNA processing: MSigDB v2023.2.Hs GOBP "mRNA processing" + descendants

Usage:
    python filter_by_category.py <input.tsv> <gene_list.txt> <output.tsv>

The input TSV must have two header rows (group label, sample ID) followed by
data rows. The gene list is one gene symbol per line.
"""
import sys
import pandas as pd
from io import StringIO


def main():
    if len(sys.argv) < 4:
        sys.exit("Usage: filter_by_category.py <input.tsv> <gene_list.txt> <output.tsv>")

    input_file, keyword_file, output_file = sys.argv[1:4]

    with open(input_file, encoding="utf-8") as f:
        lines = f.readlines()
    if len(lines) < 3:
        sys.exit("Invalid input format: expected 2 header rows plus data rows.")

    header1 = lines[0].rstrip("\n")
    header2 = lines[1].rstrip("\n")
    df = pd.read_csv(StringIO("".join(lines[2:])), sep="\t", header=None, dtype=str)

    with open(keyword_file, encoding="utf-8") as f:
        keywords = {line.strip() for line in f if line.strip()}

    filtered = df[df.iloc[:, 0].isin(keywords)]

    with open(output_file, "w", encoding="utf-8") as out:
        out.write(header1 + "\n")
        out.write(header2 + "\n")
        filtered.to_csv(out, sep="\t", index=False, header=False)

    sys.stderr.write(
        f"[OK] {len(filtered)} / {len(df)} genes kept after filtering -> {output_file}\n"
    )


if __name__ == "__main__":
    main()
