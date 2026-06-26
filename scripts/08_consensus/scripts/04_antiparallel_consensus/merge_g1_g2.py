#!/usr/bin/env python3
"""
04_antiparallel_consensus -- Merge G1 and G2 stats side-by-side with suffixes

For each species and category, joins the per-group Wilcoxon stat tables on
`Gene` so the downstream summary table has a single row per gene with all G1
and G2 values. Header layout becomes:
    Gene | MeanFC_G1 | P_less_G1 | Padj_less_G1 | P_greater_G1 | Padj_greater_G1
         | WilcoxonCall_G1
         | MeanFC_G2 | P_less_G2 | Padj_less_G2 | P_greater_G2 | Padj_greater_G2
         | WilcoxonCall_G2

Usage:
    python merge_g1_g2.py -i1 stats_g1.tsv -i2 stats_g2.tsv -o merged.tsv
"""
import argparse


# extract_calls.py output column indices
GENE_COL  = 0
USED_COLS = (0, 1, 2, 3, 4, 5, 6)   # Gene, MeanFC, P_less, Padj_less, P_greater, Padj_greater, WilcoxonCall


def read_tsv(path):
    header = None
    rows = {}
    with open(path, encoding="utf-8") as f:
        for i, line in enumerate(f):
            cols = line.rstrip("\n").split("\t")
            if i == 0:
                header = [cols[c] for c in USED_COLS]
                continue
            if max(USED_COLS) >= len(cols):
                continue
            rows[cols[GENE_COL]] = [cols[c] for c in USED_COLS]
    return header, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i1", "--input1", required=True)
    ap.add_argument("-i2", "--input2", required=True)
    ap.add_argument("-o",  "--output", required=True)
    args = ap.parse_args()

    h1, d1 = read_tsv(args.input1)
    h2, d2 = read_tsv(args.input2)

    # Keep one Gene column; suffix the rest with _G1 / _G2
    h1_suf = [h1[0]] + [h + "_G1" for h in h1[1:]]
    h2_suf = [h + "_G2" for h in h2[1:]]                     # drop G2 Gene column
    merged_header = h1_suf + h2_suf

    common = sorted(set(d1) & set(d2))
    with open(args.output, "w", encoding="utf-8") as out:
        out.write("\t".join(merged_header) + "\n")
        for gene in common:
            row = d1[gene] + d2[gene][1:]                    # drop G2 Gene column
            out.write("\t".join(row) + "\n")
    print(f"[OK] Merged {len(common)} genes -> {args.output}")


if __name__ == "__main__":
    main()
