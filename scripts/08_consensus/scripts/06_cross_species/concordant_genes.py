#!/usr/bin/env python3
"""
06_cross_species -- Cross-species concordant antiparallel consensus genes

Manuscript Methods (verbatim):
  Antiparallel consensus genes in each species were derived as defined above.
  Cross-species concordant genes were then defined as antiparallel consensus
  genes in mouse whose orthologs were also antiparallel consensus genes in
  human, with matched directionality between corresponding groups
  (mouse G1 <-> human hGroup 1; mouse G2 <-> human hGroup 2).
  Effect-size similarity between species was summarized as the absolute
  difference between species-level median log2FCs.

Inputs:
  --mouse-overlap : mouse antiparallel consensus genes (overlap_antiparallel_genes.py output)
  --human-overlap : human antiparallel consensus genes
  --ortholog      : mouse_human_ortholog.R output (mouse_symbol <-> hgnc_symbol)
  --mouse-stats   : (optional) mouse merged G1+G2 stats TSV (for median log2FC comparison)
  --human-stats   : (optional) human merged stats TSV
  --out           : output TSV path

Output columns:
  mouse_symbol, hgnc_symbol [, mouse_median_FC, human_median_FC, abs_delta]

Usage:
  python concordant_genes.py \
      --mouse-overlap results/mouse/synaptic/overlapped.txt \
      --human-overlap results/human/synaptic/overlapped.txt \
      --ortholog      results/cross/mouse_human_ortholog.tsv \
      --mouse-stats   results/mouse/synaptic/merged.tsv \
      --human-stats   results/human/synaptic/merged.tsv \
      --out           results/cross/synaptic_concordant.tsv
"""
import argparse
import csv


def load_genes(path):
    s = set()
    with open(path, encoding="utf-8") as f:
        for line in f:
            tok = line.strip().split()
            if tok:
                s.add(tok[0])
    return s


def load_ortholog(path):
    mouse_to_human = {}
    with open(path, encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ms = row.get("mouse_symbol")
            hs = row.get("hgnc_symbol")
            if ms and hs:
                mouse_to_human[ms] = hs
    return mouse_to_human


def load_median_fc(path, group_suffix):
    """Pull the MeanFC_<G1|G2> column from a merge_g1_g2.py output table."""
    fc = {}
    if not path:
        return fc
    target = f"MeanFC_{group_suffix}"
    with open(path, encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if target not in reader.fieldnames:
            return fc
        for row in reader:
            try:
                fc[row["Gene"]] = float(row[target])
            except (TypeError, ValueError):
                continue
    return fc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mouse-overlap", required=True)
    ap.add_argument("--human-overlap", required=True)
    ap.add_argument("--ortholog",      required=True)
    ap.add_argument("--mouse-stats",   default=None)
    ap.add_argument("--human-stats",   default=None)
    ap.add_argument("--mouse-group",   default="G1",
                    help="Mouse group suffix to compare (default: G1)")
    ap.add_argument("--human-group",   default="G1",
                    help="Human group suffix to compare (default: G1)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    mouse_set = load_genes(args.mouse_overlap)
    human_set = load_genes(args.human_overlap)
    ortholog  = load_ortholog(args.ortholog)

    mouse_fc = load_median_fc(args.mouse_stats, args.mouse_group)
    human_fc = load_median_fc(args.human_stats, args.human_group)

    rows = []
    for m_gene in sorted(mouse_set):
        h_gene = ortholog.get(m_gene)
        if h_gene and h_gene in human_set:
            m_val = mouse_fc.get(m_gene)
            h_val = human_fc.get(h_gene)
            delta = (abs(m_val - h_val)
                     if m_val is not None and h_val is not None else None)
            rows.append((m_gene, h_gene, m_val, h_val, delta))

    with open(args.out, "w", encoding="utf-8") as out:
        out.write("mouse_symbol\thgnc_symbol\tmouse_median_FC\thuman_median_FC\tabs_delta\n")
        for m, h, mv, hv, d in rows:
            out.write(f"{m}\t{h}\t"
                      f"{'' if mv is None else f'{mv:.4f}'}\t"
                      f"{'' if hv is None else f'{hv:.4f}'}\t"
                      f"{'' if d  is None else f'{d:.4f}'}\n")

    print(f"[OK] {len(rows)} cross-species concordant genes -> {args.out}")


if __name__ == "__main__":
    main()
