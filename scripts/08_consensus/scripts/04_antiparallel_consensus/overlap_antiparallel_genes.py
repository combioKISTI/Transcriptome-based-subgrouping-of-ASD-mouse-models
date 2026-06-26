#!/usr/bin/env python3
"""
04_antiparallel_consensus -- Define antiparallel consensus genes

Manuscript Methods (verbatim):
  Within each species, antiparallel consensus genes for a functional category
  were defined as those called significant in opposite directions in the two
  groups:
    synaptic                : (WilcoxonDown in G1) intersect (WilcoxonUp   in G2)
    chromatin / mRNA proc.  : (WilcoxonUp   in G1) intersect (WilcoxonDown in G2)

Inputs are the per-group call TSVs produced by `extract_calls.py`. Outputs:
  (1) overlapped_genes.txt        -- one gene per line (antiparallel consensus)
  (2) overlapped_genes_updown.txt -- Gene<TAB>G1_call<TAB>G2_call (for plots)

Usage:
    python overlap_antiparallel_genes.py \
        --g1 mouse_synaptic_g1_calls.tsv \
        --g2 mouse_synaptic_g2_calls.tsv \
        --category synaptic \
        --out-genes results/overlapped_synaptic.txt \
        --out-updown results/overlapped_synaptic_updown.txt
"""
import argparse
import csv

CATEGORY_RULES = {
    "synaptic":  {"G1": "WilcoxonDown", "G2": "WilcoxonUp"},
    "chromatin": {"G1": "WilcoxonUp",   "G2": "WilcoxonDown"},
    "mrna":      {"G1": "WilcoxonUp",   "G2": "WilcoxonDown"},
}


def load_calls(path: str) -> dict:
    """Return {gene: WilcoxonCall} dict from an extract_calls.py output TSV."""
    out = {}
    with open(path, encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "WilcoxonCall" not in reader.fieldnames or "Gene" not in reader.fieldnames:
            raise SystemExit(f"{path}: 'Gene' and 'WilcoxonCall' columns are required.")
        for row in reader:
            out[row["Gene"]] = row["WilcoxonCall"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--g1", required=True, help="Group1 calls TSV (extract_calls.py output)")
    ap.add_argument("--g2", required=True, help="Group2 calls TSV")
    ap.add_argument("--category", required=True, choices=CATEGORY_RULES.keys())
    ap.add_argument("--out-genes",  required=True)
    ap.add_argument("--out-updown", required=True)
    args = ap.parse_args()

    rule = CATEGORY_RULES[args.category]
    g1 = load_calls(args.g1)
    g2 = load_calls(args.g2)

    common = sorted(set(g1) & set(g2))
    overlapped = [g for g in common
                  if g1[g] == rule["G1"] and g2[g] == rule["G2"]]

    with open(args.out_genes, "w", encoding="utf-8") as f:
        f.write("\n".join(overlapped) + "\n")

    with open(args.out_updown, "w", encoding="utf-8") as f:
        f.write("Gene\tG1_call\tG2_call\n")
        for g in overlapped:
            f.write(f"{g}\t{g1[g]}\t{g2[g]}\n")

    print(f"[OK] {args.category}: {len(overlapped)} antiparallel consensus genes "
          f"(rule: G1={rule['G1']}, G2={rule['G2']}) -> {args.out_genes}")


if __name__ == "__main__":
    main()
