#!/usr/bin/env bash
################################################################################
# Mouse consensus pipeline -- Methods steps 2 (input prep) through 5
# (antiparallel consensus).
#
# Usage:
#   bash pipelines/run_mouse_consensus.sh <category> <input_log2fc.tsv>
#
#   <category> : synaptic | chromatin | mrna
#   <input>    : per-gene log2FC TSV (row 1 group label, row 2 sample ID,
#                row 3+ gene symbol + per-sample log2FC values)
################################################################################
set -euo pipefail

CATEGORY="${1:?Usage: run_mouse_consensus.sh <synaptic|chromatin|mrna> <input.tsv>}"
INPUT="${2:?missing input log2FC matrix}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPTS="$ROOT/scripts"
DATA="$ROOT/data"
OUT="$ROOT/results/mouse/$CATEGORY"
mkdir -p "$OUT"

# Pick the curated gene list for this category
case "$CATEGORY" in
  synaptic)  GENE_LIST="$DATA/list_genes_syngo.txt"    ;;
  chromatin) GENE_LIST="$DATA/list_genes_chromatin.txt" ;;
  mrna)      GENE_LIST="$DATA/list_genes_mRNA.txt"     ;;
  *) echo "Unknown category: $CATEGORY" >&2; exit 1 ;;
esac

# (2a) Reorder columns into [G1 | G2]
GROUPED="$OUT/grouped.tsv"
bash "$SCRIPTS/02_prepare_input/split_columns_by_group.sh" "$INPUT" "$GROUPED"

# (2b) Filter rows by category
FILTERED="$OUT/${CATEGORY}_filtered.tsv"
python3 "$SCRIPTS/02_prepare_input/filter_by_category.py" \
    "$GROUPED" "$GENE_LIST" "$FILTERED"

# (3) Split the filtered matrix into G1-only and G2-only matrices for the
#     per-group Wilcoxon test
INPUT_G1="$OUT/input_g1.tsv"
INPUT_G2="$OUT/input_g2.tsv"
awk -F'\t' -v g1="$INPUT_G1" -v g2="$INPUT_G2" '
NR<=2 {
  for (i=1; i<=NF; i++) {
    if (i==1) { printf "%s",$i > g1; printf "%s",$i > g2 }
    else if (NR==1) {
      if ($i ~ /Group1|G1|hGroup1|hG1/)      col[i]="g1"
      else if ($i ~ /Group2|G2|hGroup2|hG2/) col[i]="g2"
    }
    if (i>1) {
      if (col[i]=="g1") printf "\t%s",$i > g1
      if (col[i]=="g2") printf "\t%s",$i > g2
    }
  }
  printf "\n" > g1; printf "\n" > g2; next
}
{
  printf "%s",$1 > g1; printf "%s",$1 > g2
  for (i=2; i<=NF; i++) {
    if (col[i]=="g1") printf "\t%s",$i > g1
    if (col[i]=="g2") printf "\t%s",$i > g2
  }
  printf "\n" > g1; printf "\n" > g2
}' "$FILTERED"

# (3) Wilcoxon test per group
Rscript "$SCRIPTS/03_wilcoxon/wilcoxon_one_sided_BH.R" "$INPUT_G1" "$OUT/g1.pdf"
Rscript "$SCRIPTS/03_wilcoxon/wilcoxon_one_sided_BH.R" "$INPUT_G2" "$OUT/g2.pdf"

# Extract the per-gene call table
python3 "$SCRIPTS/03_wilcoxon/extract_calls.py" \
    -i "$OUT/g1.alpha0.050_eps0.00_wilcoxon_labeled_sorted.tsv" \
    -o "$OUT/g1_calls.tsv"
python3 "$SCRIPTS/03_wilcoxon/extract_calls.py" \
    -i "$OUT/g2.alpha0.050_eps0.00_wilcoxon_labeled_sorted.tsv" \
    -o "$OUT/g2_calls.tsv"

# (4-5) Antiparallel consensus
python3 "$SCRIPTS/04_antiparallel_consensus/overlap_antiparallel_genes.py" \
    --g1 "$OUT/g1_calls.tsv" --g2 "$OUT/g2_calls.tsv" \
    --category "$CATEGORY" \
    --out-genes  "$OUT/overlapped_genes.txt" \
    --out-updown "$OUT/overlapped_genes_updown.txt"

# Subset the matrix and run a clustering heatmap
python3 "$SCRIPTS/04_antiparallel_consensus/filter_input_by_overlapped.py" \
    -i "$FILTERED" -g "$OUT/overlapped_genes.txt" \
    -o "$OUT/filtered_overlapped_matrix.tsv"
Rscript "$SCRIPTS/04_antiparallel_consensus/heatmap_cluster.R" \
    "$OUT/filtered_overlapped_matrix.tsv" \
    "$OUT/overlapped_heatmap.pdf"

# Subset the per-group stats and merge them side-by-side
python3 "$SCRIPTS/04_antiparallel_consensus/extract_overlapped_stats.py" \
    -i "$OUT/g1_calls.tsv" -g "$OUT/overlapped_genes.txt" \
    -o "$OUT/g1_calls_overlapped.tsv"
python3 "$SCRIPTS/04_antiparallel_consensus/extract_overlapped_stats.py" \
    -i "$OUT/g2_calls.tsv" -g "$OUT/overlapped_genes.txt" \
    -o "$OUT/g2_calls_overlapped.tsv"
python3 "$SCRIPTS/04_antiparallel_consensus/merge_g1_g2.py" \
    -i1 "$OUT/g1_calls_overlapped.tsv" \
    -i2 "$OUT/g2_calls_overlapped.tsv" \
    -o  "$OUT/merged_g1_g2.tsv"

echo "[OK] mouse $CATEGORY pipeline done -> $OUT"
