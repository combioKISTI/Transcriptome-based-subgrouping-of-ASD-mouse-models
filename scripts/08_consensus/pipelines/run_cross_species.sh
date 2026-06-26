#!/usr/bin/env bash
################################################################################
# Cross-species concordance -- Methods step 7
#
# Usage:
#   bash pipelines/run_cross_species.sh <category>
#
# Prerequisite: `run_mouse_consensus.sh` and `run_human_consensus.sh` must
# have already been executed for the same category.
################################################################################
set -euo pipefail

CATEGORY="${1:?Usage: run_cross_species.sh <synaptic|chromatin|mrna>}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPTS="$ROOT/scripts/06_cross_species"
OUT="$ROOT/results/cross_species/$CATEGORY"
mkdir -p "$OUT"

MOUSE_DIR="$ROOT/results/mouse/$CATEGORY"
HUMAN_DIR="$ROOT/results/human/$CATEGORY"

[[ -s "$MOUSE_DIR/overlapped_genes.txt" ]] || { echo "missing $MOUSE_DIR/overlapped_genes.txt -- run the mouse pipeline first" >&2; exit 1; }
[[ -s "$HUMAN_DIR/overlapped_genes.txt" ]] || { echo "missing $HUMAN_DIR/overlapped_genes.txt -- run the human pipeline first" >&2; exit 1; }

# (a) Build (or reuse) the mouse-to-human one-to-one ortholog table from
#     Ensembl BioMart release 102
ORTHO="$OUT/mouse_human_ortholog.tsv"
if [[ ! -s "$ORTHO" ]]; then
  echo "[INFO] Querying Ensembl BioMart release 102 for one-to-one orthologs..."
  Rscript "$SCRIPTS/mouse_human_ortholog.R" \
      --input "$MOUSE_DIR/overlapped_genes.txt" \
      --output "$ORTHO"
fi

# (b) Concordant genes with matched directionality
python3 "$SCRIPTS/concordant_genes.py" \
    --mouse-overlap "$MOUSE_DIR/overlapped_genes.txt" \
    --human-overlap "$HUMAN_DIR/overlapped_genes.txt" \
    --ortholog "$ORTHO" \
    --mouse-stats "$MOUSE_DIR/merged_g1_g2.tsv" \
    --human-stats "$HUMAN_DIR/merged_g1_g2.tsv" \
    --out "$OUT/${CATEGORY}_concordant.tsv"

echo "[OK] cross-species $CATEGORY done -> $OUT"
