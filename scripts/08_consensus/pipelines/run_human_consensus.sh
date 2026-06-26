#!/usr/bin/env bash
################################################################################
# Human consensus pipeline -- Methods steps 1 (ID mapping) through 5
# (antiparallel consensus).
#
# Usage:
#   bash pipelines/run_human_consensus.sh <category> <input_log2fc.tsv>
#
# If the first column of the input matrix looks like an Ensembl gene ID
# (starts with ENSG), the Ensembl-to-symbol mapping step is run automatically.
# If the first column is already an HGNC symbol, the pipeline behaves the
# same as the mouse pipeline.
################################################################################
set -euo pipefail

CATEGORY="${1:?Usage: run_human_consensus.sh <synaptic|chromatin|mrna> <input.tsv>}"
INPUT="${2:?missing input log2FC matrix}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPTS="$ROOT/scripts"
DATA="$ROOT/data"
OUT="$ROOT/results/human/$CATEGORY"
mkdir -p "$OUT"

# Detect ID style from the first data row; if it is an Ensembl ID, run mapping
FIRST_GENE=$(awk -F'\t' 'NR==3 {print $1; exit}' "$INPUT")
WORK_INPUT="$INPUT"
if [[ "$FIRST_GENE" == ENSG* ]]; then
  echo "[INFO] First column looks like an Ensembl ID -- running ensembl_to_symbol"
  IDS_TXT="$OUT/ensembl_ids.txt"
  awk -F'\t' 'NR>2 {print $1}' "$INPUT" > "$IDS_TXT"
  MAP_TXT="$OUT/ensembl_to_symbol.tsv"
  python3 "$SCRIPTS/01_id_mapping/ensembl_to_symbol.py" "$IDS_TXT" "$MAP_TXT"

  MAPPED="$OUT/input_mapped.tsv"
  awk -F'\t' -v m="$MAP_TXT" '
    BEGIN {
      while ((getline line < m) > 0) {
        if (line ~ /^input_id\t/) continue
        n = split(line, t, "\t"); map_e[t[1]] = t[3]
      }
    }
    NR<=2 { print; next }
    {
      sym = map_e[$1]; if (sym == "" || sym == "NA") next
      $1 = sym
      OFS = "\t"; print
    }' "$INPUT" > "$MAPPED"
  WORK_INPUT="$MAPPED"
fi

# From step 2 onward the pipeline is identical to the mouse one, so reuse it
bash "$ROOT/pipelines/run_mouse_consensus.sh" "$CATEGORY" "$WORK_INPUT" \
  || true

# Move the results into the human output directory
SRC_OUT="$ROOT/results/mouse/$CATEGORY"
if [ -d "$SRC_OUT" ] && [ "$(ls -A "$SRC_OUT" 2>/dev/null)" ]; then
  mv "$SRC_OUT"/* "$OUT/" 2>/dev/null || true
  rmdir "$SRC_OUT" 2>/dev/null || true
fi
echo "[OK] human $CATEGORY pipeline done -> $OUT"
