#!/usr/bin/env bash
################################################################################
# 01. Bulk RNA-seq raw FASTQ mapping using Salmon
# Author : Hyojin Kang
# Updated: 2026-05  (re-organized layout, robust loop, logging)
#
# Required environment variables (see config/paths.sh):
#   WORKDIR, SALMON, SALMON_INDEX, GTF
# Per-run:
#   TAG  - sample identifier (paired FASTQ named ${TAG}_1.fastq.gz / _2.fastq.gz)
#
# Usage:
#   source ../../config/paths.sh
#   for TAG in $(cat sample_list.txt); do
#     export TAG && bash 01_FASTQ_mapping_by_salmon.sh
#   done
################################################################################

set -euo pipefail

# ---- Configuration ----------------------------------------------------------
# salmon version 1.1.0  (https://combine-lab.github.io/salmon/)
: "${WORKDIR:?Set WORKDIR (e.g. via config/paths.sh)}"
: "${SALMON:?Set SALMON (path to salmon binary)}"
: "${SALMON_INDEX:?Set SALMON_INDEX}"
: "${GTF:?Set GTF}"
: "${TAG:?Set TAG (sample identifier)}"

NTHREAD="${NTHREAD:-8}"

# ---- I/O --------------------------------------------------------------------
FASTQ1="$WORKDIR/fastq/${TAG}_1.fastq.gz"
FASTQ2="$WORKDIR/fastq/${TAG}_2.fastq.gz"
SALMON_OUT="$WORKDIR/salmon/${TAG}"
LOG_DIR="$WORKDIR/logs/salmon"

mkdir -p "$SALMON_OUT" "$LOG_DIR"

# ---- Sanity checks ----------------------------------------------------------
[[ -s "$FASTQ1" ]] || { echo "Missing FASTQ1: $FASTQ1" >&2; exit 1; }
[[ -s "$FASTQ2" ]] || { echo "Missing FASTQ2: $FASTQ2" >&2; exit 1; }
[[ -d "$SALMON_INDEX" ]] || { echo "Missing index: $SALMON_INDEX" >&2; exit 1; }

# ---- Run --------------------------------------------------------------------
# --libType ISR : strand-specific dUTP paired-end
#   I = Inward (paired-end facing each other)
#   S = Stranded
#   R = Second read in same direction as the transcript
"$SALMON" quant \
  -i "$SALMON_INDEX" \
  --libType ISR \
  -p "$NTHREAD" \
  --seqBias --validateMappings --mimicBT2 --gcBias \
  -1 "$FASTQ1" -2 "$FASTQ2" \
  -o "$SALMON_OUT" \
  -g "$GTF" \
  2> "$LOG_DIR/${TAG}.log"

echo "[OK] $TAG -> $SALMON_OUT"
