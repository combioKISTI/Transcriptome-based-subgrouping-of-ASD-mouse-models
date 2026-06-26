#!/usr/bin/env bash
################################################################################
# 07. Cell Ranger multi (snRNA-seq alignment + counting)
# Author : Yukyung Jun
# Updated: 2026-05  (split out as standalone step; env-var consistency)
#
# Tool   : Cell Ranger 10.0.0  (cellranger multi)
#
# Produces the per-sample feature-barcode matrices consumed by
# 07_qc_integration.R ($WORKDIR/cellranger/outputs/<sample>/...).
#
# Required environment variables:
#   WORKDIR    - project root
#   CELLRANGER - path to the cellranger install dir (contains the binary)
#   TAG        - sample id for this run (DATASET-DEPENDENT! one sample per call)
#
# Per-sample multi config CSV is expected at:
#   $WORKDIR/cellranger/${TAG}_config.csv
# That CSV points to the reference transcriptome and FASTQ paths, e.g.:
#   reference : refdata-gex-mm10-2020-A
#   fastqs    : .../fastq_snRNA/<run>_RawData_Outs
#
# Run one sample per invocation (drive across samples with a loop or a
# scheduler array), e.g.:
#   for TAG in $(cat "$WORKDIR/cellranger/sample_list.txt"); do
#     TAG="$TAG" bash 07_cellranger_multi.sh
#   done
################################################################################

set -euo pipefail

: "${WORKDIR:?}" ; : "${CELLRANGER:?}" ; : "${TAG:?}"

CR_DIR="$WORKDIR/cellranger"
CONFIG_CSV="$CR_DIR/${TAG}_config.csv"

[[ -x "$CELLRANGER/cellranger" ]] || { echo "cellranger binary not found: $CELLRANGER/cellranger" >&2; exit 1; }
[[ -f "$CONFIG_CSV" ]]            || { echo "config CSV not found: $CONFIG_CSV" >&2; exit 1; }

# ---- Run --------------------------------------------------------------------
# cellranger multi writes its run directory under the current working dir,
# so cd into the cellranger workspace first.
cd "$CR_DIR"
"$CELLRANGER/cellranger" multi \
  --id="$TAG" \
  --csv="$CONFIG_CSV"

echo "[OK] Cell Ranger multi done -> $CR_DIR/$TAG"
