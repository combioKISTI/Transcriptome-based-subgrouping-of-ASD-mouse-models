#!/usr/bin/env bash
################################################################################
# 04. Alternative splicing analysis (rMATS)
# Author : Yukyung Jun
# Updated: 2026-05  (variable name fix: OUTDIR/TMPDIR consistency)
#
# Tool   : rMATS 4.0.2  (https://rnaseq-mats.sourceforge.net/)
#
# Required environment variables (see config/paths.sh):
#   WORKDIR, RMATS, GENOME_FASTA, GTF, STAR_INDEX
# Per-run:
#   MODEL    - mouse model (e.g., "Adnp")
#   SEX      - "M" or "F"
#   READ_LEN - read length of the FASTQ data (DATASET-DEPENDENT! adjust here)
################################################################################

set -euo pipefail

: "${WORKDIR:?}" ; : "${RMATS:?}" ; : "${GTF:?}" ; : "${STAR_INDEX:?}"
: "${MODEL:?}"   ; : "${SEX:?}"
READ_LEN="${READ_LEN:-101}"

# ---- Inputs (text files containing comma-separated FASTQ paths) -------------
MUT_FASTQ="$WORKDIR/fastq_list/${MODEL}_${SEX}_Mut.txt"
WT_FASTQ="$WORKDIR/fastq_list/${MODEL}_${SEX}_WT.txt"

# ---- Outputs ----------------------------------------------------------------
OUTDIR="$WORKDIR/splicing/${MODEL}/${MODEL}_${SEX}"
TMPDIR_R="$WORKDIR/splicing/tmp/${MODEL}_${SEX}"
mkdir -p "$OUTDIR" "$TMPDIR_R"

# ---- Run --------------------------------------------------------------------
# NOTE: --readLength must match your sequencing read length (varies per batch)
python "$RMATS" \
  --nthread 8 \
  --readLength "$READ_LEN" \
  --bi "$STAR_INDEX" \
  --s1 "$MUT_FASTQ" --s2 "$WT_FASTQ" \
  --gtf "$GTF" \
  --od "$OUTDIR" \
  --tmp "$TMPDIR_R"

echo "[OK] rMATS done -> $OUTDIR"
