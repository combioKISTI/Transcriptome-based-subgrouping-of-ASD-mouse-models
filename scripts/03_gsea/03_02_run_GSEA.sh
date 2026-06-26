#!/usr/bin/env bash
################################################################################
# 03.02 Gene Set Enrichment Analysis (GSEA Preranked)
# Author : Hyojin Kang
# Updated: 2026-05  (env var documentation; loop over GO categories)
#
# Tool   : GSEA Linux 4.3.3 (gsea-cli.sh GSEAPreranked)
# DB     : MSigDB v2023.2.Hs (https://www.gsea-msigdb.org/gsea/msigdb)
# Input  : *.rnk produced by 03_01_generate_rnk.sh on the human-orthologized
#          DEG table from 03_00_mouse_to_human_ortholog.R
#
# Required environment variables (see config/paths.sh):
#   WORKDIR, GSEA_HOME, GSEA_VERSION, GSEA_RUN
# Per-run:
#   GENE     - mouse gene model (e.g., "ADNP")
#   COMPTAG  - comparison label (e.g., "ADNP_Male_HvM_vs_WvM")
#   SAMPLE   - run label
################################################################################

set -euo pipefail

: "${WORKDIR:?}" ; : "${GSEA_HOME:?}" ; : "${GSEA_VERSION:?}" ; : "${GSEA_RUN:?}"
: "${GENE:?}"    ; : "${COMPTAG:?}"   ; : "${SAMPLE:?}"

HUMAN_CHIP_DIR="msigdb_${GSEA_VERSION}_chip_files_to_download_locally"
HUMAN_GMT_DIR="msigdb_${GSEA_VERSION}_GMTs"
HUMAN_CHIP="Human_HGNC_ID_MSigDB.${GSEA_VERSION}.chip"

RNK="$WORKDIR/${COMPTAG}/${GENE}_deseq2_${COMPTAG}_all_human.txt.rnk"
[[ -s "$RNK" ]] || { echo "Missing rnk: $RNK -- run 03_00 then 03_01 first." >&2; exit 1; }

OUTPUT="$WORKDIR/gsea/${SAMPLE}"
mkdir -p "$OUTPUT"

run_gsea_category() {
  local CAT="$1"   # e.g. "c5.go.bp"
  local LABEL="${SAMPLE}_${CAT//./_}"
  local GMT="$GSEA_HOME/$HUMAN_GMT_DIR/${CAT}.${GSEA_VERSION}.symbols.gmt"

  $GSEA_RUN \
    -gmx "$GMT" -rnk "$RNK" -chip "$HUMAN_CHIP" -collapse false \
    -mode Max_probe -norm meandiv -nperm 1000 \
    -scoring_scheme classic -rpt_label "$LABEL" \
    -include_only_symbols true -make_sets true -plot_top_x 30 \
    -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false \
    -out "$OUTPUT"
}

for CAT in c5.go.bp c5.go.cc c5.go.mf; do
  echo "[GSEA] running $CAT ..."
  run_gsea_category "$CAT"
done

echo "[OK] GSEA completed -> $OUTPUT"
