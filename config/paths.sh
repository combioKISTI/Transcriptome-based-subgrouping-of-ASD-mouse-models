#!/usr/bin/env bash
################################################################################
# Common environment variables for all pipeline scripts
# Source this file at the top of each script:
#   source "$(dirname "$0")/../../config/paths.sh"
################################################################################

# Project root (override via env if needed)
export WORKDIR="${WORKDIR:-/path/to/project}"

# Reference (mouse: GRCm38)
export REF="${REF:-$WORKDIR/ref/GRCm38}"
export GENOME_FASTA="${GENOME_FASTA:-$REF/genome.fa}"
export GTF="${GTF:-$REF/gtf/genes.gtf}"
export SALMON_INDEX="${SALMON_INDEX:-$REF/salmon_index}"
export STAR_INDEX="${STAR_INDEX:-$REF/STAR}"

# Tools (set absolute paths)
export SALMON="${SALMON:-/path/salmon/salmon-1.1.0/bin/salmon}"
export RMATS="${RMATS:-$WORKDIR/tools/rMATS/rmats.py}"

# GSEA
export GSEA_HOME="${GSEA_HOME:-/path/GSEA}"
export GSEA_VERSION="${GSEA_VERSION:-v2023.2.Hs}"
export GSEA_RUN="${GSEA_RUN:-$GSEA_HOME/GSEA_Linux_4.3.3/gsea-cli.sh GSEAPreranked}"

# Per-run variables (must be exported before invoking each script)
#   TAG       : sample/run identifier (e.g., "ADNP_Male")
#   GENE      : gene model identifier  (e.g., "ADNP")
#   COMPTAG   : comparison label (e.g., "ADNP_Male_HvM_vs_WvM")
#   SAMPLE    : run label for outputs

# Threading
export NTHREAD="${NTHREAD:-8}"

# Print summary when sourced interactively
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "[paths.sh] WORKDIR=$WORKDIR"
  echo "[paths.sh] REF=$REF"
  echo "[paths.sh] SALMON=$SALMON"
  echo "[paths.sh] GSEA_HOME=$GSEA_HOME"
fi
