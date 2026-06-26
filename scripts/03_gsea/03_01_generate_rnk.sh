#!/usr/bin/env bash
################################################################################
# 03.01 Rank file (RNK) generation from DEG table
# Author : Hyojin Kang
# Updated: 2026-05  (p-value=0 guard added)
#
# Reference: https://genomespot.blogspot.com/2015/01/how-to-generate-rank-file-from-gene.html
#
# Usage:
#   bash 03_01_generate_rnk.sh <DEG_table.tsv> <gene_id_col> <log2FC_col> <pvalue_col>
#
# Input  : DESeq2 output (tab-separated, with header row)
# Output : <DEG_table.tsv>.rnk  (gene<TAB>signed_-log10P, sorted)
################################################################################

set -euo pipefail

XLS="${1:?Usage: 03_01_generate_rnk.sh <table> <id_col> <fc_col> <p_col>}"
ID="${2:?gene_id column}"
FC="${3:?fold-change column}"
P="${4:?p-value column}"

# p-value=0 -> replaced with 1e-300 to avoid log(0) = -Inf
sed 1d "$XLS" | tr -d '"' \
  | awk -v I="$ID" -v F="$FC" -v P="$P" 'BEGIN{FS="\t"; OFS="\t"} $I!="" {print $I, $F, $P}' \
  | awk 'BEGIN{FS=OFS="\t"} $1!="NA" && $2!="NA" && $3!="NA" {if ($3+0==0) $3=1e-300; print}' \
  | awk 'BEGIN{FS=OFS="\t"} {s=1} $2<0{s=-1} {print $1, s*-1*log($3)/log(10)}' \
  | awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $2*$2}' | sort -k3gr \
  | awk 'BEGIN{FS=OFS="\t"} !arr[$1]++ {print $1, $2}' \
  | sort -k2gr > "${XLS}.rnk"

echo "[OK] Wrote ${XLS}.rnk ($(wc -l < "${XLS}.rnk") genes)"
