#!/usr/bin/env bash
################################################################################
# 02_prepare_input -- Reorder sample columns into [Group1 | Group2] order
#
# Input : per-gene log2FC TSV (row 1 = group label, row 2 = sample ID,
#         row 3+ = gene name in column 1 followed by per-sample log2FC values).
#         CRLF line endings are normalised on the fly.
# Output: same layout but sample columns are sorted into Group1 first, then
#         Group2. The gene column is preserved as column 1.
#
# Manuscript-defined input dimensions:
#   - mouse  : G1 = 15 (mutant line × sex combinations), G2 = 19
#   - human  : hGroup 1 = 20 ASD individuals, hGroup 2 = 20
#
# This script identifies group membership from the labels in row 1
# (matches Group1 / G1 / hGroup1 / hG1 vs Group2 / G2 / hGroup2 / hG2).
# It generalises the dataset-specific `split_group.sh` shims used in the
# original ad-hoc scripts.
#
# Usage:
#   bash split_columns_by_group.sh <input.tsv> <output.tsv>
################################################################################
set -euo pipefail

IN="${1:?Usage: split_columns_by_group.sh <input.tsv> <output.tsv>}"
OUT="${2:?missing output path}"

# Strip CRLF (Windows line endings) before processing so awk regex matches
# Group2 / Group1 cleanly in the trailing column too.
TMP=$(mktemp)
trap 'rm -f "$TMP"' EXIT
tr -d '\r' < "$IN" > "$TMP"

awk -F'\t' -v out="$OUT" '
NR==1 {
    # Identify G1 / G2 column indices from the row 1 group labels
    n = NF
    g1_idx[1] = 1; g1_n = 1   # column 1 (gene name) is always preserved
    g2_n = 0
    for (i = 2; i <= n; i++) {
        if ($i ~ /Group1|G1|hGroup1|hG1/)      { g1_n++; g1_idx[g1_n] = i }
        else if ($i ~ /Group2|G2|hGroup2|hG2/) { g2_n++; g2_idx[g2_n] = i }
    }
}
{
    line = ""
    for (k = 1; k <= g1_n; k++) line = line (k > 1 ? "\t" : "") $(g1_idx[k])
    for (k = 1; k <= g2_n; k++) line = line "\t" $(g2_idx[k])
    print line > out
}
END {
    print "[OK] Reordered " (g1_n - 1) " G1 + " g2_n " G2 columns -> " out > "/dev/stderr"
}
' "$TMP"
