#!/usr/bin/env bash
################################################################################
# Hypergeometric tests -- assess the significance of antiparallel overlaps.
#
# Default values are the human BA9 numbers reported in the manuscript. To
# evaluate other backgrounds (e.g. the mouse pipeline) edit N, K, n, x to
# match your dataset.
################################################################################
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT="$ROOT/scripts/05_hypergeometric/hypergeometric_test.R"
OUT="$ROOT/results/hypergeometric"
mkdir -p "$OUT"

run_one() {
  local label="$1" N="$2" K="$3" n="$4" x="$5"
  Rscript "$SCRIPT" --N "$N" --K "$K" --n "$n" --x "$x" --label "$label" \
      | tee "$OUT/${label}.txt"
  echo
}

# Human BA9 -- values reported in the manuscript text
run_one "human_synaptic"  1555 218 806 195
run_one "human_chromatin"  961  64  98  23
run_one "human_mrna"       947  66  74  27

# Add additional categories or species as needed, e.g.:
# run_one "mouse_synaptic"  <N> <K> <n> <x>

echo "[OK] hypergeometric tests written to $OUT"
