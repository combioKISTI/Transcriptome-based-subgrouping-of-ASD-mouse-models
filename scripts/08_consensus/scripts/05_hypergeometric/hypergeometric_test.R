#!/usr/bin/env Rscript
################################################################################
# 05_hypergeometric -- Hypergeometric test of antiparallel overlap significance
#
# Manuscript Methods (verbatim):
#   The statistical significance of each antiparallel overlap (e.g., G1-down
#   intersect G2-up for synaptic) was assessed by a one-sided hypergeometric
#   test (phyper() in R, lower.tail = FALSE) using the expressed genes
#   retained in the WGCNA input matrix as the common background universe
#   across all three functional categories.
#
# Parameters:
#   N : background universe size (number of expressed genes used as WGCNA input)
#   K : number of genes called in one direction in one group (e.g. G2-down)
#   n : number of genes called in the opposite direction in the other group
#       (e.g. G1-up)
#   x : observed antiparallel overlap size
#
# Usage:
#   Rscript hypergeometric_test.R --N 1555 --K 218 --n 806 --x 195 \
#                                 --label "human_synaptic"
################################################################################
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(label = "category")
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    val <- args[i + 1]
    if (key %in% c("--N", "--K", "--n", "--x")) {
      out[[sub("^--", "", key)]] <- as.numeric(val)
    } else if (key == "--label") {
      out$label <- val
    }
    i <- i + 2
  }
  required <- c("N", "K", "n", "x")
  missing  <- setdiff(required, names(out))
  if (length(missing) > 0) {
    stop(sprintf("Missing argument(s): %s\nUsage: --N <int> --K <int> --n <int> --x <int> [--label <str>]",
                 paste(missing, collapse = ", ")))
  }
  out
}

args <- parse_args()
N <- args$N; K <- args$K; n <- args$n; x <- args$x

p_value    <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
expected   <- n * K / N
enrichment <- x / expected

cat(sprintf("=== Hypergeometric test [%s] ===\n", args$label))
cat(sprintf("Background N        : %d\n",   N))
cat(sprintf("Set K               : %d\n",   K))
cat(sprintf("Set n               : %d\n",   n))
cat(sprintf("Observed overlap x  : %d\n",   x))
cat(sprintf("Expected overlap    : %.2f\n", expected))
cat(sprintf("Fold enrichment     : %.2fx\n", enrichment))
cat(sprintf("Hypergeometric p    : %.3e\n", p_value))
