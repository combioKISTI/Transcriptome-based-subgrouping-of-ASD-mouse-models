################################################################################
# 02b. Export ComBat-seq corrected count matrix for downstream WGCNA (script 05)
# Author : Hyojin Kang
# Updated: 2026-05  (new helper script -- consolidates outputs across all
#                    strain x sex DGE runs into a single matrix)
#
# Required environment variables:
#   WORKDIR   - project root
#   TAGS      - whitespace-separated list of strain_sex tags (e.g.
#               "ADNP_Male ADNP_Female CHD8_Male ...")
#
# Output:
#   $WORKDIR/data/bulk/exp_batch_correction.txt
#   $WORKDIR/data/bulk/exp_batch_correction.RData
################################################################################

suppressPackageStartupMessages({
  library("Biobase")
})

wkdir <- Sys.getenv("WORKDIR", unset = "/path/to/project/")
tags  <- strsplit(Sys.getenv("TAGS", unset = "ADNP_Male"), "\\s+")[[1]]

out_dir <- file.path(wkdir, "data", "bulk")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

merged <- NULL
for (tag in tags) {
  rds <- file.path(wkdir, tag, "combatCount.RData")
  if (!file.exists(rds)) {
    warning("Missing combatCount.RData for ", tag, " -- run 02_DGE first.")
    next
  }
  m <- readRDS(rds)
  colnames(m) <- paste(tag, colnames(m), sep = "__")
  if (is.null(merged)) {
    merged <- m
  } else {
    common <- intersect(rownames(merged), rownames(m))
    merged <- cbind(merged[common, , drop = FALSE], m[common, , drop = FALSE])
  }
}

stopifnot(!is.null(merged))

write.table(merged,
            file  = file.path(out_dir, "exp_batch_correction.txt"),
            sep   = "\t", quote = FALSE)
saveRDS(merged, file.path(out_dir, "exp_batch_correction.RData"))

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_02b_export.txt"))
cat(sprintf("[OK] Wrote %d genes x %d samples to %s\n",
            nrow(merged), ncol(merged), out_dir))
