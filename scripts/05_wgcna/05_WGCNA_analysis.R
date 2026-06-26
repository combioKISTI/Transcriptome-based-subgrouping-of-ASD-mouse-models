################################################################################
# 05. Weighted Gene Co-expression Network Analysis (WGCNA)
# Author : Yukyung Jun
# Updated: 2026-05  (soft_power fallback, exp_dat fix, grey-module exclusion,
#                    explicit output saves)
#
# Pipeline:
#   (1) Load ComBat-seq corrected expression + metadata
#   (2) Network construction (signed)
#   (3) Module-trait correlation
#   (4) Odds ratio enrichment vs. external gene sets (Fisher)
#   (5) Functional enrichment via enrichR
################################################################################

suppressPackageStartupMessages({
  library("WGCNA")
  library("dplyr")
  library("tidyr")
  library("stringr")
  library("enrichR")
  library("ggplot2")
})

allowWGCNAThreads()

wkdir   <- Sys.getenv("WORKDIR", unset = "/path/to/project/")
out_dir <- file.path(wkdir, "results", "wgcna")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# (1) Load data
################################################################################
metadata <- read.table(file.path(wkdir, "data", "rnaseq_sample_metadata.txt"),
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
metadata_target <- metadata %>% filter(Drug == "Vehicle")

normdata <- read.table(file.path(wkdir, "data", "bulk", "exp_batch_correction.txt"),
                       sep = "\t", header = TRUE, row.names = 1,
                       check.names = FALSE)
normdata_target <- normdata[, colnames(normdata) %in% metadata_target$Sample]

stopifnot(ncol(normdata_target) > 0)

################################################################################
# (2) Network construction
################################################################################
norm.counts <- t(normdata_target)   # samples x genes

powers <- 1:30
sft <- pickSoftThreshold(as.data.frame(norm.counts),
                         powerVector = powers,
                         networkType = "signed",
                         verbose     = 5)

# First power crossing R^2 > 0.9 ; fall back to 12 if none reaches the threshold
candidate <- which(sft$fitIndices$SFT.R.sq > 0.9)
soft_power <- if (length(candidate) > 0) min(candidate) else 12L
message(sprintf("Soft thresholding power = %d", soft_power))

bwnet <- blockwiseModules(
  as.data.frame(norm.counts),
  power           = soft_power,
  networkType     = "signed",
  TOMType         = "signed",
  minModuleSize   = 30,
  mergeCutHeight  = 0.25,
  numericLabels   = TRUE,
  saveTOMs        = FALSE,
  verbose         = 3
)

module_eigengenes <- bwnet$MEs
saveRDS(bwnet, file.path(out_dir, "bwnet.RData"))
write.table(module_eigengenes,
            file.path(out_dir, "module_eigengenes.txt"),
            sep = "\t", quote = FALSE)
write.table(data.frame(gene = names(bwnet$colors), module = bwnet$colors),
            file.path(out_dir, "gene_module_assignment.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

################################################################################
# (3) Module-trait correlation
################################################################################
traits <- metadata_target %>%
  mutate(Sex       = ifelse(Sex == "Male", 1, 0),
         Genotype  = ifelse(Genotype == "Mutant", 1, 0),
         Group1.HT = ifelse(Group == "Group1" & Genotype == 1, 1, 0),
         Group2.HT = ifelse(Group == "Group2" & Genotype == 1, 1, 0)) %>%
  select(Sex, Genotype, Group1.HT, Group2.HT)

moduleTraitCor    <- cor(module_eigengenes, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(norm.counts))
moduleTraitFDR    <- matrix(p.adjust(moduleTraitPvalue, method = "BH"),
                            nrow     = nrow(moduleTraitPvalue),
                            dimnames = dimnames(moduleTraitPvalue))

write.table(moduleTraitCor,
            file.path(out_dir, "moduleTraitCor.txt"),
            sep = "\t", quote = FALSE)
write.table(moduleTraitFDR,
            file.path(out_dir, "moduleTraitFDR.txt"),
            sep = "\t", quote = FALSE)

################################################################################
# (4) Odds-ratio enrichment vs. published gene sets
################################################################################
read_gmt <- function(file) {
  lines <- readLines(file)
  lapply(lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    list(gene_set_name = parts[1],
         description   = parts[2],
         genes         = parts[-c(1, 2)])
  })
}

# Background = expressed genes in this dataset
background_genes <- rownames(normdata_target)

# Exclude grey (module 0 = unassigned) from enrichment
modules_to_test <- setdiff(unique(bwnet$colors), 0)

gmt_files <- c("asd.v1.symbols.gmt", "Gandal_2022_WGCNA.gmt", "Gupta_2014_WGCNA.gmt")
or_all_mat_final <- list()

for (gmt_f in gmt_files) {
  gmt_path <- file.path(wkdir, "ref", "genesets", gmt_f)
  if (!file.exists(gmt_path)) {
    warning("Skipping missing GMT: ", gmt_path); next
  }
  gmt_mat <- read_gmt(gmt_path)
  odds_matrix <- sapply(gmt_mat, function(gs) {
    sapply(modules_to_test, function(module_id) {
      module_genes <- names(bwnet$colors)[bwnet$colors == module_id]
      a <- length(intersect(module_genes, gs$genes))
      b <- length(setdiff(module_genes, gs$genes))
      c <- length(setdiff(gs$genes, module_genes))
      d <- length(setdiff(background_genes, union(module_genes, gs$genes)))
      fisher.test(matrix(c(a, b, c, d), nrow = 2))$estimate
    })
  })
  rownames(odds_matrix) <- paste0("ME", modules_to_test)
  colnames(odds_matrix) <- sapply(gmt_mat, `[[`, "gene_set_name")
  or_all_mat_final[[gmt_f]] <- odds_matrix
  write.table(odds_matrix,
              file.path(out_dir, sprintf("OR_%s.txt", gmt_f)),
              sep = "\t", quote = FALSE)
}

################################################################################
# (5) Functional enrichment with enrichR (per non-grey module)
################################################################################
selected_dbs <- c("GO_Biological_Process_2025",
                  "GO_Molecular_Function_2025",
                  "GO_Cellular_Component_2025")

enrich_summary <- list()
for (module_id in modules_to_test) {
  module_genes <- names(bwnet$colors)[bwnet$colors == module_id]
  if (length(module_genes) < 10) next
  res <- tryCatch(enrichr(module_genes, selected_dbs),
                  error = function(e) { warning(e); NULL })
  if (is.null(res)) next
  top <- lapply(names(res), function(db) {
    df <- res[[db]]
    if (nrow(df) == 0) return(NULL)
    df %>%
      arrange(Adjusted.P.value) %>%
      slice_head(n = 10) %>%
      mutate(Module   = sprintf("ME%d", module_id),
             Database = db,
             logFDR   = -log10(Adjusted.P.value))
  }) %>% bind_rows()
  enrich_summary[[as.character(module_id)]] <- top
}
enrich_all <- bind_rows(enrich_summary)
write.table(enrich_all,
            file.path(out_dir, "enrichR_top10_per_module.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

################################################################################
# (6) Reproducibility
################################################################################
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_05_WGCNA.txt"))
cat("[OK] WGCNA outputs ->", out_dir, "\n")
