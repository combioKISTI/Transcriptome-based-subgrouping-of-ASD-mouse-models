#!/usr/bin/env Rscript
################################################################################
# 06_cross_species -- Mouse to Human one-to-one ortholog mapping (BioMart)
#
# Manuscript Methods (verbatim):
#   For cross-species comparison, mouse and human log2FC matrices were
#   intersected on one-to-one orthologs using HGNC symbols mapped via
#   Ensembl BioMart release 102, retaining only genes detected in both
#   species.
#
# Output : TSV with columns
#   mouse_ensembl_gene_id | mouse_symbol | human_ensembl_gene_id | hgnc_symbol
# Only one-to-one orthologs with confidence = 1 are retained.
#
# Usage:
#   Rscript mouse_human_ortholog.R --input mouse_genes.txt --output orthologs.tsv
#
# If --input is omitted the full BioMart one-to-one ortholog table is written.
################################################################################
suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(input = NULL, output = "mouse_human_ortholog.tsv")
  i <- 1
  while (i <= length(args)) {
    key <- args[i]; val <- args[i + 1]
    if      (key == "--input")  out$input  <- val
    else if (key == "--output") out$output <- val
    i <- i + 2
  }
  out
}

args <- parse_args()

# Ensembl BioMart release 102 (manuscript specification)
mouse_mart <- useEnsembl(biomart = "genes",
                         dataset = "mmusculus_gene_ensembl",
                         version = 102)

ortho <- getBM(
  attributes = c("ensembl_gene_id",
                 "external_gene_name",
                 "hsapiens_homolog_ensembl_gene",
                 "hsapiens_homolog_associated_gene_name",
                 "hsapiens_homolog_orthology_type",
                 "hsapiens_homolog_orthology_confidence"),
  mart = mouse_mart
) %>%
  filter(hsapiens_homolog_orthology_type == "ortholog_one2one",
         hsapiens_homolog_orthology_confidence == 1,
         hsapiens_homolog_associated_gene_name != "") %>%
  rename(mouse_ensembl_gene_id  = ensembl_gene_id,
         mouse_symbol           = external_gene_name,
         human_ensembl_gene_id  = hsapiens_homolog_ensembl_gene,
         hgnc_symbol            = hsapiens_homolog_associated_gene_name) %>%
  select(mouse_ensembl_gene_id, mouse_symbol,
         human_ensembl_gene_id, hgnc_symbol) %>%
  distinct()

if (!is.null(args$input)) {
  query <- readLines(args$input)
  query <- query[query != ""]
  ortho <- ortho %>% filter(mouse_ensembl_gene_id %in% query | mouse_symbol %in% query)
}

write.table(ortho, file = args$output,
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           sub("\\.tsv$", ".sessionInfo.txt", args$output))
cat(sprintf("[OK] %d one-to-one orthologs -> %s\n", nrow(ortho), args$output))
