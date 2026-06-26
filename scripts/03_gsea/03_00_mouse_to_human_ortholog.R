################################################################################
# 03.00 Mouse -> Human ortholog mapping for DEG tables
# Author : Hyojin Kang
# Created: 2026-05  (new step -- previously the conversion was done manually
#                    and was not reproducible from the public repo)
#
# Input  : DESeq2 output table from script 02 (one row per Ensembl mouse gene)
# Output : Same table with HGNC symbol column. Used as input to
#          03_01_generate_rnk.sh and then GSEA against MSigDB human gene sets.
#
# Required environment variables:
#   WORKDIR, GENE, COMPTAG  (e.g. COMPTAG="ADNP_Male_HvM_vs_WvM")
################################################################################

suppressPackageStartupMessages({
  library("biomaRt")
  library("dplyr")
})

wkdir   <- Sys.getenv("WORKDIR", unset = "/path/to/project/")
GENE    <- Sys.getenv("GENE",    unset = "ADNP")
COMPTAG <- Sys.getenv("COMPTAG", unset = "ADNP_Male")

in_file  <- file.path(wkdir, COMPTAG,
                      sprintf("%s_deseq2_%s_all.txt", GENE, COMPTAG))
out_file <- file.path(wkdir, COMPTAG,
                      sprintf("%s_deseq2_%s_all_human.txt", GENE, COMPTAG))

deg <- read.table(in_file, header = TRUE, sep = "\t",
                  check.names = FALSE, stringsAsFactors = FALSE)
deg$ensembl_mouse_gene_id <- rownames(deg)

# --- Mouse -> Human via Ensembl Compara ortholog table ------------------------
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ortho <- getBM(
  attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene",
                 "hsapiens_homolog_associated_gene_name",
                 "hsapiens_homolog_orthology_type",
                 "hsapiens_homolog_orthology_confidence"),
  mart = mouse
)
# Keep only one2one orthologs with confidence == 1 (recommended for GSEA)
ortho <- ortho %>%
  filter(hsapiens_homolog_orthology_type    == "ortholog_one2one",
         hsapiens_homolog_orthology_confidence == 1,
         hsapiens_homolog_associated_gene_name != "") %>%
  select(ensembl_mouse_gene_id          = ensembl_gene_id,
         human_ensembl_gene_id           = hsapiens_homolog_ensembl_gene,
         hgnc_symbol                     = hsapiens_homolog_associated_gene_name) %>%
  distinct()

merged <- deg %>%
  inner_join(ortho, by = "ensembl_mouse_gene_id")

# Collapse duplicate human symbols by keeping the row with smallest p-value
if ("pvalue" %in% colnames(merged)) {
  merged <- merged %>%
    group_by(hgnc_symbol) %>%
    slice_min(pvalue, with_ties = FALSE, n = 1) %>%
    ungroup()
}

write.table(merged, file = out_file, sep = "\t",
            quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(wkdir, COMPTAG, "sessionInfo_03_00_ortholog.txt"))
cat(sprintf("[OK] %d mouse genes -> %d human one2one orthologs -> %s\n",
            nrow(deg), nrow(merged), out_file))
