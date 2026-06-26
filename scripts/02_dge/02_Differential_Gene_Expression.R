################################################################################
# 02. Bulk RNA-seq processing & Differential Expressed Genes (DEGs)
# Author : Hyojin Kang
# Updated: 2026-05  (sva import; tx2gene fix; drug length fix; DIR -> wkdir;
#                    ComBat-seq output saved for downstream WGCNA)
#
# Pipeline:
#   (1) Load count data
#   (2) Batch correction (ComBat-seq)
#   (3) Gene filtering
#   (4) Normalization (DESeq2)
#   (5) Differential expressed genes
#
# DGE is performed independently for each mouse strain x sex.
# Detailed batch info: data/rnaseq_sample_metadata.txt
################################################################################

suppressPackageStartupMessages({
  library("Biobase")
  library("DESeq2")
  library("sva")        # ComBat_seq
  library("tximport")
  library("biomaRt")
  library("readr")
  library("ggplot2")
})

################################################################################
# (0) Configuration -- override via environment / command-line as needed
################################################################################
wkdir <- Sys.getenv("WORKDIR", unset = "/path/to/project/")
GENE  <- Sys.getenv("GENE",    unset = "ADNP")
TAG   <- Sys.getenv("TAG",     unset = "ADNP_Male")

setwd(file.path(wkdir, TAG))
out_dir <- file.path(wkdir, TAG)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# (1) Load count data
################################################################################
sampleList  <- list.files(file.path(wkdir, TAG, "salmon"))
sampleFiles <- file.path(wkdir, TAG, "salmon", sampleList, "quant.genes.sf")
n_samples   <- length(sampleFiles)

# --- Sample design (example: 30 samples = 15 HT + 15 WT, 6 groups x 5 reps) ----
# Adjust these vectors when the experimental layout differs.
genotype <- c(rep("HT", 15), rep("WT", 15))
drug     <- c(rep("fluoxetine", 5), rep("LiCl", 5), rep("vehicle", 5),
              rep("fluoxetine", 5), rep("LiCl", 5), rep("vehicle", 5))   # length 30 (fixed)
group    <- c(rep("HfM", 5), rep("HLM", 5), rep("HvM", 5),
              rep("WfM", 5), rep("WLM", 5), rep("WvM", 5))
batch    <- factor(c(rep("batch1", 5), rep("batch2", 5),
                     rep("batch1", 5), rep("batch2", 5),
                     rep("batch1", 5), rep("batch2", 5)))

stopifnot(length(genotype) == n_samples,
          length(drug)     == n_samples,
          length(group)    == n_samples,
          length(batch)    == n_samples)

samples <- data.frame(
  sampleName = sampleList,
  fileName   = sampleFiles,
  genotype   = genotype,
  drug       = drug,
  condition  = group,
  batch      = batch,
  stringsAsFactors = FALSE
)

# --- Ensembl annotation -------------------------------------------------------
ensembl_mm10 <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                        dataset = "mmusculus_gene_ensembl")
tx2gene <- getBM(
  attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
  mart       = ensembl_mm10
)
names(tx2gene) <- c("TXNAME", "GENEID")

annotation <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_accession",
                 "chromosome_name", "start_position", "end_position", "description"),
  mart = ensembl_mm10
)

# Note: salmon "quant.genes.sf" already aggregates to gene level; tx2gene is
# kept for reference but tximport reads gene-level quantifications directly.
txi <- tximport(sampleFiles, type = "salmon",
                tx2gene = tx2gene, ignoreTxVersion = TRUE)

################################################################################
# (2) Batch correction (ComBat-seq)
################################################################################
ddsTxi   <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
ddsCount <- counts(ddsTxi)
saveRDS(ddsCount, file.path(out_dir, "rdsDdsCount.RData"))

combatCount <- ComBat_seq(counts = ddsCount, batch = batch)
saveRDS(combatCount, file.path(out_dir, "combatCount.RData"))
write.table(combatCount,
            file = file.path(out_dir, sprintf("%s_exp_batch_correction.txt", TAG)),
            sep = "\t", quote = FALSE)

ddsCombat <- DESeqDataSetFromMatrix(countData = combatCount,
                                    colData   = samples,
                                    design    = ~ condition)

################################################################################
# (3) Gene filtering : keep genes with total count >= 5 * n_samples
################################################################################
keep      <- rowSums(counts(ddsCombat)) >= n_samples * 5
ddsCombat <- ddsCombat[keep, ]

################################################################################
# (4) Normalization (DESeq2)
################################################################################
dds <- DESeq(ddsCombat, parallel = TRUE)

################################################################################
# (5) DEGs : example contrast HvM vs WvM (vehicle, mutant vs wildtype)
################################################################################
res <- results(dds,
               contrast      = c("condition", "HvM", "WvM"),
               pAdjustMethod = "BH",
               alpha         = 0.05)
res <- res[!is.na(res$pvalue), ]

write.table(as.data.frame(res),
            file  = file.path(out_dir, sprintf("%s_deseq2_%s_all.txt", GENE, TAG)),
            sep   = "\t",
            quote = FALSE)

################################################################################
# (6) Reproducibility
################################################################################
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_02_DGE.txt"))
