################################################################################
# 06. Human bulk RNA-seq analysis (Brodmann area BA9)
# Author : Hyojin Kang
# Updated: 2026-05  (samples-defined-before-use fix; biomaRt mapping activated;
#                    ComBat-seq output saved; clearer block ordering)
#
# Source data:
#   https://github.com/dhglab/Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/blob/master/data_provided/01_RNAseqProcessing/01_02_A_01_RawData.RData
#
# Pipeline:
#   (1) Load
#   (2) Gene filtering
#   (3) Sample filtering
#   (4) Build colData (samples)
#   (5) Batch correction (ComBat-seq)
#   (6) Normalization (DESeq2)
#   (7) Annotation (Ensembl -> HGNC via biomaRt)
#   (8) DEGs : Group1 / Group2 vs CTL (subtypes from co-expression clustering)
################################################################################

suppressPackageStartupMessages({
  library("DESeq2")
  library("sva")
  library("biomaRt")
})

wkdir <- Sys.getenv("WORKDIR", unset = ".")
out_dir <- file.path(wkdir, "results", "human")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# (1) Load
################################################################################
load(file.path(wkdir, "data", "human", "01_02_A_01_RawData.RData"))
# Provides: rsem_gene, rsem_gene_effLen, datMeta

################################################################################
# (2) Gene filtering : keep cpm > 0.1 in 30% of samples; remove effLen <= 15 bp
################################################################################
cpm <- apply(rsem_gene, 2, function(x) (x / sum(x)) * 1e6)
keep <- apply(cpm > 0.1, 1, sum)
cpm  <- cpm[keep > 0.3 * ncol(cpm), ]

eff_len <- rsem_gene_effLen[match(rownames(cpm), names(rsem_gene_effLen))]
short   <- which(eff_len <= 15)
if (length(short) > 0) cpm <- cpm[-short, ]

rsem_gene_effLen <- rsem_gene_effLen[match(rownames(cpm), names(rsem_gene_effLen))]
rsem_gene        <- rsem_gene[match(rownames(cpm), rownames(rsem_gene)), ]
rm(cpm, keep, short, eff_len)
cat(sprintf("[INFO] After gene filter: %d genes x %d samples\n",
            nrow(rsem_gene), ncol(rsem_gene)))

################################################################################
# (3) Sample filtering : BA9 + unstranded + ASD/CTL + dedup + outlier removal
################################################################################
datMeta.flt <- datMeta[datMeta$region == "BA9", ]
datMeta.flt <- datMeta.flt[datMeta.flt$SeqMethod == "unstranded", ]
datMeta.flt <- datMeta.flt[datMeta.flt$Diagnosis %in% c("ASD", "CTL"), ]

# Remove duplicated subject samples
datMeta.flt <- datMeta.flt[!grepl("_dup2", datMeta.flt$sample_id), ]
datMeta.flt <- datMeta.flt[!(datMeta.flt$sample_id %in%
                             c("UMB4334_BA9_2015-47", "UMB5302_BA9_2015-47")), ]

# Outliers excluded by sample clustering (matches original publication)
datMeta.flt <- datMeta.flt[!(datMeta.flt$sample_id %in%
                             c("AN01093_BA9_2012-150", "AN02987_BA9_2013-222")), ]

cat(sprintf("[INFO] After sample filter: %d samples\n", nrow(datMeta.flt)))
write.table(datMeta.flt,
            file = file.path(out_dir, "sample_table_flt.txt"),
            sep  = "\t", quote = FALSE)

rsem_gene.flt <- rsem_gene[, colnames(rsem_gene) %in% datMeta.flt$sample_id]
counts        <- round(rsem_gene.flt)
batch         <- factor(datMeta.flt$batch)

################################################################################
# (4) Build colData (samples) BEFORE any DESeqDataSet call
################################################################################
condition <- factor(datMeta.flt$Diagnosis)
samples   <- data.frame(sampleName = colnames(counts),
                        condition  = condition,
                        batch      = batch,
                        stringsAsFactors = FALSE)
write.table(samples,
            file = file.path(out_dir, "sample_table.txt"),
            sep  = "\t", quote = FALSE)

################################################################################
# (5) Batch correction (ComBat-seq)
################################################################################
ddsSet      <- DESeqDataSetFromMatrix(countData = counts,
                                      colData   = samples,
                                      design    = ~ condition)
combatCount <- ComBat_seq(counts = as.matrix(counts), batch = batch)
saveRDS(combatCount, file.path(out_dir, "combatCount_human.RData"))

################################################################################
# (6) Normalization
################################################################################
ddsCombat <- DESeqDataSetFromMatrix(countData = round(combatCount),
                                    colData   = samples,
                                    design    = ~ condition)
dds <- DESeq(ddsCombat)

normTable <- as.data.frame(counts(dds, normalized = TRUE))
write.table(normTable,
            file = file.path(out_dir, "rsem_gene_outlier_removed.txt"),
            sep  = "\t", quote = FALSE)

################################################################################
# (7) Annotation : Ensembl ID -> HGNC symbol
################################################################################
ens_ids        <- sub("\\..*", "", rownames(normTable))
mart           <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
genes_mapping  <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                        filters    = "ensembl_gene_id",
                        values     = unique(ens_ids),
                        mart       = mart)
saveRDS(genes_mapping, file.path(out_dir, "ensembl_to_hgnc.RData"))

normTable_anno <- normTable
normTable_anno$ensembl_gene_id <- ens_ids
normTable_anno <- merge(normTable_anno, genes_mapping,
                        by = "ensembl_gene_id", all.x = TRUE)
write.table(normTable_anno,
            file = file.path(out_dir, "rsem_gene_batch_corrected_anno.txt"),
            sep  = "\t", quote = FALSE, row.names = FALSE)

################################################################################
# (8) DEGs : Group1 / Group2 vs CTL
################################################################################
meta <- read.table(file.path(wkdir, "data", "human_rnaseq_sample_metadata.txt"),
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
condition_g <- factor(meta$group)
samples_g   <- data.frame(sampleName = colnames(combatCount),
                          condition  = condition_g,
                          stringsAsFactors = FALSE)

ddsCombatG <- DESeqDataSetFromMatrix(countData = round(combatCount),
                                     colData   = samples_g,
                                     design    = ~ condition)
ddsG <- DESeq(ddsCombatG)

write_group_degs <- function(dds, contrast_levels, label) {
  res <- results(dds,
                 contrast      = c("condition", contrast_levels[1], contrast_levels[2]),
                 pAdjustMethod = "BH",
                 alpha         = 0.05)
  res <- res[!is.na(res$pvalue), ]
  rownames(res) <- sub("\\..*", "", rownames(res))
  df <- as.data.frame(res)
  df$ensembl_gene_id <- rownames(df)
  df <- merge(df, genes_mapping, by = "ensembl_gene_id", all.x = TRUE)
  df <- df[!(is.na(df$hgnc_symbol) | df$hgnc_symbol == ""), ]
  out <- file.path(out_dir, sprintf("deseq2_%s.txt", label))
  write.table(df, file = out, sep = "\t",
              quote = FALSE, row.names = FALSE)
  cat(sprintf("[OK] %s -> %s (%d rows)\n", label, out, nrow(df)))
}

write_group_degs(ddsG, c("Group1", "CTL"), "Group1")
write_group_degs(ddsG, c("Group2", "CTL"), "Group2")

################################################################################
# (9) Reproducibility
################################################################################
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_06_human.txt"))
cat("[OK] Human RNA-seq outputs ->", out_dir, "\n")
