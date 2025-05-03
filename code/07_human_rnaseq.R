################################################################################
##### 07. Human bulk RNA-seq data analysis
##### April 2025, Hyojin Kang
################################################################################

################################################################################
### Human bulk-RNA-seq data from Brodmann area were processed analogously: 
### after batch correction, log2 fold-changes were computed against age-matched neurotypical controls, 
### pairwise Pearson correlations were calculated across all individuals, 
### and the resulting matrix was subjected to the same clustering workflow

################################################################################
##### Pipeline #####
### (1) Load data
### (2) Gene Filtering
### (3) Sample Filtering
### (4) Batch Correction
### (5) Normalization
### (6) Annotation (add gene symbol)
### (7) Differential Expressed Genes (DEGs)

################################################################################
##### Load library
library("DESeq2")
library("sva")
library("biomaRt")

################################################################################
##### (1) Load data
### RNA-seq expression count matrix and sample metadata was downloaded from github link below
### Download: https://github.com/dhglab/Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/blob/master/data_provided/01_RNAseqProcessing/01_02_A_01_RawData.RData
load("01_02_A_01_RawData.RData") # counts and metadata

################################################################################
##### (2) Gene Filtering
### Keep these samples only for the gene filtering step
### Keep genes with cpm > 0.1 in 30% of samples
cpm <- apply(rsem_gene,2, function(x) (x/sum(x))*1000000) 
keep = apply(cpm>0.1,1,sum) 
idx = which(keep > 0.3*dim(cpm)[2]) ## cpm > 0.1 in 30% of samples
cpm = cpm[idx,]
### 26,475 genes genes remain - about 44% of genes pass filter

### Remove any remaining genes with an effective gene length <= 15 bp
length=rsem_gene_effLen[match(rownames(cpm),names(rsem_gene_effLen))] 
idx = which(length <= 15)
short_remove = rownames(cpm)[idx]
idx_remove = match(short_remove,rownames(cpm))
cpm = cpm[-idx_remove,]
rsem_gene_effLen=rsem_gene_effLen[match(rownames(cpm),names(rsem_gene_effLen))]
rsem_gene = rsem_gene[match(rownames(cpm),rownames(rsem_gene)),]
rm(cpm,length,keep,short_remove,idx,idx_remove)
dim(rsem_gene)
## 26364 genes, 854 samples

################################################################################
##### (3) Sample Filtering
### Keep these samples only for the sample filtering step

### Keep brain region matched with "BA9"
datMeta.flt <- datMeta[datMeta$region == "BA9",]
dim(datMeta.flt)
### 150 samples remain

### Keep SeqMethod matched with "unstranded"
datMeta.flt <- datMeta.flt[datMeta.flt$SeqMethod == "unstranded",]
dim(datMeta.flt)
### 104 samples remain

### Keep Diagnois matched with "ASD or CTL"
datMeta.flt <- datMeta.flt[(datMeta.flt$Diagnosis == "ASD" | datMeta.flt$Diagnosis  ==  "CTL"),]
dim(datMeta.flt)
### 93 samples remain

### Remove duplicated subject samples
datMeta.flt <- datMeta.flt[!grepl("_dup2", datMeta.flt[, "sample_id"]), ]
### Remove additional duplicated samples
datMeta.flt <- datMeta.flt[!(datMeta.flt[, "sample_id"] %in% c("UMB4334_BA9_2015-47", "UMB5302_BA9_2015-47")), ]

### Remove outlier samples
### Two outlier samples were identified based on sample clustering and excluded from downstream analyses. 
### Notably, these same samples had also been excluded in the original study from which the data were obtained.
datMeta.flt <- datMeta.flt[!(datMeta.flt[, "sample_id"] %in% c("AN01093_BA9_2012-150", "AN02987_BA9_2013-222")), ]

dim(datMeta.flt)
### 85 samples remain
write.table(datMeta.flt, file="sample_table_flt.txt", sep="\t", quote=F)

################################################################################
##### (4) Batch Correction
### Make submatrix matched filtered samples in step (3)
rsem_gene.flt <- rsem_gene[, colnames(rsem_gene) %in% datMeta.flt$sample_id]
dim(rsem_gene.flt)
## 26364 genes, 84 samples
#write.table(rsem_gene.flt, file="rsem_gene_flt_R.txt", sep="\t", quote=F)

### Batch Correction using ComBat_seq
counts <- rsem_gene.flt
dim(counts)
ddsSet <- DESeqDataSetFromMatrix(countData = round(counts), colData = samples, design = ~ condition)
#[1] 26364    84
batch <- factor(datMeta.flt$batch)
combatCount <- ComBat_seq(counts = as.matrix(counts), batch = batch)

################################################################################
##### (5) Normalization using DESeq2
condition <- factor(datMeta.flt$Diagnosis)
samples <- data.frame(sampleName = colnames(counts), condition = condition, batch = batch)
write.table(samples, file="sample_table.txt", sep="\t", quote=F)

ddsCombat <- DESeqDataSetFromMatrix(countData=round(combatCount), colData = samples, design = ~condition)
dds <- DESeq(ddsCombat)
#save(dds,file="dds.RData")
#dds <- load("dds.RData")

################################################################################
##### (6) Annotation (add gene symbol)
# remove ensembl ID version info 
normTable <- as.data.frame(counts(dds, normalized=TRUE))
write.table(normTable, file="rsem_gene_outlier_removed.txt", sep="\t", quote = FALSE)

#rownames(normTable) <- sub("\\..*", "", rownames(normTable))
# Add gene symbol using Ensembl BioMart
#normTable$ensembl_gene_id <- rownames(normTable)
#mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
#genes_mapping <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),filters = 'ensembl_gene_id',values = rownames(normTable),mart = mart)
#normTable_annotated <- merge(normTable, genes_mapping, by = "ensembl_gene_id", all.x = TRUE)
#write.table(normTable_annotated, file="rsem_gene_batch_corrected_anno.txt", sep="\t", quote = FALSE)

################################################################################
##### (7) Differential Expressed Genes (DEGs)
#### ASD samples were stratified into two subtypes, Group 1 and Group 2, 
#### using hierarchical clustering of gene co-expression networks 
meta <- read.table("human_rnaseq_sample_metadata.txt", header = TRUE, sep = "\t")
condition <- factor(meta$group)
meta$sample_id

samples <- data.frame(sampleName = colnames(combatCount), condition = condition)
ddsCombat <- DESeqDataSetFromMatrix(countData=round(combatCount), colData = samples, design = ~condition)
dds <- DESeq(ddsCombat)

##### Group1 from co-expression analysis
resG1 <- results(dds, contrast=c("condition", "Group1", "CTL"), pAdjustMethod = "BH", alpha=0.05)
resG1 <- resG1[!is.na(resG1$pvalue),]
head(rownames(resG1))
rownames(resG1) <- sub("\\..*", "", rownames(resG1))
resG1_df <- as.data.frame(resG1)
# add ensembl ID 
resG1_df$ensembl_gene_id <- rownames(resG1_df)
# add gene symbol 
resG1_annotated <- merge(resG1_df, genes_mapping, by = "ensembl_gene_id", all.x = TRUE)
resG1_annotated <- resG1_annotated[!(is.na(resG1_annotated$hgnc_symbol) | resG1_annotated$hgnc_symbol == ""), ]
write.table(as.data.frame(resG1_annotated), file="deseq2_Group1.txt", sep="\t", quote = FALSE)

##### Group2 from co-expression analysis
resG2 <- results(dds, contrast=c("condition", "Group2", "CTL"), pAdjustMethod = "BH", alpha=0.05)
resG2 <- resG2[!is.na(resG2$pvalue),]
head(rownames(resG2))
rownames(resG2) <- sub("\\..*", "", rownames(resG2))
resG2_df <- as.data.frame(resG2)
# add ensembl ID 
resG2_df$ensembl_gene_id <- rownames(resG2_df)
# add gene symbol 
resG2_annotated <- merge(resG2_df, genes_mapping, by = "ensembl_gene_id", all.x = TRUE)
resG2_annotated <- resG2_annotated[!(is.na(resG2_annotated$hgnc_symbol) | resG2_annotated$hgnc_symbol == ""), ]
write.table(as.data.frame(resG2_annotated), file="deseq2_Group2.txt", sep="\t", quote = FALSE)
##### 


