############################################################
##### 06 Cell type deconvolution of bulk RNA-seq 
##### April 2025, Hyojin Kang

############################################################
### Cell type deconvolution using Bisque
### Bulk RNA sequencing counts were deconvolved using the reference cell profiles from previously published single cell data (Yao et al., 2023). 
### The 8 major cell classes (GABAPRO_NKX2-1, GABAPRO_ASCL1, CTX_GLU_LAYER1, GABAPRO_DLX1-2, CTX_GLU_LAYER5, GABA_CCK, GABA_VIP, GABA_PVALB) 
### identified in a large meta-analysis of 4 million brain single cells. 
### The reference-based decomposition was performed using the ReferenceBasedDecomposition function from the BisqueRNA package (v1.0.5) (Jew et al., 2020).

############################################################
##### Reference singlecell data
### Reference: Yao, Z. et al. A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain. 
##  Nature 624, 317-332 (2023). https://doi.org/10.1038/s41586-023-06812-z

### The single-cell gene expression matrix was generated by downloading four donor samples (79,798 single cells) 
### matched for age, brain region, and genotype from the ALLEN BRAIN cell type database
## Download: https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#expression_matrices/WMB-10Xv3/)
## donor ID: 371230, 372312, 372314, 546812

############################################################
##### BisqueRNA package (v1.0.5)
### https://cran.r-project.org/web/packages/BisqueRNA/vignettes/bisque.html

############################################################
##### load library
library(Biobase)
library(BisqueRNA)
library(plyr)
library(knitr)
library(ggplot2)
library(gplots)

##### file path
wkdir <- "/path/WMB-10X/"
sc_path <- paste(wkdir, "donors/", sep="")
mixture_path <- paste(wkdir, "mixture/", sep="")

############################################################
##### Load bulk-RNAseq
### Bulk RNA-seq data can be converted from a matrix (columns are samples, rows are genes) to an ExpressionSet as follows:
mixture_file_raw <- paste(mixture_path, "htseq_norm_raw.txt", sep="")
bulk.data.raw <- read.table(mixture_file_raw, header = TRUE, row.names = 1, sep = "\t")

###  Convert the expression data to a matrix (if it's not already)
bulk.matrix <- as.matrix(bulk.data.raw)

##########################
sample_data <- data.frame(condition = colnames(bulk.matrix), row.names = colnames(bulk.matrix))
gene_data <- data.frame(gene_type = rep("Unknown", nrow(bulk.matrix)), row.names = rownames(bulk.matrix))
bulk.eset <- Biobase::ExpressionSet(assayData = bulk.matrix, phenoData = AnnotatedDataFrame(sample_data), featureData = AnnotatedDataFrame(gene_data))

############################################################
##### Load sc-RNAseq
#Single-cell data requires additional information in the ExpressionSet, specificially cell type labels and individual labels. 
#Individual labels indicate which individual each cell originated from
#To add this information, Biobase requires it to be stored in a data frame format. 
#Assuming we have character vectors of cell type labels (cell.type.labels) and individual labels (individual.labels), 

############################################################
##### cols: celltypes, rows: ensID
DONOR.1="546812"; 
sc_file_1 <- paste(sc_path, "merged/step4-raw-",DONOR.1, ".tsv", sep="")
sc.data.1 <- read.table(sc_file_1, header = TRUE, row.names = 1, sep = "\t")
sc.counts.matrix.1 <- as.matrix(sc.data.1)
celltype_file_1 <- paste(sc_path, "merged/step2-celltype_",DONOR.1, ".txt", sep="")
celltype.data.1 <- read.table(celltype_file_1, header = TRUE, sep = "\t")
cell.type.labels.1 <- celltype.data.1[,1]
length(cell.type.labels.1)
id_file_1 <- paste(sc_path, "merged/step2-sampleid_",DONOR.1, ".txt", sep="")
id.data.1 <- read.table(id_file_1, header = TRUE, row.names = 1, sep = "\t")
sample.ids.1 <- rownames(id.data.1)
length(sample.ids.1)
individual.labels.1 <- rep(DONOR.1, length(sample.ids.1)) 
colnames(sc.counts.matrix.1) <- sample.ids.1

#############################
DONOR.2="372314"; 
sc_file_2 <- paste(sc_path, "merged/step4-raw-",DONOR.2, ".tsv", sep="")
sc.data.2 <- read.table(sc_file_2, header = TRUE, row.names = 1, sep = "\t")
sc.counts.matrix.2 <- as.matrix(sc.data.2)
celltype_file_2 <- paste(sc_path, "merged/step2-celltype_",DONOR.2, ".txt", sep="")
celltype.data.2 <- read.table(celltype_file_2, header = TRUE, sep = "\t")
cell.type.labels.2 <- celltype.data.2[,1]
length(cell.type.labels.2)
id_file_2 <- paste(sc_path, "merged/step2-sampleid_",DONOR.2, ".txt", sep="")
id.data.2 <- read.table(id_file_2, header = TRUE, row.names = 1, sep = "\t")
sample.ids.2 <- rownames(id.data.2)
length(sample.ids.2)
individual.labels.2 <- rep(DONOR.2, length(sample.ids.2)) 
colnames(sc.counts.matrix.2) <- sample.ids.2

#############################
DONOR.3="371230"; 
sc_file_3 <- paste(sc_path, "merged/step4-raw-",DONOR.3, ".tsv", sep="")
sc.data.3 <- read.table(sc_file_3, header = TRUE, row.names = 1, sep = "\t")
sc.counts.matrix.3 <- as.matrix(sc.data.3)
celltype_file_3 <- paste(sc_path, "merged/step2-celltype_",DONOR.3, ".txt", sep="")
celltype.data.3 <- read.table(celltype_file_3, header = TRUE, sep = "\t")
cell.type.labels.3 <- celltype.data.3[,1]
length(cell.type.labels.3)
id_file_3 <- paste(sc_path, "merged/step2-sampleid_",DONOR.3, ".txt", sep="")
id.data.3 <- read.table(id_file_3, header = TRUE, row.names = 1, sep = "\t")
sample.ids.3 <- rownames(id.data.3)
length(sample.ids.3)
individual.labels.3 <- rep(DONOR.3, length(sample.ids.3)) 
colnames(sc.counts.matrix.3) <- sample.ids.3

#############################
DONOR.4="372312"; 
sc_file_4 <- paste(sc_path, "merged/step4-raw-",DONOR.4, ".tsv", sep="")
sc.data.4 <- read.table(sc_file_4, header = TRUE, row.names = 1, sep = "\t")
sc.counts.matrix.4 <- as.matrix(sc.data.4)
celltype_file_4 <- paste(sc_path, "merged/step2-celltype_",DONOR.4, ".txt", sep="")
celltype.data.4 <- read.table(celltype_file_4, header = TRUE, sep = "\t")
cell.type.labels.4 <- celltype.data.4[,1]
length(cell.type.labels.4)
id_file_4 <- paste(sc_path, "merged/step2-sampleid_",DONOR.4, ".txt", sep="")
id.data.4 <- read.table(id_file_4, header = TRUE, row.names = 1, sep = "\t")
sample.ids.4 <- rownames(id.data.4)
length(sample.ids.4)
individual.labels.4 <- rep(DONOR.4, length(sample.ids.4)) 
colnames(sc.counts.matrix.4) <- sample.ids.4

#############################
sc.counts.matrix <- cbind(sc.counts.matrix.1, sc.counts.matrix.2,sc.counts.matrix.3,sc.counts.matrix.4)
sample.ids <- c(sample.ids.1, sample.ids.2, sample.ids.3, sample.ids.4)
individual.labels <- c(individual.labels.1, individual.labels.2, individual.labels.3, individual.labels.4)
cell.type.labels <- c(cell.type.labels.1, cell.type.labels.2, cell.type.labels.3, cell.type.labels.4)
############################################################
#convert scRNA-seq data (with counts also in matrix format) as follows:
sc.pheno <- data.frame(check.names=F, check.rows=F, stringsAsFactors=F,
                       row.names=sample.ids, SubjectName=individual.labels, cellType=cell.type.labels)
sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),row.names=c("SubjectName","cellType"))
sc.pdata <- new("AnnotatedDataFrame", data=sc.pheno, varMetadata=sc.meta)
sc.eset <- Biobase::ExpressionSet(assayData=sc.counts.matrix, phenoData=sc.pdata)

############################################################
##### Estimate celltype decomposition
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=FALSE)
ref.based.estimates <- res$bulk.props
knitr::kable(ref.based.estimates, digits=2)

##### Write results
write.csv(ref.based.estimates, file = paste(sc_path, "merged/out-bulk-all-props.csv", sep=""), row.names = TRUE)


