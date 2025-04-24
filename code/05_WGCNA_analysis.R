################################################################################
##### 05. Weighted Gene Co-expression Network Analysis (WGCNA)
##### April 2025, Yukyung Jun
################################################################################

################################################################################
### WGCNA-based Drug Response Transcriptomic Analysis
### This pipeline performs a co-expression network analysis using the WGCNA framework 
### to identify gene modules associated with drug response in mouse models.
### The analysis includes preprocessing, network construction, module-trait correlation,
### odds ratio enrichment with published gene sets, and functional interpretation 
### using GO term enrichment (enrichR).

################################################################################
##### Pipeline #####
### (1) Load data
### (2) WGCNA Network Construction
### (3) Correlate Modules with Sample Traits
### (4) Calculate Odds Ratio with External Gene Sets
### (5) Functional Enrichment with enrichR

################################################################################
##### load library
library("WGCNA")
library("dplyr")
library("tidyr")
library("stringr")
library('enrichR')
library('ggplot2')

##### file path
wkdir <- "/your/project/path"

################################################################################
##### (1) Load data
### Load normalized RNA-seq expression matrix and sample-level metadata.
### Expression values are expected to be gene-wise normalized (e.g., Combat-seq).
normdata <- read.table(sprintf('%s/data/bulk/exp_batch_correction.txt', wkdir), sep = '\t', header = TRUE, row.names = 1)
metadata <- read.table(sprintf('%s/data/bulk/metadata.txt', wkdir), sep = '\t', header = TRUE)

################################################################################
##### (2) WGCNA Network Construction
### Build gene co-expression modules using WGCNA's blockwiseModules function.
### The soft-thresholding power is selected based on scale-free topology fitting.

#### Transpose expression matrix (samples x genes)
norm.counts <- t(normdata)

#### Select soft-thresholding power to approximate scale-free topology
powers <- 1:30
sft <- pickSoftThreshold(as.data.frame(norm.counts), powerVector = powers, networkType = "signed", verbose = 5)

#### Choose first power where R^2 > 0.9
soft_power <- which(sft$fitIndices$SFT.R.sq > 0.9)

#### Construct the network and detect modules
soft_power <- which(sft$fitIndices$SFT.R.sq > 0.9)
bwnet <- blockwiseModules(as.data.frame(norm.counts),
                          power = soft_power,
                          networkType = "signed",
                          TOMType = "signed",
                          minModuleSize = 30,
                          mergeCutHeight = 0.25,
                          numericLabels = TRUE,
                          verbose = 3)

#### Extract module eigengenes (representative expression profiles)
module_eigengenes <- bwnet$MEs

################################################################################
##### (3) Correlate Modules with Sample Traits
### Associate module eigengenes with biological or treatment metadata using Pearson correlation.
### Traits include genotype and group assignment (e.g., mutant vs. WT).

#### Match and reformat trait data
metadata_target <- metadata[match(rownames(norm.counts), metadata$Sample),]

traits <- metadata_target %>%
  mutate(Sex = ifelse(Sex == "Male", 1, 0),
         Genotype = ifelse(Genotype == "Mutant", 1, 0),
         Group1.HT = ifelse(Group == "Group1" & Genotype == "Mutant", 1, 0),
         Group2.HT = ifelse(Group == "Group2" & Genotype == "Mutant", 1, 0)) %>%
  select(Sex, Genotype, Group1.HT, Group2.HT)

#### Compute correlations and p-values
moduleTraitCor <- cor(module_eigengenes, traits, use = 'p')
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(norm.counts))
moduleTraitFDR <- matrix(p.adjust(moduleTraitPvalue, method = "BH"), 
                         nrow = nrow(moduleTraitPvalue),
                         dimnames = dimnames(moduleTraitPvalue))

################################################################################
##### (4) Calculate Odds Ratio with External Gene Sets
### Evaluate WGCNA modules for enrichment against curated gene sets (e.g., ASD-related, WGCNA modules).
### Gene sets are provided in GMT format. Fisher's exact test is used to compute odds ratios.

#### Function to read GMT gene sets
read_gmt <- function(file) {
  lines <- readLines(file)
  lapply(lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    list(gene_set_name = parts[1], description = parts[2], genes = parts[-c(1,2)])
  })
}

#### Read and evaluate gene sets
gmt_files <- c("asd.v1.symbols.gmt", "Gandal 2022 WGCNA.gmt", "Gupta 2014 WGCNA.gmt")
or_all_mat_final <- list()

for (gmt_f in gmt_files) {
  gmt_mat <- read_gmt(sprintf('%s/ref/genesets/%s', wkdir, gmt_f))
  odds_matrix <- sapply(gmt_mat, function(gs) {
    sapply(colnames(module_eigengenes), function(module) {
      module_genes <- names(bwnet$colors)[bwnet$colors == as.numeric(gsub("ME", "", module))]
      a <- length(intersect(module_genes, gs$genes))
      b <- length(setdiff(module_genes, gs$genes))
      c <- length(setdiff(gs$genes, module_genes))
      d <- length(setdiff(rownames(exp_dat), union(module_genes, gs$genes)))
      fisher.test(matrix(c(a,b,c,d), 2))$estimate
    })
  })
  or_all_mat_final[[gmt_f]] <- odds_matrix
}

################################################################################
##### (5) Functional Enrichment with enrichR
### Perform GO term enrichment analysis for each module using enrichR.
### Databases include Biological Process, Molecular Function, and Cellular Component (2025 releases).
selected_dbs <- c("GO_Biological_Process_2025", "GO_Molecular_Function_2025", "GO_Cellular_Component_2025")

# Example: enrichment for Module 0 genes
module_id <- 0
module_genes <- names(bwnet$colors)[bwnet$colors == module_id]
enrich_results <- enrichr(module_genes, selected_dbs)

# Extract top GO terms
top_terms_out <- lapply(names(enrich_results), function(db) {
  enrich_results[[db]] %>%
    arrange(Adjusted.P.value) %>%
    slice_head(n = 5) %>%
    mutate(logFDR = -log10(Adjusted.P.value),
           Database = db)
}) %>% bind_rows()

# Visualization
ggplot(top_terms_out, aes(x = logFDR, y = reorder(Term, logFDR), fill = Database)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~Database, scales = "free_y") +
  theme_minimal() +
  labs(x = "-log10(FDR)", y = NULL)
