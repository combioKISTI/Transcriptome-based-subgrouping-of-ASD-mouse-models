# Data

This folder contains sample-level metadata only. Raw FASTQ / count files must
be downloaded separately as described below.

## Files in this folder

| File | Description | Rows |
|------|-------------|------|
| `rnaseq_sample_metadata.txt` | Mouse bulk RNA-seq sample sheet (1,008 samples; Mouse_Model, Sex, Genotype, Drug, Batch, Group) | 1,008 |
| `human_rnaseq_sample_metadata.txt` | Human BA9 cohort sample sheet (Group1 / Group2 / CTL labels from co-expression clustering) | 83 |

## External datasets

### 1. Mouse bulk RNA-seq (FASTQ)
- **Source:** GEO accession (TBD when manuscript is publicly indexed)
- **Place under:** `$WORKDIR/fastq/`
- Naming convention: `<sample>_1.fastq.gz` / `<sample>_2.fastq.gz`
  where `<sample>` corresponds to the `Sample` column of `rnaseq_sample_metadata.txt`.

### 2. Human bulk RNA-seq (counts + metadata)
- **Source:** [dhglab/Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD](https://github.com/dhglab/Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD)
- **File:** `data_provided/01_RNAseqProcessing/01_02_A_01_RawData.RData`
- **Place under:** `data/human/01_02_A_01_RawData.RData`
- Contains `rsem_gene`, `rsem_gene_effLen`, `datMeta` objects.

### 3. MSigDB GMT files (human)
- **Source:** https://www.gsea-msigdb.org/gsea/msigdb
- **Version:** v2023.2.Hs
- **Place under:** `$GSEA_HOME/msigdb_v2023.2.Hs_GMTs/`
- Categories used: `c5.go.bp`, `c5.go.cc`, `c5.go.mf`
- Chip file: `Human_HGNC_ID_MSigDB.v2023.2.Hs.chip` under
  `$GSEA_HOME/msigdb_v2023.2.Hs_chip_files_to_download_locally/`

### 4. WGCNA reference gene sets (used in script 05 OR analysis)
- `asd.v1.symbols.gmt`
- `Gandal_2022_WGCNA.gmt`
- `Gupta_2014_WGCNA.gmt`
- **Place under:** `$WORKDIR/ref/genesets/`

### 5. Step 08 functional category gene lists (already shipped)
- `scripts/08_consensus/data/list_genes_syngo.txt`     -- SynGO release 2023-12-01
- `scripts/08_consensus/data/list_genes_chromatin.txt` -- MSigDB GOBP "chromatin organization" descendants
- `scripts/08_consensus/data/list_genes_mRNA.txt`      -- MSigDB GOBP "mRNA processing" descendants
