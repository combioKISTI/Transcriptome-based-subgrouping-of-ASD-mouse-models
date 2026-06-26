# Tools & versions

| Component | Version | Notes |
|-----------|---------|-------|
| OS | Linux (Ubuntu 22.04 / RHEL 7) | Tested on KISTI HPC |
| R | 4.2.x | Required for recent DESeq2 / sva / WGCNA releases |
| Python | 3.8+ | rMATS dependency |
| Salmon | 1.1.0 | Mapping; library type ISR (dUTP, paired-end) |
| STAR | 2.7.x | Used internally by rMATS via `--bi` |
| rMATS | 4.0.2 | Alternative splicing |
| GSEA CLI | 4.3.3 | `gsea-cli.sh GSEAPreranked` |
| MSigDB | v2023.2.Hs | c5.go.bp / c5.go.cc / c5.go.mf collections |

## R / Bioconductor packages

Captured per-script via `sessionInfo()`. Key dependencies:

- `Biobase`, `DESeq2`, `sva` (ComBat-seq)
- `tximport`, `biomaRt`, `readr`, `ggplot2`, `dplyr`, `tidyr`, `stringr`
- `WGCNA`, `enrichR`
- `ComplexHeatmap`, `circlize` (used in step 08 visualisations)

## External datasets

| Source | Used in | Access |
|--------|---------|--------|
| MSigDB v2023.2.Hs (Human) | 03 GSEA, 08 consensus | https://www.gsea-msigdb.org/gsea/msigdb |
| SynGO release 2023-12-01 | 08 consensus | https://www.syngoportal.org |
| Mouse GRCm38 reference | 01 / 04 | Ensembl |
| Ensembl BioMart release 102 | 03 / 08 ortholog mapping | https://www.ensembl.org/biomart |
| dhglab Broad transcriptomic dysregulation (BA9) | 07 human RNA-seq | https://github.com/dhglab/Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD |
