# 07_snrnaseq — Single-nucleus RNA-seq analysis

Cell-type-resolved snRNA-seq pipeline for the ASD mouse-model cortex: from
Cell Ranger alignment through QC, Harmony integration, reference-based cell-type
annotation, cell-type composition analysis, and per-cell-type hdWGCNA.

## Analysis steps

| Step | Description | Script |
|------|-------------|--------|
| **(0)** Alignment + counting | `cellranger multi` per sample → feature-barcode matrices | `07_cellranger_multi.sh` |
| **(1–4)** QC, merge, integration, clustering | Per-nucleus QC (nCount / nFeature / percent.mt) → per-batch merge → normalize / 3,000 HVG / PCA → Harmony (batch) → SNN clustering (res 0.3) → UMAP → drop small clusters (retain 0–18) | `07_qc_integration.R` |
| **(5–10)** Annotation + composition | FindAllMarkers → Allen Brain Cell Atlas label transfer (FindTransferAnchors → MapQuery) → prediction-score thresholding → manual cluster→cell-type mapping (19 clusters → 10 types) → cell-type proportion analysis (Wilcoxon; B6/G1/G2 and Hv/Hf/HL) | `07_annotation.R` |
| **(11)** Cell-type-specific co-expression networks | Mouse→human ortholog → metacells (k=25) → soft power → signed-hybrid hdWGCNA (minModuleSize 50) → module eigengenes / hub genes / module-trait correlation → GMT overlap (Fisher) | `07_hdWGCNA.R` |
| — | Shared logging + plotting utilities (sourced by all R steps) | `07_common_functions.R` |

## Folder layout

```
07_snrnaseq/
├── README.md
├── 07_cellranger_multi.sh      # (0) Cell Ranger 10.0.0 multi, one sample per call
├── 07_qc_integration.R         # (1-4) QC → merge → PCA → Harmony → clustering
├── 07_annotation.R             # (5-10) markers, label transfer, manual annotation, proportions
├── 07_hdWGCNA.R                # (11) per-cell-type signed-hybrid co-expression networks
└── 07_common_functions.R       # shared helpers (logger, palettes, scattermore plotters, GMT reader)
```

## Script details

### `07_cellranger_multi.sh`
Runs `cellranger multi` for one sample per invocation (drive across samples with
a loop or scheduler array). Reads a per-sample config CSV at
`$WORKDIR/cellranger/${TAG}_config.csv` pointing to the reference transcriptome
and FASTQ paths. Produces the per-sample matrices consumed by `07_qc_integration.R`.

- **Tool:** Cell Ranger 10.0.0
- **Env vars:** `WORKDIR`, `CELLRANGER`, `TAG`
- **Output:** `$WORKDIR/cellranger/outputs/<sample>/...`

### `07_qc_integration.R`
Steps 1–4. Loads each Cell Ranger sample, applies per-nucleus QC, merges by batch
(disk-backed to bound peak memory), then normalizes, runs PCA, corrects batch with
Harmony, clusters (SNN, res 0.3), and builds UMAPs. Small clusters above index 18
are dropped, keeping 19 clusters (0–18) for annotation.

- **QC thresholds (cortex snRNA-seq):** nCount 1,000–30,000 · nFeature 750–8,000 · percent.mt < 5% · samples with median nFeature < 1,200 dropped
- **Tools:** Seurat v5, harmony
- **Input:** `$WORKDIR/cellranger/outputs/<sample>/`, `sample_batch_map.tsv`
- **Output:** `01_*_PRE.rds`, `01_*_QCed.rds`, `02_merged.rds`, `03_PCA_full.rds`, `04_Harmony_full_clean.rds`

### `07_annotation.R`
Steps 5–10. Per-cluster marker discovery, reference-based label transfer against
the Allen Brain Cell Atlas (FindTransferAnchors → MapQuery), prediction-score
thresholding (0 / 0.3 / 0.5) with QC visuals, manual mapping of 19 clusters to 10
cell types, and cell-type proportion analysis (Wilcoxon) across B6 / G1 / G2 and
the Hv / Hf / HL conditions.

- **Cell types (10):** CTX-CGE GABA, CTX-MGE GABA, CNU-LGE GABA, IT-ET Glut, NP-CT-L6b Glut, OB-CR Glut, Astro-Epen, Immune, OPC-Oligo, Vascular
- **Tools:** Seurat v5
- **Input:** `04_Harmony_full_clean.rds`, `$WORKDIR/allen/Allen_ref_<level>_min<n>.rds`
- **Output:** `05_annotation_full.rds`, marker tables, composition plots/stats

### `07_hdWGCNA.R`
Step 11. Per-cell-type signed-hybrid co-expression networks on vehicle nuclei
(WT B6 Wv + mutant G1/G2 Hv). Mouse genes are mapped to human orthologs, metacells
are built (k=25), soft power is selected (scale-free R² ≥ 0.80, slope < 0, power ≥ 5),
and modules (minModuleSize 50, mergeCutHeight 0.20) are summarized via eigengenes,
hub genes, module-trait correlation, and Fisher overlap with external gene sets.

- **Tools:** Seurat v5, hdWGCNA, WGCNA, ComplexHeatmap, enrichR (offline)
- **Input:** `05_annotation_full.rds`, `mouse2human_orthologs.rds`, `ref/genesets/*.gmt`
- **Output:** `$WORKDIR/seurat/full/annotation/hdWGCNA/<celltype>/...`

> GO enrichment (enrichR) is run separately in `07_hdWGCNA_enrichR_local.R`.

## Environment variables

| Variable | Used by | Description |
|----------|---------|-------------|
| `WORKDIR` | all | Project root (contains `cellranger/`, `seurat/`, `scripts/`, `allen/`, `ref/`) |
| `CELLRANGER` | `07_cellranger_multi.sh` | Cell Ranger install dir (contains the binary) |
| `TAG` | `07_cellranger_multi.sh` | Sample id for the run (one sample per call) |
| `TMP_DIR` | R steps | Scratch tmp directory (default `$WORKDIR/tmp`) |
| `R_LIBS_USER` | `07_hdWGCNA.R` | Optional; prepended to `.libPaths()` if set |

All scripts read paths from environment variables and default to the original
project tree, so no hard-coded user paths remain in source.

## Dependencies

- **Cell Ranger 10.0.0**
- **R ≥ 4.3** with: `Seurat` (v5), `harmony`, `hdWGCNA`, `WGCNA`, `ComplexHeatmap`,
  `circlize`, `enrichR`, `dplyr`, `tidyr`, `tibble`, `ggplot2`, `ggrepel`,
  `patchwork`, `scattermore`, `scales`, `RColorBrewer`, `pheatmap`, `Matrix`,
  `data.table`, `future`

## Execution order

```bash
# (0) per sample (loop / array)
for TAG in $(cat "$WORKDIR/cellranger/sample_list.txt"); do
  TAG="$TAG" bash 07_cellranger_multi.sh
done

# (1-4) → (5-10) → (11)
Rscript 07_qc_integration.R
Rscript 07_annotation.R
Rscript 07_hdWGCNA.R
```

`07_common_functions.R` is sourced automatically by each R step (not run directly).
