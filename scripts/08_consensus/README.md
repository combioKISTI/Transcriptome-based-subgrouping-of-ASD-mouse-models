# 08_consensus — Consensus and cross-species concordance gene analysis


## Analysis steps

| Method step | Description | Scripts |
|-------------|-------------|---------|
| **(0)** Input preparation | Per-gene log2FC matrix per species × group × category | (upstream DEG analysis) |
| **(1)** ID mapping (human only) | Ensembl gene ID → HGNC symbol | `scripts/01_id_mapping/ensembl_to_symbol.py` |
| **(2)** Group split + category filter | (a) Reorder columns into [G1 | G2] (b) Filter rows by SynGO / chromatin / mRNA gene list | `scripts/02_prepare_input/{split_columns_by_group.sh, filter_by_category.py}` |
| **(3)** One-sample one-sided Wilcoxon + BH | For each gene, run `wilcox.test(x, mu=0, alternative="less")` and `..., alternative="greater"`; BH-adjust separately; call WilcoxonDown / WilcoxonUp / NoChange | `scripts/03_wilcoxon/wilcoxon_one_sided_BH.R` (heatmap + sorted TSV) plus `extract_calls.py` |
| **(4–5)** Antiparallel consensus definition | Synaptic: G1-Down ∩ G2-Up; chromatin & mRNA: G1-Up ∩ G2-Down | `scripts/04_antiparallel_consensus/overlap_antiparallel_genes.py` and four supporting scripts |
| **(6)** Hypergeometric test | `phyper(x-1, K, N-K, n, lower.tail=FALSE)` for the antiparallel overlap; background = WGCNA-input expressed genes | `scripts/05_hypergeometric/hypergeometric_test.R` |
| **(7)** Cross-species concordance | (a) Ensembl BioMart release 102 one-to-one ortholog mapping via HGNC symbols (b) Intersect mouse antiparallel ∩ human antiparallel with matched directionality (mouse G1 ↔ human hGroup 1) | `scripts/06_cross_species/{mouse_human_ortholog.R, concordant_genes.py}` |

## Folder layout

```
08_consensus/
├── README.md
├── data/
│   ├── list_genes_syngo.txt         # SynGO release 2023-12-01 (synaptic curation)
│   ├── list_genes_chromatin.txt     # MSigDB GOBP "chromatin organization" + descendants
│   └── list_genes_mRNA.txt          # MSigDB GOBP "mRNA processing" + descendants
├── scripts/
│   ├── 01_id_mapping/
│   │   └── ensembl_to_symbol.py
│   ├── 02_prepare_input/
│   │   ├── split_columns_by_group.sh
│   │   └── filter_by_category.py
│   ├── 03_wilcoxon/
│   │   ├── wilcoxon_one_sided_BH.R          # main: heatmap + labeled sorted TSV
│   │   └── extract_calls.py                 # extract WilcoxonDown/Up call table
│   ├── 04_antiparallel_consensus/
│   │   ├── overlap_antiparallel_genes.py    # define G1-call ∩ G2-call genes
│   │   ├── filter_input_by_overlapped.py    # subset matrix to overlapped genes
│   │   ├── extract_overlapped_stats.py      # subset stats per group
│   │   ├── merge_g1_g2.py                   # combine G1+G2 stats with suffixes
│   │   └── heatmap_cluster.R                # 2-way clustering heatmap
│   ├── 05_hypergeometric/
│   │   └── hypergeometric_test.R            # generic phyper with N K n x args
│   └── 06_cross_species/
│       ├── mouse_human_ortholog.R           # BioMart release 102, one-to-one
│       └── concordant_genes.py              # mouse↔human antiparallel intersection
├── pipelines/
│   ├── run_mouse_consensus.sh               # mouse end-to-end (single category)
│   ├── run_human_consensus.sh               # human end-to-end (with ID mapping)
│   ├── run_hypergeometric.sh                # batch phyper for all categories
│   └── run_cross_species.sh                 # mouse-human concordant extraction
└── results/                                  # output directory (gitignore recommended)
```

## Input file format

Each per-species × per-category TSV must follow this layout:

```
group_label    Group1   Group1   ...   Group2   Group2   ...
sample_id      g1_s1    g1_s2    ...   g2_s1    g2_s2    ...
GeneA          0.43     -0.12    ...   -0.55    -0.31    ...
GeneB          ...
```

- Row 1: group label (`Group1` / `Group2`)
- Row 2: sample identifier (mouse: genotype × sex; human: subject ID)
- Row 3+: gene symbol (column 1) and per-sample log2 fold-change

Sample sizes used in the manuscript: mouse G1 = 15, mouse G2 = 19; human hGroup 1 = 20, human hGroup 2 = 20.

## Antiparallel direction by category

| Category | G1 call | G2 call | Curation source |
|----------|---------|---------|-----------------|
| Synaptic | WilcoxonDown | WilcoxonUp | SynGO release 2023-12-01 |
| Chromatin | WilcoxonUp | WilcoxonDown | MSigDB v2023.2.Hs GOBP "chromatin organization" + descendants |
| mRNA processing | WilcoxonUp | WilcoxonDown | MSigDB v2023.2.Hs GOBP "mRNA processing" + descendants |

## Antiparallel consensus counts (manuscript values)

| Category | Mouse G1∩G2 | Human hG1∩hG2 | Cross-species concordant |
|----------|-------------|----------------|--------------------------|
| Synaptic | 93 | 195 | 14 |
| Chromatin | 51 | 23 | — |
| mRNA processing | 20 | 27 | — |

## Hypergeometric tests (manuscript values, human BA9)

`phyper(x-1, K, N-K, n, lower.tail=FALSE)`:

| Category (human) | N (background) | K (G2 opposite-direction) | n (G1 one-direction) | x (overlap) | p-value |
|------------------|----------------|---------------------------|----------------------|-------------|---------|
| Synaptic | 1,555 | 218 | 806 | 195 | (compute via script) |
| Chromatin | 961 | 64 | 98 | 23 | 6.507 × 10⁻⁹ |
| mRNA processing | 947 | 66 | 74 | 27 | 2.660 × 10⁻¹⁵ |

The background N is the gene universe retained in the WGCNA input matrix.

## Quick start

```bash
cd 08_consensus

# 1. Mouse pipeline (synaptic example)
bash pipelines/run_mouse_consensus.sh synaptic data/gene_fc_mouse_group.txt

# 2. Human pipeline (synaptic, includes Ensembl→symbol mapping if needed)
bash pipelines/run_human_consensus.sh synaptic data/gene_fc_mapped_human_group.txt

# 3. Hypergeometric tests (uses values reported in the manuscript by default)
bash pipelines/run_hypergeometric.sh

# 4. Cross-species concordant genes
bash pipelines/run_cross_species.sh synaptic
```

## Dependencies

- **R ≥ 4.3** with packages: `ComplexHeatmap`, `circlize`, `grid`, `biomaRt`, `dplyr`
- **Python ≥ 3.8** with `pandas`, `requests`
- Standard Unix tools: `bash`, `awk`, `tr`, `sort`, `comm`

## Reproducibility

The pipeline has been dry-run end-to-end against `data/gene_fc_mouse_group.txt`
and `data/gene_fc_mapped_human_group.txt`, reproducing the manuscript values
exactly:

| Run | Antiparallel count | Manuscript |
|-----|-------------------|------------|
| Mouse synaptic | 93 | 93 |
| Mouse chromatin | 51 | 51 |
| Mouse mRNA | 20 | 20 |
| Human synaptic | 195 | 195 |
| Human chromatin | 23 | 23 |
| Human mRNA | 27 | 27 |
| Cross-species synaptic | 14 | 14 |
