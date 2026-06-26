################################################################################
# 07. snRNA-seq annotation (label transfer + manual curation + proportions)
# Author : Yukyung Jun
# Updated: 2026-05  (merged marker/label-transfer with manual annotation;
#                    env-var paths, English comments, shared plotters in
#                    07_common_functions.R)
#
# Pipeline (full dataset):
#   Step 5  : FindAllMarkers - per-cluster marker discovery (annotation basis)
#   Step 6  : load Allen Brain Cell Atlas reference
#   Step 7  : reference-based label transfer (FindTransferAnchors -> MapQuery)
#   Step 8  : prediction-score thresholding + annotation QC visuals
#   Step 9  : manual cluster -> cell-type mapping (19 clusters -> 10 types)
#   Step 10 : cell-type proportion analysis (Wilcoxon; B6/G1/G2 and Hv/Hf/HL)
#
# Tools   : Seurat v5
#
# Required environment variables (defaults match the original project tree):
#   WORKDIR   - project root
#   TMP_DIR   - scratch tmp directory (default $WORKDIR/tmp)
# Input :  $WORKDIR/seurat/04_Harmony_full_clean.rds   (from 07_qc_integration.R)
#          $WORKDIR/allen/Allen_ref_<level>_min<n>.rds  (prebuilt reference)
# Output:  $WORKDIR/seurat/05_annotation_full.rds       (consumed by 07_hdWGCNA.R)
################################################################################

set.seed(777)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(data.table)
  library(Matrix)
  library(pheatmap)
  library(scattermore)
  library(future)
})

options(bitmapType = "cairo")
options(future.globals.maxSize = 1024 * 1024^3)
setDTthreads(threads = 8)

# =============================================================================
# Config (paths from environment; analysis parameters inline)
# =============================================================================
workdir   <- Sys.getenv("WORKDIR", unset = "/blues/scratch/yukyung/EJKim/Drug_Response")
tmp_dir   <- Sys.getenv("TMP_DIR", unset = file.path(workdir, "tmp"))
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_dir)

source(file.path(workdir, "scripts", "07_common_functions.R"))

outdir    <- file.path(workdir, "seurat")
full_dir  <- file.path(outdir, "full")
annot_dir <- file.path(full_dir, "annotation")
allendir  <- file.path(workdir, "allen")
dir.create(annot_dir, showWarnings = FALSE, recursive = TRUE)

# Allen reference selection
level     <- "class"
min_cells <- 50
key       <- sprintf("%s_min%d", level, min_cells)
lt_dir    <- file.path(full_dir, sprintf("label_transfer_%s", key))
allen_ref_path <- file.path(allendir, sprintf("Allen_ref_%s.rds", key))
dir.create(lt_dir, showWarnings = FALSE, recursive = TRUE)

thr_scores <- c(0, 0.3, 0.5)   # prediction-score thresholds to evaluate
config <- list(dims = 1:30)     # PCA/UMAP dims for the Allen reference build

write_log <- make_write_log(file.path(outdir, "analysis_log.txt"))
write_log("07 annotation started")

# =============================================================================
# Script-local helpers
# =============================================================================

# ---- annotation QC visualization for one score threshold (Step 8) -----------
run_visualize <- function(full_obj, out_dir, key, thr_score) {
  tag       <- sprintf("%s_thr%.2f", key, thr_score)
  ct_col    <- sprintf("celltype_%s",    key)
  lt_col    <- sprintf("lt_celltype_%s", key)
  score_col <- sprintf("lt_score_%s",    key)
  lc_col    <- sprintf("low_conf_%s",    tag)

  ct_in_data  <- unique(c(full_obj[[ct_col, drop = TRUE]], full_obj[[lt_col, drop = TRUE]]))
  plot_colors <- make_celltype_colors(ct_in_data)

  p_raw <- dimplot_scattermore(full_obj, group.by = lt_col, reduction = "umap",
                               colors = plot_colors,
                               title = sprintf("Label transfer raw [%s] - All", key))
  p_flt <- dimplot_scattermore(full_obj, group.by = ct_col, reduction = "umap",
                               colors = plot_colors,
                               title = sprintf("Annotation (score>=%.2f) [%s] - All", thr_score, key))
  png(file.path(out_dir, sprintf("UMAP_annotation_full_%s.png", tag)),
      width = 2500, height = 1100, res = 200)
  print(p_raw | p_flt); dev.off()

  png(file.path(out_dir, sprintf("Featureplot_pred_score_full_%s.png", tag)),
      width = 1000, height = 900, res = 200)
  print(score_scattermore(full_obj, score_col = score_col,
                          title = sprintf("Label transfer score [%s] - All", key)))
  dev.off()

  df_vln <- FetchData(full_obj, vars = c(score_col, ct_col))
  p_vln  <- ggplot(df_vln, aes(x = reorder(.data[[ct_col]], -.data[[score_col]], median),
                               y = .data[[score_col]])) +
    geom_violin(scale = "width", trim = TRUE, fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = thr_score, linetype = "dashed", color = "red") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "Prediction score",
         title = sprintf("Label transfer score by cell type [%s, thr=%.2f] - All", key, thr_score))
  ggsave(file.path(out_dir, sprintf("Vlnplot_pred_score_full_%s.png", tag)),
         p_vln, width = 9, height = 4, dpi = 200)

  full_obj[[lc_col]] <- full_obj[[score_col, drop = TRUE]] < thr_score
  full_obj[["lc_label_tmp"]] <- factor(
    ifelse(full_obj[[lc_col, drop = TRUE]], "Low confidence", "High confidence"))
  p_lc <- dimplot_scattermore(full_obj, group.by = "lc_label_tmp", reduction = "umap",
                              colors = c("High confidence" = "#4575b4", "Low confidence" = "#d73027"),
                              label = FALSE, pt.size = 1.0,
                              title = sprintf("Low confidence [%s, thr=%.2f] - All", key, thr_score))
  png(file.path(out_dir, sprintf("UMAP_low_confidence_full_%s.png", tag)),
      width = 1000, height = 900, res = 200)
  print(p_lc); dev.off()

  lowconf_stats <- full_obj@meta.data %>%
    group_by(across(all_of(lt_col))) %>%
    summarise(total = n(), low_conf = sum(.data[[lc_col]]),
              low_conf_frac = round(low_conf / total, 3),
              score_median  = round(median(.data[[score_col]]), 3), .groups = "drop") %>%
    arrange(desc(low_conf_frac))
  write.csv(lowconf_stats, file.path(out_dir, sprintf("lowconf_stats_full_%s.csv", tag)),
            row.names = FALSE)

  meta    <- full_obj@meta.data
  frac_df <- meta %>%
    count(batch, sample, .data[[ct_col]], name = "n") %>%
    rename(celltype = all_of(ct_col)) %>%
    group_by(batch, sample) %>% mutate(fraction = n / sum(n)) %>% ungroup()
  sample_order   <- frac_df %>% distinct(batch, sample) %>% arrange(batch, sample) %>% pull(sample)
  frac_df$sample <- factor(frac_df$sample, levels = sample_order)

  p_stack <- ggplot(frac_df, aes(x = sample, y = fraction, fill = celltype)) +
    geom_col(width = 0.95) +
    facet_grid(. ~ batch, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = plot_colors, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Fraction", fill = "Cell type",
         title = sprintf("Cell type composition by batch [%s] - All", tag)) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          strip.text = element_text(size = 10), panel.grid.major.x = element_blank())
  ggsave(file.path(out_dir, sprintf("celltype_fraction_stacked_batch_full_%s.png", tag)),
         p_stack, width = 20, height = 6, dpi = 200)

  p_stack_alpha <- frac_df %>%
    mutate(sample = factor(sample, levels = sort(unique(as.character(sample))))) %>%
    ggplot(aes(x = sample, y = fraction, fill = celltype)) +
    geom_col(width = 0.95) +
    scale_fill_manual(values = plot_colors, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Fraction", fill = "Cell type",
         title = sprintf("Cell type composition by sample [%s, thr=%.2f] - All", key, thr_score)) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          panel.grid.major.x = element_blank())
  ggsave(file.path(out_dir, sprintf("celltype_fraction_stacked_sample_full_%s.png", tag)),
         p_stack_alpha, width = 20, height = 6, dpi = 200)

  frac_wide <- frac_df %>% select(sample, batch, celltype, fraction) %>%
    pivot_wider(names_from = celltype, values_from = fraction, values_fill = 0)
  mat_hm <- as.matrix(frac_wide[, -(1:2)]); rownames(mat_hm) <- frac_wide$sample
  ann_row <- data.frame(batch = frac_wide$batch, row.names = frac_wide$sample)
  png(file.path(out_dir, sprintf("Heatmap_celltype_fraction_full_%s.png", tag)),
      width = 1400, height = 1000, res = 150)
  pheatmap(mat_hm, scale = "row", annotation_row = ann_row,
           fontsize_row = 7, fontsize_col = 9,
           main = sprintf("Cell type fraction (row-scaled) [%s, thr=%.2f] - All", key, thr_score))
  dev.off()

  stats_df <- frac_df %>%
    group_by(celltype) %>%
    summarise(p_value = tryCatch(kruskal.test(fraction ~ batch)$p.value,
                                 error = function(e) NA_real_), .groups = "drop") %>%
    mutate(FDR = p.adjust(p_value, method = "fdr")) %>% arrange(p_value)
  write.csv(stats_df, file.path(out_dir, sprintf("batch_fraction_kruskal_full_%s.csv", tag)),
            row.names = FALSE)

  full_obj
}

# ---- genotype name normalization and Group 1 / Group 2 definitions ----------
replacements <- c(
  "ADNP" = "Adnp",  "CHN"  = "Chd8-N",   "G2b"  = "Grin2b",
  "Arid" = "Arid1b","Tbl"  = "Tbl1xr1",  "Dyrk" = "Dyrk1a",
  "Na"   = "Naa15", "Scn"  = "Scn2a",    "ASH"  = "Ash1l",
  "Mt1"  = "Myt1l", "SdC"  = "Shank3-D", "T12"  = "Trip12",
  "K5B"  = "Kmt5b", "P131" = "Pten",     "Tanc" = "Tanc2",
  "TBR1" = "Tbr1",  "L68P" = "Shank3-L"
)
normalize_genotype <- function(raw_gt) {
  result <- raw_gt
  for (k in names(replacements)) result <- ifelse(result == k, replacements[[k]], result)
  result
}
g1 <- c('Naa15_M','Trip12_F','Shank3-D_F','Tanc2_M','Shank3-L_M',
        'Kmt5b_F','Kmt5b_M','Pten_M','Pten_F','Dyrk1a_M',
        'Arid1b_M','Shank3-D_M','Tanc2_F','Chd8-N_M','Ash1l_F')
g2 <- c('Adnp_F','Scn2a_M','Arid1b_F','Scn2a_F','Myt1l_M',
        'Myt1l_F','Dyrk1a_F','Shank3-L_F','Chd8-N_F','Trip12_M',
        'Grin2b_M','Adnp_M','Ash1l_M','Tbl1xr1_F','Grin2b_F',
        'Naa15_F','Tbl1xr1_M','Tbr1_M','Tbr1_F')

# ---- proportion-analysis helpers (unit_id-level fractions + Wilcoxon) -------
# NOTE: these reference the global `ct_colors`, defined in Step 9 before use.
calc_frac <- function(meta_df, grp_col) {
  meta_df %>%
    count(.data[[grp_col]], unit_id, celltype, name = "n") %>%
    group_by(unit_id) %>% mutate(fraction = n / sum(n)) %>% ungroup() %>%
    mutate(celltype = factor(celltype, levels = names(ct_colors)),
           !!grp_col := factor(.data[[grp_col]]))
}
calc_frac_mean <- function(frac_df, grp_col) {
  frac_df %>%
    group_by(.data[[grp_col]], celltype) %>%
    summarise(fraction = mean(fraction, na.rm = TRUE), .groups = "drop") %>%
    group_by(.data[[grp_col]]) %>% mutate(fraction = fraction / sum(fraction)) %>% ungroup()
}
plot_fraction_bar <- function(frac_df, grp_col, title = "", subtitle = "") {
  frac_df %>% mutate(unit_id = factor(unit_id)) %>%
    ggplot(aes(x = unit_id, y = fraction, fill = celltype)) +
    geom_bar(stat = "identity") +
    facet_grid(~ .data[[grp_col]], scales = "free_x", space = "free_x") +
    scale_fill_manual(values = ct_colors, name = "Cell type") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Fraction") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          strip.text = element_text(face = "bold", size = 11),
          legend.key.size = unit(0.4, "cm"),
          plot.subtitle = element_text(size = 8, color = "gray40"))
}
plot_fraction_mean <- function(frac_mean_df, x_var, title = "") {
  ggplot(frac_mean_df, aes(x = .data[[x_var]], y = fraction, fill = celltype)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = ct_colors, name = "Cell type") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = title, x = NULL, y = "Fraction") +
    theme_classic(base_size = 11) + theme(legend.key.size = unit(0.4, "cm"))
}
wilcox_pairwise <- function(df, grp_col, pairs) {
  bind_rows(lapply(pairs, function(p) {
    df %>% group_by(celltype) %>%
      summarise(
        p_value = tryCatch(wilcox.test(fraction[.data[[grp_col]] == p[1]],
                                       fraction[.data[[grp_col]] == p[2]],
                                       exact = FALSE)$p.value, error = function(e) NA_real_),
        mean_a = mean(fraction[.data[[grp_col]] == p[1]], na.rm = TRUE),
        mean_b = mean(fraction[.data[[grp_col]] == p[2]], na.rm = TRUE),
        n_a = sum(.data[[grp_col]] == p[1]), n_b = sum(.data[[grp_col]] == p[2]),
        .groups = "drop") %>%
      mutate(comparison = paste0(p[1], "_vs_", p[2]),
             FDR = p.adjust(p_value, method = "fdr"),
             signif = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                                p_value < 0.05 ~ "*", p_value < 0.1 ~ ".", TRUE ~ "ns"))
  })) %>% arrange(comparison, FDR)
}
plot_boxplot <- function(frac_df, grp_col, stat_df, grp_colors, title = "", subtitle = "") {
  stat_label <- stat_df %>%
    mutate(celltype = factor(celltype, levels = names(ct_colors)),
           label = paste0(sub("_vs_", " vs ", comparison), ": p=", signif(p_value, 2)),
           y_rank = as.integer(factor(comparison, levels = rev(unique(comparison)))))
  frac_df %>% mutate(celltype = factor(celltype, levels = names(ct_colors))) %>%
    ggplot(aes(x = .data[[grp_col]], y = fraction, fill = .data[[grp_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5) +
    geom_jitter(aes(color = .data[[grp_col]]), width = 0.15, size = 1.5, alpha = 0.8) +
    facet_wrap(~ celltype, scales = "free_y", ncol = 5) +
    geom_text(data = stat_label, aes(x = 2, y = Inf, label = label, vjust = y_rank * 1.8),
              inherit.aes = FALSE, size = 2.5, hjust = 0.5) +
    scale_fill_manual(values = grp_colors) + scale_color_manual(values = grp_colors) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Fraction") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none", strip.text = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 8, color = "gray40"))
}

# =============================================================================
# Load clustered object
# =============================================================================
obj_h <- readRDS(file.path(outdir, "04_Harmony_full_clean.rds"))

# =============================================================================
# Step 5: FindAllMarkers (per cluster)
# =============================================================================
write_log("=== Step 5: FindAllMarkers ===")
marker_dir <- file.path(full_dir, "cluster_markers")
dir.create(marker_dir, showWarnings = FALSE, recursive = TRUE)

DefaultAssay(obj_h) <- "RNA"
n_sample <- 300   # cap cells per cluster for speed

cells_sub <- obj_h@meta.data %>% rownames_to_column("cell") %>%
  group_by(seurat_clusters) %>% slice_sample(n = n_sample, replace = FALSE) %>% pull(cell)
obj_sub          <- subset(obj_h, cells = cells_sub)
obj_sub[["RNA"]] <- JoinLayers(obj_sub[["RNA"]])
Idents(obj_sub)  <- Idents(obj_h)[cells_sub]
write_log(sprintf("subset for markers: %d cells (<= %d/cluster)", ncol(obj_sub), n_sample))

plan(sequential)   # force non-parallel for reproducible FindAllMarkers
markers_all <- FindAllMarkers(obj_sub, assay = "RNA", test.use = "wilcox",
                              logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE,
                              max.cells.per.ident = n_sample, verbose = TRUE)
rm(obj_sub); gc()

if (is.null(markers_all) || nrow(markers_all) == 0) {
  write_log("WARNING: FindAllMarkers returned no results; skipping Step 5b plots")
} else {
  markers_sig <- markers_all %>% filter(p_val_adj < 0.05) %>% arrange(cluster, desc(avg_log2FC))
  write.csv(markers_all, file.path(marker_dir, "FindAllMarkers_full_all.csv"), row.names = FALSE)
  write.csv(markers_sig, file.path(marker_dir, "FindAllMarkers_full_sig.csv"), row.names = FALSE)
  for (cl in unique(markers_sig$cluster))
    write.csv(markers_sig %>% filter(cluster == cl),
              file.path(marker_dir, sprintf("markers_cluster%s.csv", cl)), row.names = FALSE)
  write.csv(markers_sig %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup(),
            file.path(marker_dir, "FindAllMarkers_full_top10.csv"), row.names = FALSE)
  write_log(sprintf("FindAllMarkers done: total=%d, sig=%d", nrow(markers_all), nrow(markers_sig)))

  top5_genes <- markers_sig %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5) %>%
    pull(gene) %>% unique()
  if (length(top5_genes)) {
    p <- DotPlot(obj_h, features = top5_genes) + RotatedAxis() +
      ggtitle("Top 5 markers per cluster") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
    ggsave(file.path(marker_dir, "Dotplot_top5_markers_per_cluster.png"),
           p, width = max(12, length(top5_genes) * 0.35), height = 7, dpi = 200)
  }

  top3_genes <- markers_sig %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 3) %>%
    pull(gene) %>% unique()
  if (length(top3_genes)) {
    cells_hm <- obj_h@meta.data %>% rownames_to_column("cell") %>%
      group_by(seurat_clusters) %>% slice_sample(n = n_sample, replace = FALSE) %>% pull(cell)
    obj_hm          <- subset(obj_h, cells = cells_hm)
    obj_hm[["RNA"]] <- JoinLayers(obj_hm[["RNA"]])
    p_hm <- DoHeatmap(obj_hm, features = top3_genes, size = 3) +
      ggtitle("Top 3 markers per cluster (downsampled)") +
      theme(axis.text.y = element_text(size = 6))
    ggsave(file.path(marker_dir, "Heatmap_top3_markers_per_cluster.png"),
           p_hm, width = 16, height = max(8, length(top3_genes) * 0.15), dpi = 200)
    rm(obj_hm); gc()
  }

  top1_genes <- markers_sig %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 1) %>%
    pull(gene) %>% unique() %>% intersect(rownames(obj_h))
  if (length(top1_genes)) {
    plots_top1 <- lapply(top1_genes, feature_scattermore, obj = obj_h)
    ncol_t <- min(5, length(top1_genes)); nrow_t <- ceiling(length(top1_genes) / ncol_t)
    png(file.path(marker_dir, "Featureplot_top1_markers_per_cluster.png"),
        width = 500 * ncol_t, height = 500 * nrow_t, res = 150)
    print(wrap_plots(plots_top1, ncol = ncol_t)); dev.off()
  }
}
write_log("Step 5 done")

# =============================================================================
# Step 6: load Allen reference
# =============================================================================
write_log(sprintf("=== Step 6: load Allen reference [%s] ===", key))
if (!file.exists(allen_ref_path)) stop("Allen reference not found: ", allen_ref_path)

# =============================================================================
# Step 7: MapQuery label transfer
# =============================================================================
write_log(sprintf("=== Step 7: MapQuery annotation [%s] ===", key))
ref <- readRDS(allen_ref_path)
write_log(sprintf("[%s] Allen ref loaded: cells=%d", key, ncol(ref)))

if (inherits(ref[["RNA"]], "Assay5")) ref <- JoinLayers(ref, assay = "RNA")
DefaultAssay(ref) <- "RNA"; DefaultAssay(obj_h) <- "RNA"

ref   <- FindVariableFeatures(ref,   nfeatures = 3000, selection.method = "vst", verbose = FALSE)
obj_h <- FindVariableFeatures(obj_h, nfeatures = 3000, selection.method = "vst", verbose = FALSE)
transfer_features <- intersect(VariableFeatures(ref), VariableFeatures(obj_h))
if (length(transfer_features) > 3000) transfer_features <- transfer_features[seq_len(3000)]
write_log(sprintf("[%s] transfer features: %d", key, length(transfer_features)))
if (length(transfer_features) < 200) stop(sprintf("[%s] HVG intersection too small: %d", key, length(transfer_features)))

ref <- ref[transfer_features, ]
ref <- NormalizeData(ref, verbose = FALSE)
ref <- ScaleData(ref, features = transfer_features, verbose = FALSE)
ref <- RunPCA(ref, features = transfer_features, npcs = max(config$dims), verbose = FALSE)
ref <- RunUMAP(ref, dims = config$dims, reduction = "pca", return.model = TRUE, verbose = FALSE)

anchors <- FindTransferAnchors(reference = ref, query = obj_h,
                               normalization.method = "LogNormalize",
                               reduction = "rpca", reference.reduction = "pca",
                               features = transfer_features, verbose = TRUE)
obj_h <- MapQuery(anchorset = anchors, query = obj_h, reference = ref,
                  refdata = list(celltype = "celltype"),
                  reference.reduction = "pca", reduction.model = "umap", verbose = TRUE)

lt_col    <- sprintf("lt_celltype_%s", key)
score_col <- sprintf("lt_score_%s",    key)
obj_h[[lt_col]]    <- obj_h$predicted.celltype
obj_h[[score_col]] <- obj_h$predicted.celltype.score
write_log(sprintf("[%s] predictions extracted", key))

# reference + query UMAP overlays
umap_ref <- Embeddings(ref, "umap"); umap_qry <- Embeddings(obj_h, "ref.umap")
df_ref <- data.frame(umap_ref, celltype = ref$celltype, origin = "reference")
df_qry <- data.frame(umap_qry, celltype = obj_h$predicted.celltype, origin = "query")
colnames(df_ref)[1:2] <- colnames(df_qry)[1:2] <- c("UMAP1", "UMAP2")
lt_ct_colors <- make_celltype_colors(unique(c(df_ref$celltype, df_qry$celltype)))

p_ref_qry <- ggplot() +
  geom_scattermore(data = df_ref, aes(x = UMAP1, y = UMAP2, color = celltype),
                   pointsize = 1.5, pixels = c(1000, 1000), alpha = 0.3) +
  geom_scattermore(data = df_qry, aes(x = UMAP1, y = UMAP2, color = celltype),
                   pointsize = 1.5, pixels = c(1000, 1000), alpha = 0.8) +
  scale_color_manual(values = lt_ct_colors) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
  ggtitle(sprintf("Reference (faint) + Query (bold) UMAP [%s] - All", key)) +
  theme_classic(base_size = 11) +
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
ggsave(file.path(lt_dir, sprintf("UMAP_ref_query_full_%s.png", key)), p_ref_qry, width = 7, height = 6, dpi = 200)

p_qry_only <- ggplot(df_qry, aes(x = UMAP1, y = UMAP2, color = celltype)) +
  geom_scattermore(pointsize = 1.5, pixels = c(1000, 1000), alpha = 0.8) +
  scale_color_manual(values = lt_ct_colors) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  ggtitle(sprintf("Query projected onto Reference UMAP [%s] - All", key)) +
  theme_classic(base_size = 11) +
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
ggsave(file.path(lt_dir, sprintf("UMAP_query_projected_full_%s.png", key)), p_qry_only, width = 7, height = 6, dpi = 200)

p_ref_only <- ggplot(df_ref, aes(x = UMAP1, y = UMAP2, color = celltype)) +
  geom_scattermore(pointsize = 1.5, pixels = c(1000, 1000), alpha = 0.4) +
  scale_color_manual(values = lt_ct_colors) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  ggtitle(sprintf("Reference UMAP [%s] - All", key)) +
  theme_classic(base_size = 11) +
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
ggsave(file.path(lt_dir, sprintf("UMAP_ref_only_full_%s.png", key)), p_ref_only, width = 7, height = 6, dpi = 200)

raw_tab   <- sort(table(obj_h[[lt_col, drop = TRUE]]), decreasing = TRUE)
score_sum <- summary(obj_h[[score_col, drop = TRUE]])
summary_path <- file.path(lt_dir, sprintf("lt_summary_full_%s.txt", key))
sink(summary_path)
cat(sprintf("=== MapQuery annotation (raw) [%s] ===\n", key)); print(raw_tab)
cat(sprintf("\n=== Score summary [%s] ===\n", key));            print(score_sum)
sink()
write_log(sprintf("[%s] Step 7 done; summary: %s", key, summary_path))

# =============================================================================
# Step 8: score thresholding + annotation QC visuals
# =============================================================================
write_log("=== Step 8: score thresholding & visualization ===")
ct_col <- sprintf("celltype_%s", key)

for (thr_score in thr_scores) {
  obj_h[[ct_col]] <- ifelse(obj_h[[score_col, drop = TRUE]] >= thr_score,
                            obj_h[[lt_col, drop = TRUE]], "Unassigned")
  out_dir <- file.path(lt_dir, sprintf("thr%.2f", thr_score))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  obj_h <- run_visualize(obj_h, out_dir, key, thr_score)
  write_log(sprintf("[%s, thr=%.2f] visualization done", key, thr_score))
}
write_log("Step 8 done (label transfer)")

# =============================================================================
# Step 9: manual cluster -> cell-type annotation
# =============================================================================
write_log("=== Step 9: manual annotation ===")

# parse sample name -> genotype / condition / replicate / sex / group / unit_id
obj_h@meta.data <- obj_h@meta.data %>%
  mutate(
    raw_genotype = sub("^(.+)_([^_]+)_([^_]+)$", "\\1", sample),
    condition    = sub("^(.+)_([^_]+)_([^_]+)$", "\\2", sample),
    replicate    = sub("^(.+)_([^_]+)_([^_]+)$", "\\3", sample),
    genotype     = normalize_genotype(raw_genotype),
    sample       = paste(genotype, condition, replicate, sep = "_"),
    sex          = case_when(grepl("F$", condition) ~ "F",
                             grepl("M$", condition) ~ "M",
                             TRUE ~ NA_character_),
    genotype_sex = paste(genotype, sex, sep = "_"),
    group        = case_when(genotype_sex %in% g1 ~ "G1",
                             genotype_sex %in% g2 ~ "G2",
                             TRUE ~ "B6"),
    unit_id      = case_when(group == "B6" ~ sample, TRUE ~ genotype_sex)
  )

# manual cluster -> cell-type mapping (19 clusters -> 10 cell types)
cluster_map <- c(
  "0"  = "IT-ET Glut",   "1"  = "IT-ET Glut",   "2"  = "CTX-CGE GABA",
  "3"  = "IT-ET Glut",   "4"  = "Astro-Epen",   "5"  = "CTX-CGE GABA",
  "6"  = "IT-ET Glut",   "7"  = "IT-ET Glut",   "8"  = "CTX-MGE GABA",
  "9"  = "Vascular",     "10" = "Immune",       "11" = "OPC-Oligo",
  "12" = "NP-CT-L6b Glut","13" = "OPC-Oligo",   "14" = "Astro-Epen",
  "15" = "CNU-LGE GABA", "16" = "Vascular",     "17" = "OB-CR Glut",
  "18" = "Astro-Epen"
)
obj_h@meta.data$celltype <- cluster_map[as.character(obj_h@meta.data$seurat_clusters)]

all_ct_colors <- c(
  "CTX-CGE GABA" = "#A6CEE3", "CTX-MGE GABA" = "#1F78B4", "CNU-LGE GABA" = "#33A02C",
  "IT-ET Glut"   = "#FDBF6F", "NP-CT-L6b Glut" = "#FF7F00", "OB-CR Glut" = "#E31A1C",
  "Astro-Epen"   = "#CAB2D6", "Immune" = "#6A3D9A", "OPC-Oligo" = "#FFFF99",
  "Vascular"     = "#B15928"
)
present_ct <- unique(na.omit(obj_h@meta.data$celltype))
ct_levels  <- intersect(names(all_ct_colors), present_ct)
ct_colors  <- all_ct_colors[ct_levels]               # used by proportion helpers
obj_h$celltype <- factor(obj_h$celltype, levels = ct_levels)

dot_features <- intersect(
  c("Gad2", "Dlx1", "Sp9", "Lhx6", "Pvalb", "Drd1", "Isl1", "Slc17a7",
    "Satb2", "Tfap2d", "Tbr1", "Lhx5", "Aqp4", "Cx3cr1", "Tmem119",
    "Pdgfra", "Mog", "Cldn5"),
  rownames(obj_h))

meta <- obj_h@meta.data
write.csv(meta, file.path(annot_dir, "meta.csv"), row.names = TRUE)
Idents(obj_h) <- obj_h$celltype

ggsave(file.path(annot_dir, "Vlnplot_markers_all_celltypes.pdf"),
       VlnPlot(obj_h, features = dot_features, idents = ct_levels, stack = TRUE, flip = TRUE,
               fill.by = "ident", cols = ct_colors, pt.size = 0) +
         theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
               axis.text.y = element_text(size = 8), strip.text = element_text(size = 8),
               legend.position = "none"),
       width = 8, height = max(10, length(dot_features) * 0.6), device = "pdf")

p_dot <- DotPlot(obj_h, features = dot_features, idents = ct_levels,
                 scale = TRUE, col.min = -2, col.max = 2) +
  coord_flip() +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 0, limits = c(-2, 2), name = "Average\nExpression") +
  scale_size(range = c(0, 6), name = "Percent\nExpressed") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        panel.grid = element_blank(),
        legend.title = element_text(size = 9), legend.text = element_text(size = 8),
        plot.title = element_text(size = 11, hjust = 0.5))
ggsave(file.path(annot_dir, "Dotplot_markers_all_celltypes.pdf"), p_dot,
       width = max(6, length(ct_levels) * 0.4), height = max(7, length(dot_features) * 0.35),
       device = "pdf")

ggsave(file.path(annot_dir, "UMAP_celltype_all.pdf"),
       dimplot_scattermore(obj_h, group.by = "celltype", reduction = "umap",
                           colors = ct_colors, pt.size = 1.2, pixels = c(1000, 1000)),
       width = 6, height = 5, device = "pdf")
write_log("Step 9 done: cell types assigned")

# =============================================================================
# Step 10: per-sample QC + cell-type proportion analysis
# =============================================================================
write_log("=== Step 10: proportion analysis ===")

cell_count  <- meta %>% count(sample, name = "n_cells") %>% arrange(desc(n_cells))
meta_sorted <- meta %>% mutate(sample = factor(sample, levels = cell_count$sample))

qc_df <- meta_sorted %>%
  group_by(sample) %>%
  summarise(Count_RNA = mean(nCount_RNA, na.rm = TRUE),
            nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
            percent.mt = mean(percent.mt, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Count_RNA", "nFeature_RNA", "percent.mt")))

p_count <- ggplot(cell_count %>% mutate(sample = factor(sample, levels = cell_count$sample)),
                  aes(x = sample, y = n_cells)) +
  geom_bar(stat = "identity", fill = "grey40") +
  geom_text(aes(label = scales::comma(n_cells)), vjust = -0.3, size = 2) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(y = "Cell count") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())

p_qc <- ggplot(qc_df, aes(x = sample, y = value)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  facet_wrap(~ metric, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), strip.placement = "outside",
        strip.text = element_text(size = 9))

p_frac_all <- meta_sorted %>%
  count(sample, celltype, name = "n") %>%
  group_by(sample) %>% mutate(fraction = n / sum(n)) %>% ungroup() %>%
  mutate(celltype = factor(celltype, levels = names(ct_colors))) %>%
  ggplot(aes(x = sample, y = fraction, fill = celltype)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = ct_colors, drop = FALSE, name = "Cell type") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Fraction", title = "Cell type fraction by sample") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        panel.grid.major.x = element_blank(), legend.key.size = unit(0.4, "cm"))

ggsave(file.path(annot_dir, "Fraction_sample_all_qc.pdf"),
       p_count / p_qc / p_frac_all + plot_layout(heights = c(1, 3, 3)),
       width = 20, height = 10, device = "pdf")

ggsave(file.path(annot_dir, "Fraction_sample_batch_all.pdf"),
       meta %>%
         count(batch, sample, celltype, name = "n") %>%
         group_by(sample) %>% mutate(fraction = n / sum(n)) %>% ungroup() %>%
         mutate(celltype = factor(celltype, levels = names(ct_colors))) %>%
         ggplot(aes(x = sample, y = fraction, fill = celltype)) +
         geom_col(width = 0.9) +
         facet_grid(. ~ batch, scales = "free_x", space = "free_x") +
         scale_fill_manual(values = ct_colors, drop = FALSE, name = "Cell type") +
         scale_y_continuous(labels = percent_format(accuracy = 1)) +
         labs(x = NULL, y = "Fraction", title = "Cell type fraction by sample (facet: batch)") +
         theme_bw(base_size = 11) +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
               strip.text = element_text(size = 9, face = "bold"),
               panel.grid.major.x = element_blank(), legend.key.size = unit(0.4, "cm")),
       width = 20, height = 4.5, device = "pdf")

# ---- proportion analysis 1: vehicle B6 vs G1 vs G2 --------------------------
annot_dir_hv <- file.path(annot_dir, "Hv_G1vsG2")
dir.create(annot_dir_hv, recursive = TRUE, showWarnings = FALSE)
GRP_COLORS_HV <- c("B6" = "gray", "G1" = "yellowgreen", "G2" = "orange")
GRP_PAIRS_HV  <- list(c("B6", "G1"), c("B6", "G2"), c("G1", "G2"))
SUBTITLE_HV   <- "B6: per sample  |  G1/G2: per genotype_sex"

meta_hv <- meta %>%
  filter((grepl("Hv", condition) & group %in% c("G1", "G2")) |
           (grepl("Wv", condition) & group == "B6")) %>%
  mutate(group = factor(group, levels = c("B6", "G1", "G2")))
frac_hv      <- calc_frac(meta_hv, "group")
frac_hv_mean <- calc_frac_mean(frac_hv, "group")
stat_hv      <- wilcox_pairwise(frac_hv, "group", GRP_PAIRS_HV)
write.csv(stat_hv, file.path(annot_dir_hv, "Wilcoxon_G1vsG2_Hv.csv"), row.names = FALSE)

ggsave(file.path(annot_dir_hv, "1_Fraction_Hv_G1vsG2.pdf"),
       plot_fraction_bar(frac_hv, "group", "Cell type fraction - Hv (G1 vs G2)", SUBTITLE_HV),
       width = 18, height = 4.5, device = "pdf")
ggsave(file.path(annot_dir_hv, "2_Fraction_Hv_G1vsG2_stat.pdf"),
       plot_boxplot(frac_hv, "group", stat_hv, GRP_COLORS_HV,
                    "Cell type fraction - Hv (Wilcoxon)", SUBTITLE_HV),
       width = 10, height = 6, device = "pdf")
ggsave(file.path(annot_dir_hv, "3_Fraction_group_mean_Hv_G1vsG2.pdf"),
       plot_fraction_mean(frac_hv_mean, "group", "Cell type fraction by group mean - Hv (G1 vs G2)"),
       width = 5, height = 4.5, device = "pdf")

# ---- proportion analysis 2: within G1 / G2, Hv vs Hf vs HL ------------------
annot_dir_cond <- file.path(annot_dir, "Het_condition")
dir.create(annot_dir_cond, recursive = TRUE, showWarnings = FALSE)
COND_LEVELS <- c("Hv", "Hf", "HL")
COND_COLORS <- c("Hv" = "#4E79A7", "Hf" = "#F28E2B", "HL" = "#E15759")
COND_PAIRS  <- list(c("Hv", "Hf"), c("Hv", "HL"))
SUBTITLE_HT <- "HT: per genotype_sex"

for (grp in c("G1", "G2")) {
  meta_grp <- meta %>%
    filter(group == grp, grepl("^H", condition)) %>%
    mutate(condition = sub("[MF]$", "", condition),
           condition = factor(condition, levels = COND_LEVELS),
           unit_id = paste(genotype_sex, condition, sep = "_")) %>%
    filter(!is.na(condition))
  write_log(sprintf("[%s] Hv=%d, Hf=%d, HL=%d cells", grp,
                    sum(meta_grp$condition == "Hv"), sum(meta_grp$condition == "Hf"),
                    sum(meta_grp$condition == "HL")))

  frac_grp      <- calc_frac(meta_grp, "condition")
  frac_grp_mean <- calc_frac_mean(frac_grp, "condition")
  stat_grp      <- wilcox_pairwise(frac_grp, "condition", COND_PAIRS)
  write.csv(stat_grp, file.path(annot_dir_cond, sprintf("Wilcoxon_condition_%s.csv", grp)),
            row.names = FALSE)

  ggsave(file.path(annot_dir_cond, sprintf("1_Fraction_%s.pdf", grp)),
         plot_fraction_bar(frac_grp, "condition",
                           sprintf("Cell type fraction - %s (Hv/Hf/HL)", grp), SUBTITLE_HT),
         width = 20, height = 4.5, device = "pdf")
  ggsave(file.path(annot_dir_cond, sprintf("2_Fraction_%s_stat.pdf", grp)),
         plot_boxplot(frac_grp, "condition", stat_grp, COND_COLORS,
                      sprintf("Cell type fraction - %s (Wilcoxon)", grp), SUBTITLE_HT),
         width = 10, height = 6, device = "pdf")
  ggsave(file.path(annot_dir_cond, sprintf("3_Fraction_mean_%s.pdf", grp)),
         plot_fraction_mean(frac_grp_mean, "condition",
                            sprintf("Cell type fraction by condition mean - %s", grp)),
         width = 5, height = 4.5, device = "pdf")
}
write_log("Step 10 done")

# =============================================================================
# Save + reproducibility
# =============================================================================
saveRDS(obj_h, file.path(outdir, "05_annotation_full.rds"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_07_annotation.txt"))
write_log("=== 07 annotation complete: 05_annotation_full.rds ===")
cat("[OK] 07 annotation done ->", outdir, "\n")
