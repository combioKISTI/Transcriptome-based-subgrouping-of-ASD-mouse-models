################################################################################
# 07. snRNA-seq QC, merge, integration, and clustering
# Author : Yukyung Jun
# Updated: 2026-05  (merged QC/merge with normalization/Harmony/clustering;
#                    env-var paths, English comments)
#
# Pipeline:
#   Step 1 : load each Cell Ranger sample -> QC filter
#            (nCount / nFeature / percent.mt) -> save PRE + QCed RDS
#   Step 2 : drop low-quality samples -> per-batch merge -> 02_merged.rds
#   Step 3 : NormalizeData -> 3,000 HVG -> ScaleData -> PCA -> 03_PCA_full.rds
#   Step 4 : Harmony (batch) -> SNN clustering (res 0.3) -> UMAP
#            -> drop small clusters (>18, too few cells) -> parse metadata
#            -> 04_Harmony_full_clean.rds
#
# Upstream: Cell Ranger 10.0.0 `multi` (see 07_cellranger_multi.sh)
# Tools   : Seurat v5, harmony
#
# Required environment variables (defaults point to the original project tree):
#   WORKDIR   - project root (contains cellranger/, seurat/, scripts/)
#   TMP_DIR   - scratch tmp directory   (default $WORKDIR/tmp)
# Inputs:
#   $WORKDIR/cellranger/outputs/<sample>/           (10x matrices)
#   $WORKDIR/cellranger/sample_batch_map.tsv        (Batch, Sample.Name)
################################################################################

set.seed(777)
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(data.table)
})

options(bitmapType = "cairo")
options(future.globals.maxSize = 1024 * 1024^3)

# =============================================================================
# Config (paths from environment; analysis parameters inline)
# =============================================================================
workdir     <- Sys.getenv("WORKDIR", unset = "/blues/scratch/yukyung/EJKim/Drug_Response")
tmp_dir     <- Sys.getenv("TMP_DIR", unset = file.path(workdir, "tmp"))
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_dir)

source(file.path(workdir, "scripts", "07_common_functions.R"))

outdir      <- file.path(workdir, "seurat")
qc_dir      <- file.path(outdir, "01_QC")
qc_plot_dir <- file.path(qc_dir, "plots")
dir.create(qc_plot_dir, showWarnings = FALSE, recursive = TRUE)

config <- list(
  # ---- QC thresholds (snRNA-seq cortex) -------------------------------------
  # nCount   1000-30000 : cortex mean UMI ~3,000-8,000; upper bound trims
  #                       doublets / debris
  # nFeature  750-8000  : cortical neurons express ~2,000-5,000 genes
  # percent.mt  < 5%    : nuclei carry little mitochondrial RNA
  min_count        = 1000,
  max_count        = 30000,
  min_features     = 750,
  max_features     = 8000,
  mito_threshold   = 5,
  # ---- normalization / dimensionality reduction -----------------------------
  n_variable_features   = 3000,
  npcs                  = 30,
  # ---- clustering (Step 4) --------------------------------------------------
  resolution            = 0.3,
  n_neighbors           = 20,
  min_dist              = 0.2,
  # Samples whose median nFeature falls below this are dropped before merge.
  min_median_nfeature   = 1200
)

cellranger_dir <- file.path(workdir, "cellranger", "outputs")
sample_names   <- list.files(cellranger_dir, full.names = FALSE)
data_dirs      <- file.path(cellranger_dir, sample_names)

# ---- batch mapping ----------------------------------------------------------
batch_map_df <- read.delim(file.path(workdir, "cellranger", "sample_batch_map.tsv"),
                           sep = "\t", header = TRUE, stringsAsFactors = FALSE)
if ("Sample Name" %in% names(batch_map_df))
  names(batch_map_df)[names(batch_map_df) == "Sample Name"] <- "Sample.Name"
stopifnot(all(c("Batch", "Sample.Name") %in% names(batch_map_df)))
sample_to_batch <- setNames(batch_map_df$Batch, batch_map_df$Sample.Name)

missing <- setdiff(sample_names, names(sample_to_batch))
if (length(missing)) stop("batch mapping missing for: ", paste(missing, collapse = ", "))
for (d in data_dirs) if (!dir.exists(d)) stop("directory not found: ", d)

# ---- logging ----------------------------------------------------------------
write_log <- make_write_log(file.path(outdir, "analysis_log.txt"))
write_log("07 QC + integration started")
write_log(sprintf("QC thresholds - nCount: [%d, %d] | nFeature: [%d, %d] | mito < %g%%",
                  config$min_count, config$max_count,
                  config$min_features, config$max_features, config$mito_threshold))

# =============================================================================
# Script-local helpers
# =============================================================================

# ---- per-sample QC summary row ----------------------------------------------
qc_summary <- function(meta_df, sample, state) {
  data.frame(
    sample          = sample, state = state,
    cells           = nrow(meta_df),
    nCount_mean     = mean(meta_df$nCount_RNA),
    nCount_median   = median(meta_df$nCount_RNA),
    nFeature_mean   = mean(meta_df$nFeature_RNA),
    nFeature_median = median(meta_df$nFeature_RNA),
    mito_mean       = mean(meta_df$percent.mt),
    mito_median     = median(meta_df$percent.mt)
  )
}

# ---- before/after QC plots --------------------------------------------------
save_qc_plots <- function(pre_obj, post_obj, sname, outdir) {
  feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

  png(file.path(outdir, sprintf("%s_qc_before_vln.png", sname)),
      width = 1400, height = 500, res = 100, type = "cairo")
  print(VlnPlot(pre_obj, features = feats, layer = "counts", ncol = 3))
  dev.off()

  png(file.path(outdir, sprintf("%s_qc_before_scatter.png", sname)),
      width = 1200, height = 400, res = 125, type = "cairo")
  p1 <- FeatureScatter(pre_obj, "nCount_RNA", "percent.mt")   & NoLegend()
  p2 <- FeatureScatter(pre_obj, "nCount_RNA", "nFeature_RNA") & NoLegend()
  print(p1 | p2)
  dev.off()

  png(file.path(outdir, sprintf("%s_qc_after_vln.png", sname)),
      width = 1400, height = 500, res = 100, type = "cairo")
  print(VlnPlot(post_obj, features = feats, layer = "counts", ncol = 3))
  dev.off()

  png(file.path(outdir, sprintf("%s_qc_after_scatter.png", sname)),
      width = 1200, height = 400, res = 125, type = "cairo")
  p1 <- FeatureScatter(post_obj, "nCount_RNA", "percent.mt")   & NoLegend()
  p2 <- FeatureScatter(post_obj, "nCount_RNA", "nFeature_RNA") & NoLegend()
  print(p1 | p2)
  dev.off()

  png(file.path(outdir, sprintf("%s_nCount_hist.png", sname)),
      width = 900, height = 400, res = 125, type = "cairo")
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  hist(log10(pre_obj$nCount_RNA + 1), breaks = 60,
       main = paste0(sname, " - nCount (log10, pre)"),
       xlab = "log10(nCount_RNA + 1)", col = "steelblue")
  abline(v = log10(c(config$min_count, config$max_count) + 1), col = "red", lty = 2)
  hist(log10(post_obj$nCount_RNA + 1), breaks = 60,
       main = paste0(sname, " - nCount (log10, post)"),
       xlab = "log10(nCount_RNA + 1)", col = "steelblue")
  dev.off()

  invisible(NULL)
}

# ---- disk-backed merge helpers ----------------------------------------------
merge_many_disk <- function(files, sids, project = "snrna") {
  x <- RenameCells(readRDS(files[1]), add.cell.id = sids[1])
  if (length(files) > 1) for (k in 2:length(files)) {
    y <- RenameCells(readRDS(files[k]), add.cell.id = sids[k])
    x <- merge(x, y, merge.data = FALSE, project = project); rm(y); gc()
  }
  x
}

prep_for_merge <- function(x) {
  DefaultAssay(x) <- "RNA"
  if (inherits(x[["RNA"]], "Assay5")) x <- JoinLayers(x, assay = "RNA")
  DietSeurat(x, assays = "RNA", counts = TRUE, data = FALSE,
             scale.data = FALSE, dimreducs = NULL, graphs = NULL)
}

# =============================================================================
# Step 1: per-sample QC
# =============================================================================
write_log("=== Step 1: per-sample QC ===")
qc_summary_list <- list()

for (i in seq_along(sample_names)) {
  sname <- sample_names[i]
  write_log(sprintf("[%d/%d] %s", i, length(sample_names), sname))

  mtx <- Read10X(data.dir = data_dirs[i])
  obj <- CreateSeuratObject(counts = mtx, project = sname, min.cells = 3)
  obj$sample          <- sname
  obj$batch           <- sample_to_batch[[sname]]
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^(MT-|mt-)", assay = "RNA")
  obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]",    assay = "RNA")

  obj_pre <- obj

  # ---- QC filter (nCount + nFeature bounds + mito) --------------------------
  obj_post <- subset(obj_pre,
    subset = nCount_RNA   >  config$min_count    &
             nCount_RNA   <  config$max_count    &
             nFeature_RNA >  config$min_features &
             nFeature_RNA <  config$max_features &
             percent.mt   <  config$mito_threshold)

  write_log(sprintf("  pre=%d -> post=%d (removed=%d, %.1f%%)",
                    ncol(obj_pre), ncol(obj_post),
                    ncol(obj_pre) - ncol(obj_post),
                    100 * (ncol(obj_pre) - ncol(obj_post)) / ncol(obj_pre)))

  saveRDS(obj_pre,  file.path(qc_dir, sprintf("01_%s_PRE.rds",  sname)))
  saveRDS(obj_post, file.path(qc_dir, sprintf("01_%s_QCed.rds", sname)))
  save_qc_plots(obj_pre, obj_post, sname, qc_plot_dir)

  pre_s   <- qc_summary(obj_pre@meta.data,  sname, "pre")
  post_s  <- qc_summary(obj_post@meta.data, sname, "post")
  r_rate  <- round(100 * (pre_s$cells - post_s$cells) / pre_s$cells, 2)
  removed <- cbind(
    data.frame(sample = sname, state = "removed",
               cells = pre_s$cells - post_s$cells,
               nCount_mean = NA, nCount_median = NA,
               nFeature_mean = NA, nFeature_median = NA,
               mito_mean = NA, mito_median = NA),
    removed_rate = r_rate)
  qc_summary_list[[sname]] <- rbind(
    cbind(pre_s,  removed_rate = 0),
    cbind(post_s, removed_rate = r_rate),
    removed[, c(names(pre_s), "removed_rate")])
  rm(obj, obj_pre, obj_post, mtx); gc()
}

qc_df <- do.call(rbind, qc_summary_list)
write.csv(qc_df, file.path(qc_dir, "QC_summary_all_samples.csv"), row.names = FALSE)
write_log("Step 1 done")

# =============================================================================
# Step 2: full merge -> 02_merged.rds
# =============================================================================
write_log("=== Step 2: full merge ===")

rds_files  <- list.files(qc_dir, pattern = "^01_.*_QCed\\.rds$", full.names = TRUE)
sample_ids <- sub("^01_(.*)_QCed\\.rds$", "\\1", basename(rds_files))
names(rds_files) <- sample_ids

# ---- detect low-quality samples directly from saved metadata ----------------
write_log("computing QC stats...")
sample_stats <- lapply(names(rds_files), function(s) {
  meta <- readRDS(rds_files[[s]])@meta.data
  data.frame(sample = s, n_cells = nrow(meta),
             nFeature_med = median(meta$nFeature_RNA),
             nCount_med   = median(meta$nCount_RNA),
             mt_mean      = round(mean(meta$percent.mt), 2))
}) %>% bind_rows() %>% arrange(nFeature_med)
write.csv(sample_stats, file.path(outdir, "sample_qc_stats.csv"), row.names = FALSE)

bad_samples_final <- sample_stats %>%
  filter(nFeature_med < config$min_median_nfeature) %>%
  pull(sample)
write_log(sprintf("excluding %d samples: %s",
                  length(bad_samples_final), paste(bad_samples_final, collapse = ", ")))

rds_files  <- rds_files[!names(rds_files) %in% bad_samples_final]
sample_ids <- names(rds_files)

files_by_batch     <- split(rds_files, unname(sample_to_batch[sample_ids]))
chunk_size         <- 20
batch_merged_paths <- c()

for (b in names(files_by_batch)) {
  files_b <- files_by_batch[[b]]; sids_b <- names(files_b)
  write_log(sprintf("batch=%s: %d samples", b, length(files_b)))

  if (length(files_b) > chunk_size) {
    idx <- split(seq_along(files_b), ceiling(seq_along(files_b) / chunk_size))
    chunk_paths <- c()
    for (k in seq_along(idx)) {
      ii    <- idx[[k]]
      obj_k <- merge_many_disk(files_b[ii], sids_b[ii],
                               project = paste0("snrna_", b, "_c", k))
      obj_k$batch <- b
      out_k <- file.path(outdir, sprintf("02_merged_%s_chunk%03d.rds", b, k))
      saveRDS(obj_k, out_k, compress = FALSE); chunk_paths[k] <- out_k
      rm(obj_k); gc()
    }
    chunk_objs <- lapply(chunk_paths, readRDS)
    obj_b <- Reduce(function(x, y) merge(x, y), chunk_objs)
    rm(chunk_objs); gc()
  } else {
    obj_b <- merge_many_disk(files_b, sids_b, project = paste0("snrna_", b))
  }
  obj_b$batch <- b
  out_b <- file.path(outdir, sprintf("02_merged_%s.rds", b))
  saveRDS(obj_b, out_b, compress = FALSE)
  batch_merged_paths[b] <- out_b
  rm(obj_b); gc()
}

# ---- merge batches (largest first), then re-layer on common genes -----------
sizes   <- file.info(unname(batch_merged_paths))$size
batches <- rev(names(batch_merged_paths)[order(sizes)])

combined <- prep_for_merge(readRDS(batch_merged_paths[[batches[1]]]))
for (j in 2:length(batches)) {
  write_log(sprintf("merging batch: %s", batches[j]))
  tmp      <- prep_for_merge(readRDS(batch_merged_paths[[batches[j]]]))
  combined <- merge(combined, tmp, merge.data = FALSE, project = "snrna_combined")
  rm(tmp); gc()
}

DefaultAssay(combined) <- "RNA"
layer_names <- Layers(combined, assay = "RNA")
meta        <- combined@meta.data

genes_common <- Reduce(intersect, lapply(layer_names, function(ly) {
  rownames(LayerData(combined, assay = "RNA", layer = ly))
}))
write_log(sprintf("common genes: %d", length(genes_common)))

tmp_layer_dir <- file.path(outdir, "tmp_layers")
dir.create(tmp_layer_dir, showWarnings = FALSE)

# Write each layer's counts (common genes) to disk to bound peak memory.
rds_paths <- c()
for (ly in layer_names) {
  write_log(sprintf("saving layer: %s", ly))
  m        <- LayerData(combined, assay = "RNA", layer = ly)[genes_common, ]
  out_path <- file.path(tmp_layer_dir, paste0(gsub("[^A-Za-z0-9]", "_", ly), ".rds"))
  saveRDS(m, out_path, compress = FALSE)
  rds_paths <- c(rds_paths, out_path)
  rm(m); gc()
}
rm(combined); gc()

obj_list <- list()
for (i in seq_along(rds_paths)) {
  m <- readRDS(rds_paths[i])
  obj_list[[i]] <- CreateSeuratObject(counts = m, meta.data = meta[colnames(m), , drop = FALSE])
  write_log(sprintf("layer object %s: %d cells", basename(rds_paths[i]), ncol(obj_list[[i]])))
  rm(m); gc()
}

obj <- merge(obj_list[[1]], y = obj_list[2:length(obj_list)],
             merge.data = FALSE, project = "snrna_combined")
rm(obj_list); gc()

write_log(sprintf("final merged object: %d cells, %d features", ncol(obj), nrow(obj)))
saveRDS(obj, file.path(outdir, "02_merged.rds"), compress = FALSE)
write_log("Step 2 done: 02_merged.rds")

# ---- clean up Step 2 intermediates (keep only 02_merged.rds) ----------------
unlink(unname(batch_merged_paths))
unlink(list.files(outdir, pattern = "^02_merged_.*_chunk[0-9]+\\.rds$", full.names = TRUE))
unlink(tmp_layer_dir, recursive = TRUE)
write_log("cleaned Step 2 intermediates (per-batch merges, chunks, tmp_layers)")

# =============================================================================
# Step 3: normalization and PCA  (continues with the in-memory merged object)
# =============================================================================
write_log("=== Step 3: normalization & PCA ===")
full_dir <- file.path(outdir, "full")
dir.create(full_dir, showWarnings = FALSE, recursive = TRUE)

# With Seurat v5 layers, NormalizeData/FindVariableFeatures operate per layer,
# so HVGs are selected per sample as described in the methods.
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst",
                            nfeatures = config$n_variable_features, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, npcs = config$npcs, verbose = FALSE)

saveRDS(obj, file.path(outdir, "03_PCA_full.rds"))
write_log("Step 3 done: 03_PCA_full.rds")

# =============================================================================
# Step 4: Harmony integration, clustering, UMAP
# =============================================================================
write_log("=== Step 4: Harmony & UMAP ===")

dims_use <- choose_dims_elbow(
  obj, max_pcs = 50, min_pcs = 10,
  out_png = file.path(full_dir, "PCA_elbow_full.png"),
  title   = "PCA elbow - All samples"
)

# ---- pre-Harmony clustering & UMAP (for batch-effect comparison) ------------
obj <- FindNeighbors(obj, reduction = "pca", dims = dims_use,
                     graph.name = "pca_nn", verbose = FALSE)
obj <- FindClusters(obj, resolution = config$resolution, graph.name = "pca_nn",
                    random.seed = 777, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = dims_use,
               n.neighbors = config$n_neighbors, min.dist = config$min_dist,
               reduction.name = "umap_before", reduction.key = "UMAPbefore_",
               verbose = FALSE)

png(file.path(full_dir, "UMAP_preh_cluster_full.png"), width = 1200, height = 900, res = 150)
print(DimPlot(obj, reduction = "umap_before", label = TRUE, pt.size = 0.05) +
        ggtitle("Pre-Harmony - All"))
dev.off()

png(file.path(full_dir, "UMAP_preh_batch_full.png"), width = 1200, height = 900, res = 150)
print(DimPlot(obj, reduction = "umap_before", group.by = "batch", pt.size = 0.05) +
        ggtitle("Pre-Harmony batch - All"))
dev.off()

# ---- Harmony correction on sequencing batch ---------------------------------
obj_h <- harmony::RunHarmony(
  obj,
  group.by.vars    = "batch",
  reduction        = "pca",
  dims.use         = dims_use,
  reduction.save   = "harmony",
  plot_convergence = FALSE
)
rm(obj); gc()

# ---- post-Harmony clustering & UMAP -----------------------------------------
obj_h <- FindNeighbors(obj_h, reduction = "harmony", dims = dims_use,
                       graph.name = "harmony_nn", verbose = FALSE)
obj_h <- FindClusters(obj_h, resolution = config$resolution, graph.name = "harmony_nn",
                      random.seed = 777, verbose = FALSE)
obj_h <- RunUMAP(obj_h, reduction = "harmony", dims = dims_use,
                 n.neighbors = config$n_neighbors, min.dist = config$min_dist, verbose = FALSE)

# Post-clustering inspection showed that clusters above index 18 each contained
# very few cells (too small to annotate reliably), so they are dropped here.
# The retained 19 clusters (0-18) are carried into annotation.
max_cluster_id <- 18
keep_idents <- as.character(0:max_cluster_id)
dropped     <- setdiff(levels(droplevels(Idents(obj_h))), keep_idents)
if (length(dropped)) {
  drop_sizes <- table(Idents(obj_h))[dropped]
  write_log(sprintf("dropping %d small cluster(s): %s",
                    length(dropped),
                    paste(sprintf("%s(n=%d)", names(drop_sizes), drop_sizes), collapse = ", ")))
}
cells_keep <- WhichCells(obj_h, idents = keep_idents)
obj_h <- subset(obj_h, cells = cells_keep)
write_log(sprintf("after cluster filtering: cells=%d, clusters=%d",
                  ncol(obj_h), nlevels(droplevels(Idents(obj_h)))))

# ---- color palettes ---------------------------------------------------------
n_clusters            <- length(unique(Idents(obj_h)))
cluster_colors        <- make_colors(n_clusters)
names(cluster_colors) <- levels(Idents(obj_h))
batch_colors          <- c("tomato", "seagreen", "steelblue")
names(batch_colors)   <- unique(obj_h@meta.data$batch)

# ---- parse sample name: {genotype}_{condition}_{replicate} ------------------
obj_h@meta.data <- obj_h@meta.data %>%
  mutate(
    genotype  = sub("^([^_]+)_.*",        "\\1", sample),
    condition = sub("^[^_]+_([^_]+)_.*",  "\\1", sample),
    replicate = sub("^.*_([^_]+)$",       "\\1", sample),
    group     = paste(genotype, condition, sep = "_")
  )

# ---- UMAP: clusters ---------------------------------------------------------
png(file.path(full_dir, "UMAP_harmony_cluster_scattermore.png"),
    width = 1500, height = 1300, res = 300)
print(umap_scattermore(obj_h, reduction = "umap",
                       colors = cluster_colors, title = "Harmony - Cluster"))
dev.off()

# ---- UMAP: per-batch highlight, before vs after Harmony ---------------------
batches <- unique(obj_h@meta.data$batch)
plots_before <- mapply(plot_batch_highlight, batch_id = batches, color = batch_colors,
                       MoreArgs = list(obj = obj_h, reduction = "umap_before"), SIMPLIFY = FALSE)
plots_after  <- mapply(plot_batch_highlight, batch_id = batches, color = batch_colors,
                       MoreArgs = list(obj = obj_h, reduction = "umap"), SIMPLIFY = FALSE)
png(file.path(full_dir, "UMAP_batch_before_after.png"),
    width = 400 * length(batches), height = 900, res = 150)
row_before <- wrap_plots(plots_before, nrow = 1) +
  plot_annotation(title = "Before Harmony",
                  theme = theme(plot.title = element_text(size = 12, face = "bold")))
row_after  <- wrap_plots(plots_after, nrow = 1) +
  plot_annotation(title = "After Harmony",
                  theme = theme(plot.title = element_text(size = 12, face = "bold")))
print(row_before / row_after)
dev.off()

# ---- UMAP: group / condition / genotype -------------------------------------
p_group <- DimPlot(obj_h, reduction = "umap", group.by = "group",
                   pt.size = 0.05, raster = FALSE) +
  ggtitle("Harmony - Genotype x Condition") +
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
p_cond  <- DimPlot(obj_h, reduction = "umap", group.by = "condition",
                   pt.size = 0.05, raster = FALSE) + ggtitle("Harmony - Condition")
p_geno  <- DimPlot(obj_h, reduction = "umap", group.by = "genotype",
                   pt.size = 0.05, raster = FALSE) + ggtitle("Harmony - Genotype")
png(file.path(full_dir, "UMAP_harmony_group_full.png"), width = 1800, height = 1400, res = 150)
print((p_group / (p_cond | p_geno)) + plot_layout(heights = c(2, 1)))
dev.off()

saveRDS(obj_h, file.path(outdir, "04_Harmony_full_clean.rds"))
write_log("Step 4 done: 04_Harmony_full_clean.rds")

# =============================================================================
# Reproducibility
# =============================================================================
writeLines(capture.output(sessionInfo()),
           file.path(outdir, "sessionInfo_07_qc_integration.txt"))
write_log("=== 07 snRNA-seq pipeline complete ===")
cat("[OK] 07 snRNA-seq QC -> clustering done ->", outdir, "\n")
