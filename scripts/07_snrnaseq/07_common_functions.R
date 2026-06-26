################################################################################
# 07_common_functions.R
# Shared logging and plotting utilities for the snRNA-seq pipeline
# (07_qc_integration.R, 07_annotation.R, 07_hdWGCNA.R).
# Author : Yukyung Jun
# Updated: 2026-05  (extracted shared helpers, English comments, env-var paths)
#
# This file is SOURCED, not run directly. Each pipeline script does:
#   source(file.path(Sys.getenv("WORKDIR", unset = "."), "scripts",
#                    "07_common_functions.R"))
#
# Provides:
#   make_write_log()        - timestamped logger factory
#   make_colors()           - qualitative palette of arbitrary length
#   make_celltype_colors()  - stable cell-type palette (Unassigned = grey)
#   choose_dims_elbow()     - automatic PC selection (elbow + 85% variance)
#   read_gmt()              - minimal GMT reader
#   dimplot_scattermore()   - DimPlot replacement backed by scattermore
#   feature_scattermore()   - FeaturePlot replacement (gene expression)
#   score_scattermore()     - FeaturePlot replacement (prediction score)
#   umap_scattermore()      - UMAP with per-group centroid labels
#   plot_batch_highlight()  - highlight a single batch on a reduction
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scattermore)
  library(scales)
  library(RColorBrewer)
})

# =============================================================================
# Logging
# =============================================================================

# Return a logger that appends timestamped messages to `path` and echoes them
# to stderr. Usage: write_log <- make_write_log(file.path(outdir, "log.txt"))
make_write_log <- function(path) {
  force(path)
  function(msg) {
    cat(sprintf("[%s] %s\n", Sys.time(), msg), file = path, append = TRUE)
    message(msg)
  }
}

# =============================================================================
# Palettes
# =============================================================================

# Qualitative palette recycled to length n (used for clusters).
make_colors <- function(n) {
  base_cols <- c(brewer.pal(8, "Set1"),  brewer.pal(8, "Set2"),
                 brewer.pal(8, "Set3"),  brewer.pal(8, "Dark2"),
                 brewer.pal(8, "Paired"))
  rep_len(unique(base_cols), n)
}

# Stable cell-type palette; "Unassigned" is always grey and listed last.
make_celltype_colors <- function(celltypes) {
  celltypes      <- sort(unique(as.character(celltypes)))
  non_unassigned <- setdiff(celltypes, "Unassigned")
  n              <- length(non_unassigned)
  if (n == 0) return(c("Unassigned" = "#CCCCCC"))
  palettes   <- c("Paired", "Set1", "Set2", "Dark2", "Set3")
  pal_colors <- unique(unlist(lapply(palettes, function(p) {
    brewer.pal(brewer.pal.info[p, "maxcolors"], p)
  })))
  cols <- if (n <= length(pal_colors)) pal_colors[seq_len(n)]
          else colorRampPalette(pal_colors)(n)
  c(setNames(cols, non_unassigned), "Unassigned" = "#CCCCCC")
}

# =============================================================================
# Automatic PC selection
# =============================================================================

# Elbow (max distance from the first-to-last eigenvalue line) combined with the
# first PC reaching 85% cumulative variance; never fewer than `min_pcs`.
choose_dims_elbow <- function(obj, reduction = "pca", max_pcs = 50,
                              min_pcs = 10, out_png = NULL, title = "PCA elbow") {
  sdev <- obj@reductions[[reduction]]@stdev
  eig  <- sdev^2
  k    <- min(length(eig), max_pcs)
  x <- 1:k; y <- eig[1:k]
  x1 <- 1; y1 <- y[1]; x2 <- k; y2 <- y[k]
  denom <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
  d     <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) / max(denom, 1e-10)
  elbow <- max(min_pcs, min(which.max(d), k))
  cum85 <- which(cumsum(eig[1:k] / sum(eig[1:k])) >= 0.85)[1]
  dims_use <- 1:max(elbow, cum85, min_pcs)
  if (!is.null(out_png)) {
    png(out_png, width = 900, height = 600, res = 130)
    op <- par(mar = c(4.5, 4.5, 3, 1))
    plot(x, y, type = "b", pch = 16, main = title, xlab = "PC", ylab = "Eigenvalue")
    abline(a = y1 - (y2 - y1) / (x2 - x1) * x1, b = (y2 - y1) / (x2 - x1), lty = 2)
    points(elbow, y[elbow], pch = 21, bg = "red", cex = 1.5)
    text(elbow, y[elbow], paste0(" elbow=", elbow), pos = 3)
    par(op); dev.off()
  }
  message(sprintf("dims selected: 1:%d (elbow=%d, cum85=%d)",
                  max(dims_use), elbow, cum85))
  dims_use
}

# =============================================================================
# GMT reader
# =============================================================================

# Read a .gmt file into a list of {gene_set_name, description, genes}.
read_gmt <- function(file) {
  lapply(readLines(file), function(line) {
    parts <- strsplit(line, "\t")[[1]]
    list(gene_set_name = parts[1], description = parts[2], genes = parts[-c(1, 2)])
  })
}

# =============================================================================
# scattermore-based plotting (fast rasterized point layers for large objects)
# =============================================================================

# DimPlot replacement. Groups by `group.by` (or active idents); optional
# centroid labels. `colors` defaults to make_celltype_colors().
dimplot_scattermore <- function(obj, group.by = NULL, reduction = "umap",
                                colors = NULL, label = TRUE, pt.size = 1.2,
                                pixels = c(2000, 2000), title = "", alpha = 0.6) {
  emb <- as.data.frame(Embeddings(obj, reduction))
  colnames(emb)[1:2] <- c("DIM1", "DIM2")
  emb$group <- if (is.null(group.by)) as.character(Idents(obj))
               else as.character(obj@meta.data[[group.by]])
  emb$group <- factor(emb$group)
  if (is.null(colors)) colors <- make_celltype_colors(levels(emb$group))

  # Random shuffle so no single group is always drawn on top.
  emb <- emb[sample(nrow(emb)), ]

  p <- ggplot(emb, aes(x = DIM1, y = DIM2, color = group)) +
    geom_scattermore(pointsize = pt.size, pixels = pixels, alpha = alpha) +
    scale_color_manual(values = colors, drop = FALSE,
                       name = if (!is.null(group.by)) group.by else "Cluster") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
    ggtitle(title) +
    theme_classic(base_size = 11) +
    theme(plot.title      = element_text(face = "bold", size = 12),
          legend.text     = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"))

  if (label) {
    centroids <- emb %>%
      group_by(group) %>%
      summarise(x = median(DIM1), y = median(DIM2), .groups = "drop")
    p <- p + geom_label_repel(
      data = centroids, aes(x = x, y = y, label = group), inherit.aes = FALSE,
      size = 3, fontface = "bold", fill = scales::alpha("white", 0.7),
      label.size = 0.2, max.overlaps = 30)
  }
  p
}

# FeaturePlot replacement for a single gene; zero cells drawn grey underneath.
feature_scattermore <- function(obj, gene, reduction = "umap",
                                pixels = c(1500, 1500), pt.size = 0.8) {
  emb <- as.data.frame(Embeddings(obj, reduction))
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  emb$expr <- FetchData(obj, vars = gene)[[1]]
  emb_zero <- emb[emb$expr == 0, ]
  emb_expr <- emb[emb$expr >  0, ]
  emb_expr <- emb_expr[order(emb_expr$expr), ]

  ggplot() +
    geom_scattermore(data = emb_zero, aes(x = UMAP_1, y = UMAP_2),
                     color = "grey90", pointsize = pt.size, pixels = pixels) +
    geom_scattermore(data = emb_expr, aes(x = UMAP_1, y = UMAP_2, color = expr),
                     pointsize = pt.size * 1.5, pixels = pixels) +
    scale_color_gradientn(colors = c("#d0e8f5", "#2171b5", "#08306b"),
                          name = "Expression") +
    ggtitle(gene) +
    theme_classic(base_size = 10) +
    theme(plot.title        = element_text(face = "bold.italic", size = 11),
          legend.key.height = unit(0.5, "cm"))
}

# FeaturePlot replacement for a continuous prediction/QC score column.
score_scattermore <- function(obj, score_col, reduction = "umap",
                              pixels = c(1500, 1500), pt.size = 0.8, title = NULL) {
  emb <- as.data.frame(Embeddings(obj, reduction))
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  emb$score <- FetchData(obj, vars = score_col)[[1]]
  emb <- emb[order(emb$score), ]

  ggplot(emb, aes(x = UMAP_1, y = UMAP_2, color = score)) +
    geom_scattermore(pointsize = pt.size, pixels = pixels) +
    scale_color_viridis_c(name = "Score", option = "plasma") +
    ggtitle(if (is.null(title)) score_col else title) +
    theme_classic(base_size = 10) +
    theme(plot.title        = element_text(face = "bold", size = 11),
          legend.key.height = unit(0.5, "cm"))
}

# UMAP with median centroid text labels per group (used pre-annotation).
umap_scattermore <- function(obj, reduction = "umap", group.by = NULL,
                             colors = NULL, title = "", pt.size = 1.5,
                             pixels = c(3000, 2500)) {
  emb <- as.data.frame(Embeddings(obj, reduction = reduction))
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  emb$group <- if (is.null(group.by)) as.character(Idents(obj))
               else as.character(obj@meta.data[[group.by]])
  emb$group <- factor(emb$group)

  if (is.null(colors)) {
    colors <- make_colors(nlevels(emb$group))
    names(colors) <- levels(emb$group)
  }

  centroids <- emb %>%
    group_by(group) %>%
    summarise(x = median(UMAP_1), y = median(UMAP_2), .groups = "drop")

  ggplot(emb, aes(x = UMAP_1, y = UMAP_2, color = group)) +
    geom_scattermore(pointsize = pt.size, pixels = pixels, alpha = 0.6) +
    geom_text(data = centroids, aes(x = x, y = y, label = group),
              color = "black", size = 3, fontface = "bold", inherit.aes = FALSE) +
    scale_color_manual(values = colors) +
    ggtitle(title) +
    theme_classic(base_size = 12) +
    theme(legend.text     = element_text(size = 7),
          legend.key.size = unit(0.35, "cm"),
          plot.title      = element_text(face = "bold"))
}

# Highlight a single batch (grey background, colored foreground) on a reduction.
plot_batch_highlight <- function(obj, reduction, batch_id,
                                 color = "tomato", pt.size = 0.05) {
  cells_hl <- WhichCells(obj, expression = batch == batch_id)
  DimPlot(obj, reduction = reduction, cells.highlight = cells_hl,
          cols.highlight = color, cols = "lightgrey",
          pt.size = pt.size, sizes.highlight = pt.size, raster = FALSE) +
    ggtitle(batch_id) +
    theme(legend.position = "none",
          plot.title = element_text(size = 9, face = "bold"))
}
