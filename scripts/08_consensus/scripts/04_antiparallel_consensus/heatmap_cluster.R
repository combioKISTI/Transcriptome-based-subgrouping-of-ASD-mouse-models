#!/usr/bin/env Rscript
################################################################################
# 04_antiparallel_consensus -- 2-way clustering heatmap of antiparallel genes
#
# Input : matrix produced by `filter_input_by_overlapped.py`
#         (row 1 group / row 2 sample / row 3+ gene + log2FC values).
# Output: a 2-way clustered heatmap PDF with G1 and G2 column blocks
#         separated by a thin gap.
#
# Usage:
#   Rscript heatmap_cluster.R <input.tsv> <output.pdf>
################################################################################
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript heatmap_cluster.R <input.tsv> <output.pdf>")
}
input_file  <- args[1]
output_file <- args[2]

dt <- read.delim(input_file, header = FALSE, sep = "\t",
                 stringsAsFactors = FALSE, quote = "")
if (nrow(dt) < 3 || ncol(dt) < 3) {
  stop("Invalid input format: expected 3 rows x 3 columns or more")
}

groups_raw  <- as.character(dt[1, -1])
samples_raw <- as.character(dt[2, -1])
genes       <- as.character(dt[-c(1, 2), 1])
mat_raw     <- as.matrix(dt[-c(1, 2), -1, drop = FALSE])

mat <- suppressWarnings(apply(mat_raw, 2, as.numeric))
rownames(mat) <- make.unique(genes)
colnames(mat) <- samples_raw

groups <- factor(groups_raw, levels = unique(groups_raw))
names(groups) <- samples_raw
group_levels <- levels(groups)

# Cluster columns within each group (preserve the inter-group split)
cluster_within_group <- function(m, g) {
  out <- character(0)
  for (grp in levels(g)) {
    cols <- names(g)[g == grp]
    sub_m <- m[, cols, drop = FALSE]
    hc <- hclust(dist(t(sub_m)), method = "complete")
    out <- c(out, cols[hc$order])
  }
  out
}
ord_cols <- cluster_within_group(mat, groups)
mat <- mat[, ord_cols]
groups <- groups[ord_cols]

col_fun <- colorRamp2(
  c(-2, -0.5, -0.2, 0, 0.2, 0.5, 2),
  c("#053061", "#2166ac", "#92c5de", "#FFFFBF", "#f4a582", "#d6604d", "#67001f")
)
clamp <- function(m, lo = -2, hi = 2) { m[m < lo] <- lo; m[m > hi] <- hi; m }
mat_vis <- clamp(mat)

group_colors <- c(Group1 = "#8dba43", Group2 = "#e99436")
miss <- setdiff(group_levels, names(group_colors))
if (length(miss) > 0) group_colors <- c(group_colors, setNames(rep("#bbbbbb", length(miss)), miss))

mat_list <- lapply(group_levels, function(g) mat_vis[, groups == g, drop = FALSE])
names(mat_list) <- group_levels

heatmaps <- list()
for (i in seq_along(mat_list)) {
  sub_mat <- mat_list[[i]]
  sub_grp <- group_levels[i]

  top_annot <- HeatmapAnnotation(
    df  = data.frame(Group = factor(rep(sub_grp, ncol(sub_mat)), levels = group_levels)),
    col = list(Group = group_colors),
    annotation_name_gp = gpar(fontsize = 0),
    annotation_legend_param = list(Group = list(at = group_levels, labels = group_levels)),
    show_legend = (i == 1)
  )
  ht_i <- Heatmap(
    sub_mat, name = "Fold Change", col = col_fun,
    top_annotation = top_annot,
    cluster_rows = TRUE, cluster_columns = TRUE,
    show_row_dend = TRUE, show_column_dend = TRUE,
    show_row_names = TRUE, row_names_gp = gpar(fontsize = 6),
    show_column_names = TRUE,
    column_title = NULL, show_heatmap_legend = (i == 1),
    heatmap_legend_param = list(title = "Fold Change",
                                at = c(-2, -1, 0, 1, 2),
                                labels = c("-2", "-1", "0", "1", "2"))
  )
  heatmaps[[i]] <- ht_i
}
ht_all <- Reduce(`+`, heatmaps)

pdf(output_file, width = 10, height = 11, onefile = TRUE)
draw(ht_all, heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE, gap = unit(2, "mm"))
dev.off()
cat(sprintf("[OK] heatmap PDF: %s\n", output_file))
