#!/usr/bin/env Rscript
################################################################################
# 03_wilcoxon -- One-sample one-sided Wilcoxon signed-rank test (BH-adjusted)
#
# Manuscript Methods (verbatim):
#   For every gene, two one-sided one-sample Wilcoxon signed-rank tests were
#   performed against a null median log2FC of 0 -- one with alternative = "less"
#   (testing for downregulation) and one with alternative = "greater" (testing
#   for upregulation) -- using wilcox.test() in R v4.3.3.
#   P-values from the two directions were independently adjusted across genes
#   by the Benjamini-Hochberg procedure, yielding Padj_less and Padj_greater;
#   a gene was called WilcoxonDown if Padj_less < 0.05 and WilcoxonUp if
#   Padj_greater < 0.05.
#
# This script processes one group's input matrix (e.g. G1) and emits:
#   (a) a heatmap PDF (full matrix + a Down/Up-only 2-way clustered matrix), and
#   (b) a labelled, sorted TSV with per-gene Wilcoxon stats and call.
#
# Usage:
#   Rscript wilcoxon_one_sided_BH.R <input.tsv> <output_prefix.pdf> [alpha=0.05] [epsilon=0.0]
#
# Input format: row 1 group label, row 2 sample ID, row 3+ gene name plus per-
# sample log2 fold-change values.
################################################################################
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript wilcoxon_one_sided_BH.R <input.tsv> <output.pdf> [alpha=0.05] [epsilon=0.0]")
}
input_file  <- args[1]
output_file <- args[2]

alpha   <- 0.05
epsilon <- 0.0
if (length(args) > 2) {
  for (kv in args[3:length(args)]) {
    if (!grepl("=", kv)) next
    key <- sub("=.*$", "", kv)
    val <- sub("^.*?=", "", kv)
    if (key == "alpha")   suppressWarnings(alpha   <- as.numeric(val))
    if (key == "epsilon") suppressWarnings(epsilon <- as.numeric(val))
  }
}
stopifnot(alpha > 0, alpha < 1, epsilon >= 0)

# -- Load --------------------------------------------------------------------
dt <- read.delim(input_file, header = FALSE, sep = "\t",
                 check.names = FALSE, stringsAsFactors = FALSE, quote = "")
if (nrow(dt) < 3 || ncol(dt) < 3) {
  stop("Invalid input format: expected 2 header rows + data rows, and >=3 columns")
}

header1_first <- as.character(dt[1, 1])
header2_first <- as.character(dt[2, 1])
groups_raw    <- as.character(dt[1, -1])
samples_raw   <- as.character(dt[2, -1])
genes         <- as.character(dt[-c(1, 2), 1])
mat_raw       <- as.matrix(dt[-c(1, 2), -1, drop = FALSE])

mat_num <- suppressWarnings(apply(mat_raw, 2, as.numeric))
rownames(mat_num) <- make.unique(genes)
colnames(mat_num) <- samples_raw

groups <- factor(groups_raw, levels = unique(groups_raw))
names(groups) <- samples_raw

# -- Per-sample direction labels (informational only) ------------------------
label_value <- function(x, eps = 0) {
  ifelse(x < -eps, "down", ifelse(x > eps, "up", "no_change"))
}
label_mat <- apply(mat_num, 2, label_value, eps = epsilon)

# -- Wilcoxon signed-rank tests (one-sided, both alternatives) ---------------
wtest <- function(x, alt) {
  x <- x[!is.na(x)]
  if (length(x) == 0 || all(x == 0)) return(1)
  out <- tryCatch(wilcox.test(x, mu = 0, alternative = alt)$p.value,
                  error = function(e) 1)
  ifelse(is.finite(out), out, 1)
}
p_less    <- apply(mat_num, 1, wtest, alt = "less")
p_greater <- apply(mat_num, 1, wtest, alt = "greater")

padj_less    <- p.adjust(p_less,    method = "BH")
padj_greater <- p.adjust(p_greater, method = "BH")

mean_fc <- rowMeans(mat_num, na.rm = TRUE)

wilcoxon_call <- ifelse(padj_less    < alpha & mean_fc < 0, "WilcoxonDown",
                 ifelse(padj_greater < alpha & mean_fc > 0, "WilcoxonUp",   "NoChange"))

# -- Sort: Down -> Up -> NoChange, within each group by mean_fc --------------
priority <- factor(wilcoxon_call,
                   levels = c("WilcoxonDown", "WilcoxonUp", "NoChange"),
                   ordered = TRUE)
ord <- order(priority, mean_fc)

mat_sorted          <- mat_num[ord, , drop = FALSE]
label_mat_sorted    <- label_mat[ord, , drop = FALSE]
mean_fc_sorted      <- mean_fc[ord]
call_sorted         <- wilcoxon_call[ord]
p_less_sorted       <- p_less[ord]
padj_less_sorted    <- padj_less[ord]
p_greater_sorted    <- p_greater[ord]
padj_greater_sorted <- padj_greater[ord]

# -- Visualisation -----------------------------------------------------------
col_fun <- colorRamp2(
  c(-2, -0.5, -0.2, 0, 0.2, 0.5, 2),
  c("#053061", "#2166ac", "#92c5de", "#FFFFBF", "#f4a582", "#d6604d", "#67001f")
)
clamp <- function(m, lo = -2, hi = 2) { m[m < lo] <- lo; m[m > hi] <- hi; m }
mat_vis_all <- clamp(mat_sorted)

call_colors  <- c(WilcoxonDown = "#377eb8", WilcoxonUp = "#e41a1c", NoChange = "#999999")
group_colors <- c(Group1 = "#8dba43", Group2 = "#e99436")
miss <- setdiff(levels(groups), names(group_colors))
if (length(miss) > 0) group_colors <- c(group_colors, setNames(rep("#bbbbbb", length(miss)), miss))

top_annot <- HeatmapAnnotation(
  df  = data.frame(Group = groups),
  col = list(Group = group_colors),
  annotation_name_side = "left",
  show_legend = TRUE
)
row_annot <- rowAnnotation(
  Call = factor(call_sorted, levels = names(call_colors)),
  col  = list(Call = call_colors),
  annotation_name_side = "top",
  width = unit(4, "mm")
)

ht_all <- Heatmap(
  mat_vis_all, name = "FC", col = col_fun,
  top_annotation = top_annot, left_annotation = row_annot,
  cluster_rows = FALSE, cluster_columns = TRUE,
  clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
  show_row_names = FALSE, show_column_names = TRUE,
  heatmap_legend_param = list(title = "Fold Change",
                              at = c(-2, -1, 0, 1, 2),
                              labels = c("-2", "-1", "0", "1", "2"))
)

tag        <- sprintf("alpha%.3f_eps%.2f", alpha, epsilon)
base_noext <- sub("\\.pdf$", "", output_file, ignore.case = TRUE)
pdf_all    <- sprintf("%s.%s.pdf", base_noext, tag)
tsv_out    <- sprintf("%s.%s_wilcoxon_labeled_sorted.tsv", base_noext, tag)

pdf(pdf_all, width = 9, height = 11, onefile = TRUE)
draw(ht_all, heatmap_legend_side = "right",
     annotation_legend_side = "right", merge_legend = TRUE)
dev.off()

# Down/Up-only 2-way clustering
keep_idx <- which(call_sorted != "NoChange")
if (length(keep_idx) >= 2) {
  mat_dnup     <- mat_sorted[keep_idx, , drop = FALSE]
  mat_vis_dnup <- clamp(mat_dnup)
  call_dnup    <- call_sorted[keep_idx]

  row_annot_dnup <- rowAnnotation(
    Call = factor(call_dnup, levels = names(call_colors)),
    col  = list(Call = call_colors),
    annotation_name_side = "top",
    width = unit(4, "mm")
  )
  ht_dnup <- Heatmap(
    mat_vis_dnup, name = "FC", col = col_fun,
    top_annotation = top_annot, left_annotation = row_annot_dnup,
    cluster_rows = TRUE, cluster_columns = TRUE,
    clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
    clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
    show_row_names = FALSE, show_column_names = TRUE,
    heatmap_legend_param = list(title = "Fold Change",
                                at = c(-2, -1, 0, 1, 2),
                                labels = c("-2", "-1", "0", "1", "2"))
  )
  pdf_dnup <- sprintf("%s.%s.onlyDownUp_clustered.pdf", base_noext, tag)
  pdf(pdf_dnup, width = 9, height = 11, onefile = TRUE)
  draw(ht_dnup, heatmap_legend_side = "right",
       annotation_legend_side = "right", merge_legend = TRUE)
  dev.off()
} else {
  warning("Fewer than 2 Down/Up calls -- skipping the 2-way clustering heatmap.")
}

# -- TSV output --------------------------------------------------------------
samples_final <- colnames(mat_sorted)
groups_final  <- as.character(groups[samples_final])

con <- file(tsv_out, open = "wt", encoding = "UTF-8")
writeLines(paste(c(header1_first, groups_final),  collapse = "\t"), con = con)
writeLines(paste(c(header2_first, samples_final), collapse = "\t"), con = con)
write.table(cbind(rownames(mat_sorted), mat_sorted),
            file = con, sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)
writeLines("", con = con)

label_cols <- as.data.frame(label_mat_sorted, stringsAsFactors = FALSE)
colnames(label_cols) <- paste0(colnames(label_cols), "_label")
df_out <- data.frame(
  Gene         = rownames(mat_sorted),
  MeanFC       = mean_fc_sorted,
  P_less       = p_less_sorted,
  Padj_less    = padj_less_sorted,
  P_greater    = p_greater_sorted,
  Padj_greater = padj_greater_sorted,
  WilcoxonCall = call_sorted,
  label_cols,
  check.names = FALSE
)
write.table(df_out, file = con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

writeLines(capture.output(sessionInfo()),
           sprintf("%s.sessionInfo.txt", base_noext))
message(sprintf("[OK] heatmap     : %s", pdf_all))
message(sprintf("[OK] sorted TSV  : %s", tsv_out))
