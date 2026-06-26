################################################################################
# 07. snRNA-seq cell-type-specific hdWGCNA
# Author : Yukyung Jun
# Updated: 2026-05  (env-var paths, English comments, standalone run,
#                    minModuleSize unified to 50, removed unused helpers,
#                    per-cell-type RDS save made explicit via config flag)
#
# Pipeline (per cell type):
#   network cells = WT vehicle (B6 Wv) + mutant vehicle (G1/G2 Hv)
#   mouse -> human ortholog -> metacells (k=25) -> soft power (R2>=0.80,
#   slope<0, power>=5) -> signed hybrid network (minModuleSize 50,
#   mergeCutHeight 0.20) -> module eigengenes / connectivity -> hub genes
#   -> module-trait correlation -> GMT gene-set overlap (Fisher)
#
# Tools   : Seurat v5, hdWGCNA, WGCNA, ComplexHeatmap, enrichR (offline)
#           Mouse->human orthologs via a precomputed biomaRt cache
#
# Required environment variables (defaults match the original project tree):
#   WORKDIR      - project root
#   R_LIBS_USER  - optional; prepended to .libPaths() if set
# Inputs:
#   $WORKDIR/seurat/05_annotation_full.rds          (from 07_annotation.R)
#   $WORKDIR/seurat/mouse2human_orthologs.rds        (biomaRt cache)
#   $WORKDIR/ref/genesets/*.gmt
# Output:
#   $WORKDIR/seurat/full/annotation/hdWGCNA/<celltype>/...
################################################################################

# Optional: prepend a user library path (kept out of source for portability).
r_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(r_lib)) .libPaths(r_lib)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(WGCNA)
  library(circlize)
  library(tidyr)
  library(methods)
})

# enrichR attempts a network call on attach; on an isolated HPC node this hangs.
# Patch its .onAttach to a no-op so the namespace can load offline.
unlockBinding(".onAttach", asNamespace("enrichR"))
assignInNamespace(".onAttach", function(libname, pkgname) {
  message("Welcome to enrichR (offline mode)")
}, ns = "enrichR")

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(hdWGCNA)
})

filter <- dplyr::filter   # guard against masking by WGCNA / ComplexHeatmap
enableWGCNAThreads(nThreads = 8)

# =============================================================================
# Config (paths from environment; analysis parameters inline)
# =============================================================================
workdir   <- Sys.getenv("WORKDIR", unset = "/blues/scratch/yukyung/EJKim/Drug_Response")
source(file.path(workdir, "scripts", "07_common_functions.R"))

outdir    <- file.path(workdir, "seurat")
full_dir  <- file.path(outdir, "full")
annot_dir <- file.path(full_dir, "annotation")
dir_wgcna <- file.path(annot_dir, "hdWGCNA")
dir.create(dir_wgcna, recursive = TRUE, showWarnings = FALSE)

config <- list(
  metacell_k        = 25,
  min_module_size   = 50,     # uniform across all cell types
  merge_cut_height  = 0.20,
  network_type      = "signed hybrid"
)

write_log <- make_write_log(file.path(outdir, "analysis_log.txt"))
write_log("07 hdWGCNA started")

# =============================================================================
# 1. Load annotated object
# =============================================================================
obj_h <- readRDS(file.path(outdir, "05_annotation_full.rds"))
obj_h@meta.data$cond2 <- sub("[MF]$", "", obj_h@meta.data$condition)  # strip sex suffix
meta <- obj_h@meta.data

CT_LEVELS <- sort(unique(meta$celltype))
write_log(sprintf("cell types (%d): %s", length(CT_LEVELS), paste(CT_LEVELS, collapse = ", ")))

# =============================================================================
# 2. Mouse -> human ortholog conversion (from precomputed biomaRt cache)
# =============================================================================
human_rds <- file.path(outdir, "05_annotation_full_human.rds")

if (file.exists(human_rds)) {
  write_log("human Seurat object loaded from cache")
  seurat_human <- readRDS(human_rds)
} else {
  ortho_cache <- file.path(outdir, "mouse2human_orthologs.rds")
  if (!file.exists(ortho_cache)) stop("ortholog cache not found: ", ortho_cache)

  write_log("building human Seurat object from ortholog table")
  ortho_table  <- readRDS(ortho_cache)
  mouse_genes  <- rownames(obj_h)
  human_mapped <- ortho_table[ortho_table$mouse_symbol %in% mouse_genes, ]
  human_mapped <- human_mapped[!duplicated(human_mapped$mouse_symbol), ]

  seurat_human <- obj_h[human_mapped$mouse_symbol, ]
  rownames(seurat_human) <- human_mapped$human_symbol[match(rownames(seurat_human),
                                                            human_mapped$mouse_symbol)]
  write_log(sprintf("mouse genes: %d -> human mapped: %d", length(mouse_genes), nrow(human_mapped)))
  saveRDS(seurat_human, human_rds)
}

# =============================================================================
# 3. Network cells: WT vehicle (B6 Wv) + mutant vehicle (G1/G2 Hv)
# =============================================================================
get_cells <- function(grp_val = NULL, cond_val = NULL) {
  idx <- rep(TRUE, nrow(meta))
  if (!is.null(grp_val))  idx <- idx & meta$group %in% grp_val
  if (!is.null(cond_val)) idx <- idx & meta$cond2 %in% cond_val
  rownames(meta)[idx]
}
cells_network <- c(get_cells("B6", "Wv"), get_cells(c("G1", "G2"), "Hv"))
write_log(sprintf("network cells (total): %d", length(cells_network)))

# =============================================================================
# 4. Gene-set files (for module overlap; read_gmt from 07_common_functions.R)
# =============================================================================
gmt_files <- c("asd.v1.symbols.gmt", "asdrisk.jung.v1.gmt",
               "Gandal 2022 WGCNA.gmt", "Gupta 2014 WGCNA.gmt",
               "Cortex p25 WGCNA.gmt", "Cortex p40 WGCNA.gmt")

# =============================================================================
# 5. hdWGCNA per cell type
# =============================================================================
run_hdWGCNA <- function(seurat_full, celltype, cells_network, out_base, soft_power = NULL) {

  ct_safe <- gsub("[/ ]", "_", celltype)
  out_dir <- file.path(out_base, ct_safe)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_log(sprintf("===== hdWGCNA: %s =====", celltype))

  # ---- subset to this cell type ---------------------------------------------
  cells_ct <- intersect(cells_network, rownames(meta)[meta$celltype == celltype])
  write_log(sprintf("  network cells: %d", length(cells_ct)))
  if (length(cells_ct) < 200) {
    write_log(sprintf("  SKIP: too few cells (%d)", length(cells_ct)))
    return(NULL)
  }
  seurat_sub <- subset(seurat_full, cells = cells_ct)
  seurat_sub <- JoinLayers(seurat_sub)

  # ---- metacells (grouped by cell type + sample identity) -------------------
  seurat_sub <- SetupForWGCNA(seurat_sub, gene_select = "fraction",
                              fraction = 0.05, wgcna_name = celltype)
  seurat_sub <- MetacellsByGroups(
    seurat_sub,
    group.by         = c("celltype", "unit_id"),
    k                = config$metacell_k,
    target_metacells = 500,
    max_shared       = 10,
    ident.group      = "celltype",
    wgcna_name       = celltype
  )
  seurat_sub <- NormalizeMetacells(seurat_sub, wgcna_name = celltype)
  write_log(sprintf("  metacells: %d", ncol(GetMetacellObject(seurat_sub, celltype))))

  # ---- soft-thresholding power ----------------------------------------------
  seurat_sub <- SetDatExpr(seurat_sub, group_name = celltype, group.by = "celltype",
                           assay = "RNA", layer = "data", wgcna_name = celltype)
  seurat_sub <- TestSoftPowers(seurat_sub, networkType = config$network_type,
                               wgcna_name = celltype)
  sp_table <- GetPowerTable(seurat_sub, wgcna_name = celltype)

  if (is.null(soft_power)) {
    # primary criterion: scale-free fit R2 >= 0.80, slope < 0, power >= 5
    sp_ok <- sp_table %>% filter(SFT.R.sq >= 0.80, slope < 0, Power >= 5) %>% arrange(Power)
    if (nrow(sp_ok) > 0) {
      soft_power <- sp_ok$Power[1]
    } else {
      sp_ok <- sp_table %>% filter(SFT.R.sq >= 0.80, Power >= 5) %>% arrange(Power)
      if (nrow(sp_ok) > 0) {
        soft_power <- sp_ok$Power[1]
      } else {
        best <- sp_table %>% filter(Power >= 5) %>% arrange(desc(SFT.R.sq)) %>% slice_head(n = 1)
        soft_power <- if (nrow(best) > 0) best$Power[1] else 12
        write_log(sprintf("  WARNING: no power with R2 >= 0.8; using power=%d", soft_power))
      }
    }
  }
  write_log(sprintf("  soft power: %d", soft_power))

  p_sp <- PlotSoftPowers(seurat_sub, wgcna_name = celltype)
  p_sp <- lapply(p_sp, function(p) {
    p + geom_vline(xintercept = soft_power, color = "red", linetype = "dashed", linewidth = 0.8)
  })
  pdf(file.path(out_dir, "pickSoftThreshold.pdf"), width = 10, height = 6)
  print(wrap_plots(p_sp, ncol = 2)); dev.off()

  # ---- network construction -------------------------------------------------
  seurat_sub <- ConstructNetwork(
    seurat_sub,
    soft_power     = soft_power,
    networkType    = config$network_type, TOMType = "signed",
    minModuleSize  = config$min_module_size,
    mergeCutHeight = config$merge_cut_height, detectCutHeight = 0.995,
    wgcna_name = celltype, overwrite_tom = TRUE
  )
  pdf(file.path(out_dir, "clusterDendrogram.pdf"), width = 8, height = 5)
  PlotDendrogram(seurat_sub, wgcna_name = celltype, main = celltype); dev.off()

  # ---- module eigengenes & connectivity -------------------------------------
  seurat_sub <- ScaleData(seurat_sub, features = VariableFeatures(seurat_sub))
  seurat_sub <- ModuleEigengenes(seurat_sub, group.by.vars = "group", wgcna_name = celltype)
  seurat_sub <- ModuleConnectivity(seurat_sub, group.by = "celltype",
                                   group_name = celltype, wgcna_name = celltype)

  MEs   <- GetMEs(seurat_sub, harmonized = FALSE, wgcna_name = celltype)
  ME_df <- as.data.frame(MEs)
  write.table(cbind(Sample = rownames(ME_df), ME_df),
              file.path(out_dir, "moduleEigenGenes.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)

  # ---- module membership ----------------------------------------------------
  modules     <- GetModules(seurat_sub, wgcna_name = celltype)
  gene_col    <- intersect(c("gene_name", "gene"), colnames(modules))[1]
  module_list <- setdiff(unique(modules$module), "grey")
  module_mat  <- data.frame(Gene = modules[[gene_col]], Module = modules$module,
                            stringsAsFactors = FALSE)

  # number modules by descending gene count; grey -> ME0
  gene_counts    <- sort(table(module_mat$Module[module_mat$Module != "grey"]), decreasing = TRUE)
  module_colors  <- names(gene_counts)
  module_num_map <- setNames(paste0(ct_safe, "-ME", seq_along(module_colors)), module_colors)
  module_num_map[["grey"]] <- paste0(ct_safe, "-ME0")
  module_mat$Module_num    <- module_num_map[as.character(module_mat$Module)]

  write.table(module_mat, file.path(out_dir, "moduleGenes.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  module_summary <- aggregate(Gene ~ Module + Module_num, data = module_mat, FUN = length)
  colnames(module_summary)[3] <- "GeneCount"
  module_summary <- module_summary[order(as.integer(sub(".*-ME", "", module_summary$Module_num))), ]
  write.table(module_summary, file.path(out_dir, "moduleGenesNum.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  write_log(sprintf("  modules: %d", length(module_list)))

  # ---- hub genes ------------------------------------------------------------
  hub_genes <- GetHubGenes(seurat_sub, n_hubs = 10, wgcna_name = celltype)
  hub_genes$module_num <- module_num_map[as.character(hub_genes$module)]
  write.table(hub_genes, file.path(out_dir, "hubGenes.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)

  # ---- module-trait correlation (Sex / Genotype / Group1 / Group2) ----------
  seurat_sub$Genotype_num <- ifelse(seurat_sub@meta.data$group %in% c("G1", "G2"), 1, 0)
  seurat_sub$Group1_num   <- ifelse(seurat_sub@meta.data$group == "G1", 1, 0)
  seurat_sub$Group2_num   <- ifelse(seurat_sub@meta.data$group == "G2", 1, 0)
  seurat_sub$Sex_num      <- ifelse(seurat_sub@meta.data$sex == "M", 1, 0)

  seurat_sub <- ModuleTraitCorrelation(
    seurat_sub,
    traits     = c("Genotype_num", "Group1_num", "Group2_num", "Sex_num"),
    group.by   = "celltype", wgcna_name = celltype
  )
  mt_cor  <- GetModuleTraitCorrelation(seurat_sub, wgcna_name = celltype)
  cor_df  <- as.data.frame(t(mt_cor$cor[[celltype]]))
  fdr_df  <- as.data.frame(t(mt_cor$fdr[[celltype]]))
  pval_df <- as.data.frame(t(mt_cor$pval[[celltype]]))
  colnames(cor_df)  <- gsub("_num$", "", colnames(cor_df))
  colnames(fdr_df)  <- gsub("_num$", "", colnames(fdr_df))
  colnames(pval_df) <- gsub("_num$", "", colnames(pval_df))

  col_order <- c("Sex", "Genotype", "Group1", "Group2")
  cor_df$module_num  <- module_num_map[as.character(rownames(cor_df))]
  cor_df$celltype    <- celltype
  cor_df  <- cor_df[, c(col_order, "module_num", "celltype")]
  pval_df$module_num <- module_num_map[as.character(rownames(pval_df))]
  pval_df <- pval_df[, c(col_order, "module_num")]
  fdr_df$module_num  <- module_num_map[as.character(rownames(fdr_df))]
  fdr_df  <- fdr_df[, c(col_order, "module_num")]

  write.table(cor_df,  file.path(out_dir, "moduleTraitCor.txt"),  quote = FALSE, sep = "\t")
  write.table(pval_df, file.path(out_dir, "moduleTraitPval.txt"), quote = FALSE, sep = "\t")
  write.table(fdr_df,  file.path(out_dir, "moduleTraitFDR.txt"),  quote = FALSE, sep = "\t")

  cor_wide <- cor_df[, col_order, drop = FALSE]; rownames(cor_wide) <- cor_df$module_num
  fdr_wide <- fdr_df[, col_order, drop = FALSE]; rownames(fdr_wide) <- fdr_df$module_num

  sig_mat <- matrix("", nrow(cor_wide), ncol(cor_wide), dimnames = dimnames(cor_wide))
  sig_mat[!is.na(fdr_wide) & fdr_wide < 0.05]  <- "*"
  sig_mat[!is.na(fdr_wide) & fdr_wide < 0.01]  <- "**"
  sig_mat[!is.na(fdr_wide) & fdr_wide < 0.001] <- "***"
  mods <- rownames(cor_wide); mods <- mods[order(as.integer(sub(".*-ME", "", mods)))]

  pdf(file.path(out_dir, "corModuleTraits.pdf"), width = 9, height = 10)
  ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(
      as.matrix(cor_wide[mods, ]), name = "Pearson r",
      col = circlize::colorRamp2(c(-0.4, 0, 0.4), c("#4575B4", "white", "#D73027")),
      cluster_columns = FALSE, cluster_rows = FALSE,
      row_names_side = "left", column_names_rot = 45, na_col = "grey90",
      cell_fun = function(j, i, x, y, w, h, fill) {
        s <- sig_mat[mods[i], colnames(cor_wide)[j]]
        if (!is.na(s) && s != "")
          grid::grid.text(s, x, y, gp = grid::gpar(fontsize = 10, fontface = "bold"))
      }),
    heatmap_legend_side = "right")
  dev.off()

  # ---- Fisher overlap with external gene sets (GMT) -------------------------
  all_genes <- module_mat$Gene
  or_final  <- NULL
  for (gmt_f in gmt_files) {
    gmt_path <- file.path(workdir, "ref/genesets", gmt_f)
    if (!file.exists(gmt_path)) { write_log(paste("  GMT not found:", gmt_f)); next }
    gmt_mat <- read_gmt(gmt_path)
    or_mat  <- NULL
    for (gs in gmt_mat) {
      or_vec <- sapply(module_list, function(mod) {
        mg <- module_mat$Gene[module_mat$Module == mod]
        a  <- length(intersect(mg, gs$genes))
        b  <- length(setdiff(mg, gs$genes))
        cc <- length(setdiff(gs$genes, mg))
        d  <- length(setdiff(all_genes, union(mg, gs$genes)))
        fisher.test(matrix(c(a, b, cc, d), 2, byrow = TRUE))$estimate
      })
      or_mat <- cbind(or_mat, matrix(or_vec, ncol = 1,
                                     dimnames = list(module_list, gs$gene_set_name)))
    }
    or_final <- cbind(or_final, or_mat)

    or_mat_num <- or_mat; rownames(or_mat_num) <- module_num_map[as.character(rownames(or_mat))]
    traits_use <- col_order[col_order %in% colnames(cor_wide)]
    common_mods <- intersect(rownames(cor_wide), rownames(or_mat_num))
    if (length(common_mods) == 0) next
    diff_sub <- cor_wide[common_mods, traits_use, drop = FALSE]
    fdr_sub  <- fdr_wide[common_mods, traits_use, drop = FALSE]
    sig2 <- matrix("", nrow(fdr_sub), ncol(fdr_sub), dimnames = dimnames(fdr_sub))
    sig2[!is.na(fdr_sub) & fdr_sub < 0.05]  <- "*"
    sig2[!is.na(fdr_sub) & fdr_sub < 0.01]  <- "**"
    sig2[!is.na(fdr_sub) & fdr_sub < 0.001] <- "***"

    pdf(file.path(out_dir, paste0("heatmap_", gsub(".gmt", "", gmt_f), ".pdf")),
        width = 14, height = 8)
    ComplexHeatmap::draw(
      ComplexHeatmap::Heatmap(
        diff_sub, name = "ME diff",
        col = circlize::colorRamp2(c(-0.5, 0, 0.5), c("#4575B4", "white", "#D73027")),
        cluster_columns = FALSE, row_names_side = "left", column_names_rot = 45,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if (!is.na(sig2[i, j]) && sig2[i, j] != "")
            grid::grid.text(sig2[i, j], x, y, gp = grid::gpar(fontsize = 10, fontface = "bold"))
        }) +
      ComplexHeatmap::Heatmap(
        or_mat_num[common_mods, , drop = FALSE], name = "Odds ratio",
        col = circlize::colorRamp2(c(0, 3, 10), c("white", "skyblue", "blue")),
        cluster_columns = FALSE, cluster_rows = FALSE),
      ht_gap = grid::unit(5, "mm"), heatmap_legend_side = "right", auto_adjust = FALSE)
    dev.off()
  }
  if (!is.null(or_final)) {
    rownames(or_final) <- module_num_map[as.character(rownames(or_final))]
    write.table(or_final, file.path(out_dir, "module_supp2.txt"), quote = FALSE, sep = "\t")
  }

  # NOTE: enrichR GO enrichment is run separately (07_hdWGCNA_enrichR_local.R).

  write_log(sprintf("[hdWGCNA] %s done", celltype))
  list(cor_df = cor_df, hub_genes = hub_genes, modules = modules)
}

# =============================================================================
# 6. Run all cell types (IT-ET Glut first)
# =============================================================================
results_all <- list()
ct_order <- c("IT-ET Glut", setdiff(CT_LEVELS, "IT-ET Glut"))
for (ct in ct_order) {
  results_all[[ct]] <- tryCatch(
    run_hdWGCNA(seurat_human, ct, cells_network, dir_wgcna),
    error = function(e) { write_log(sprintf("ERROR: %s - %s", ct, e$message)); NULL })
}

# =============================================================================
# 7. Combined outputs
# =============================================================================
summary_stats <- lapply(names(results_all), function(ct) {
  if (!is.null(results_all[[ct]])) results_all[[ct]]$cor_df
}) %>% bind_rows()
write.table(summary_stats, file.path(dir_wgcna, "ALL_moduleTraitCor.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(
  lapply(names(results_all), function(ct) {
    if (!is.null(results_all[[ct]])) results_all[[ct]]$hub_genes %>% mutate(celltype = ct)
  }) %>% bind_rows(),
  file.path(dir_wgcna, "ALL_hubGenes.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

# summary dotplot of strongly correlated modules (|r| >= 0.3)
sig_summary <- summary_stats %>%
  pivot_longer(cols = any_of(c("Sex", "Genotype", "Group1", "Group2")),
               names_to = "trait", values_to = "cor") %>%
  filter(!is.na(cor), abs(cor) >= 0.3) %>%
  count(celltype, trait, name = "n_sig") %>%
  mutate(celltype = factor(celltype, levels = rev(sort(unique(celltype)))))

if (nrow(sig_summary) > 0) {
  p_summary <- ggplot(sig_summary, aes(x = trait, y = celltype, size = n_sig, color = trait)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(name = "# sig modules", range = c(2, 10), breaks = c(1, 3, 5, 10)) +
    labs(title = "Significant modules per cell type & trait (hdWGCNA)", x = NULL, y = "Cell type") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 9),
          plot.title  = element_text(size = 13, face = "bold"),
          panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3))
  ggsave(file.path(dir_wgcna, "hdWGCNA_summary_dotplot.pdf"), p_summary,
         width = 7, height = max(8, length(unique(sig_summary$celltype)) * 0.6 + 3))
}
write_log("[hdWGCNA] all cell types done")

# =============================================================================
# Reproducibility
# =============================================================================
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_07_hdWGCNA.txt"))
write_log("=== 07 hdWGCNA complete ===")
cat("[OK] 07 hdWGCNA done ->", dir_wgcna, "\n")
