# ============================================================
# TAM analysis pipeline (DE -> GSEA -> plotting)
# - Seurat object: tissue_Myeloid
# - Required meta columns:
#     subtypes (cell subtype)
#     T_R      (e.g., "T0_Responder", "T0_Non_Responder", ...)
#     sample_name, Response, Timepoint (for proportion plots)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(msigdbr)
  library(fgsea)
  library(tibble)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(tidytext)
  library(patchwork)
  library(EnhancedVolcano)
  library(rstatix)
  library(ggpubr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# -----------------------------
# 0) Inputs / settings
# -----------------------------
cell_types <- c("C1QC_TAM", "GPNMB_TAM", "CCL3_TAM")
timepoints <- c("T0", "T2")

pct_cut <- 0.1
min_size <- 15
max_size <- 500
top_n <- 5

# -----------------------------
# 1) Hallmark gene sets (MSigDB H)
# -----------------------------
hallmark_list <- msigdbr(species = "Homo sapiens", category = "H") |>
  distinct(gs_name, gene_symbol) |>
  group_by(gs_name) |>
  summarise(genes = list(gene_symbol), .groups = "drop") |>
  tibble::deframe()

pretty_hallmark <- function(x, wrap = 26) {
  x |>
    str_replace_all("^HALLMARK_", "") |>
    str_replace_all("_", " ") |>
    str_wrap(wrap)
}

# -----------------------------
# 2) DE helper (one cell type + one timepoint)
# -----------------------------
run_de_one <- function(seu, celltype, tp, pct_cut = 0.1) {
  obj <- subset(seu, subtypes == celltype)
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  
  de <- FindMarkers(
    obj,
    ident.1 = paste0(tp, "_Responder"),
    ident.2 = paste0(tp, "_Non_Responder"),
    group.by = "T_R"
  ) |>
    as.data.frame()
  
  # Apply pct filter when available (Seurat FindMarkers output)
  if (all(c("pct.1", "pct.2") %in% colnames(de))) {
    de <- de |>
      filter(pct.1 > pct_cut, pct.2 > pct_cut)
  }
  
  de
}

# -----------------------------
# 3) fgsea helper (one cell type + one timepoint)
# -----------------------------
run_fgsea_one <- function(seu, celltype, tp, pathways,
                          pct_cut = 0.1,
                          minSize = 15, maxSize = 500) {
  de <- run_de_one(seu, celltype, tp, pct_cut = pct_cut)
  
  if (!"avg_log2FC" %in% colnames(de) || nrow(de) < minSize) {
    return(tibble())
  }
  
  ranks <- de$avg_log2FC
  names(ranks) <- rownames(de)
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  
  if (length(ranks) < minSize) return(tibble())
  
  set.seed(1)
  fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = minSize,
    maxSize = maxSize
  ) |>
    as_tibble() |>
    mutate(
      logFDR = -log10(padj + 1e-300)
    )
}

# -----------------------------
# 4) Build GSEA plot table (Top/Bottom N per cell_type x tp)
# -----------------------------
df_gsea <- crossing(cell_type = cell_types, tp = timepoints) |>
  mutate(res = map2(cell_type, tp, ~run_fgsea_one(
    tissue_Myeloid, .x, .y, hallmark_list,
    pct_cut = pct_cut,
    minSize = min_size,
    maxSize = max_size
  ))) |>
  unnest(res) |>
  filter(!is.na(NES), !is.na(padj))

df_plot <- df_gsea |>
  group_by(cell_type, tp) |>
  slice_max(order_by = NES, n = top_n, with_ties = FALSE) |>
  bind_rows(
    df_gsea |>
      group_by(cell_type, tp) |>
      slice_min(order_by = NES, n = top_n, with_ties = FALSE)
  ) |>
  ungroup() |>
  mutate(
    pathway_lab = pretty_hallmark(pathway, wrap = 26),
    cell_type = factor(cell_type, levels = cell_types),
    tp = factor(tp, levels = timepoints)
  )

max_abs <- max(abs(df_plot$NES), na.rm = TRUE)

# -----------------------------
# 5) Plot: GSEA for a single timepoint (facet labels on top)
#    - Dot is shown only for significant pathways (padj < 0.05)
#    - Pathways are ordered by NES within each facet
# -----------------------------
plot_gsea_timepoint <- function(df_plot, tp_use, nes_lim = NULL,
                                title_prefix = "Ranked GSEA (NR vs R)") {
  df_tp <- df_plot |>
    filter(tp == tp_use) |>
    mutate(
      # Order y-axis by NES within each facet (cell_type)
      pathway_lab = tidytext::reorder_within(pathway_lab, NES, cell_type)
    )
  
  lim_use <- nes_lim %||% c(-max_abs, max_abs)
  
  ggplot(df_tp, aes(x = NES, y = pathway_lab)) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey55") +
    geom_col(aes(fill = NES), width = 0.72, alpha = 0.95) +
    geom_point(
      data = df_tp |> filter(padj < 0.05),
      aes(size = logFDR),
      shape = 21, fill = "white", color = "grey10", stroke = 0.5
    ) +
    facet_wrap(~ cell_type, ncol = 1, scales = "free_y", strip.position = "top") +
    tidytext::scale_y_reordered() +
    scale_fill_gradient2(
      low = "#2C6AA6", mid = "white", high = "#B34A4A",
      midpoint = 0, limits = lim_use, oob = scales::squish,
      name = "NES",
      guide = guide_colorbar(barheight = unit(35, "mm"), barwidth = unit(4, "mm"))
    ) +
    scale_size_binned(
      name = expression(-log[10]("FDR")),
      range = c(1.6, 4.2),
      guide = guide_legend(override.aes = list(fill = "white"))
    ) +
    labs(
      title = paste0(title_prefix, " | ", tp_use),
      x = "NES",
      y = NULL,
      caption = "Dots are shown only for pathways with FDR < 0.05."
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, margin = margin(b = 6)),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(margin = margin(t = 6)),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(face = "bold", size = 11),
      panel.spacing.y = unit(7, "mm"),
      legend.position = "right",
      legend.box = "vertical"
    )
}

p_T0 <- plot_gsea_timepoint(df_plot, "T0", nes_lim = c(-max_abs, max_abs))
p_T2 <- plot_gsea_timepoint(df_plot, "T2", nes_lim = c(-max_abs, max_abs))

p_gsea <- (p_T0 | p_T2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Print / save
p_gsea

# -----------------------------
# 6) Volcano plot helper (one cell type + one timepoint)
# -----------------------------
plot_volcano <- function(seu, celltype, tp,
                         pct_cut = 0.1,
                         fc_cut = 0.5,
                         p_cut = 0.01,
                         select_lab = NULL) {
  de <- run_de_one(seu, celltype, tp, pct_cut = pct_cut)
  
  EnhancedVolcano(
    de,
    lab = rownames(de),
    x = "avg_log2FC",
    y = "p_val_adj",
    FCcutoff = fc_cut,
    pCutoff = p_cut,
    drawConnectors = TRUE,
    selectLab = select_lab
  ) +
    ggtitle(paste(celltype, tp, sep = " | "))
}

# Example:
# p_vol <- plot_volcano(
#   tissue_Myeloid, "C1QC_TAM", "T0",
#   select_lab = c("ALOX5AP","NFKBIZ","FCGR2B","NEAT1")
# )
# p_vol

# -----------------------------
# 7) Proportion plots: NR vs R within each timepoint (T0 and T2)
# -----------------------------
make_prop_df <- function(seu,
                         keep_timepoints = c("T0", "T2"),
                         exclude_subtypes = NULL) {
  meta <- seu@meta.data |>
    filter(Timepoint %in% keep_timepoints)
  
  if (!is.null(exclude_subtypes)) {
    meta <- meta |> filter(!subtypes %in% exclude_subtypes)
  }
  
  sample_counts <- meta |>
    group_by(sample_name) |>
    summarise(total_cells = n(), .groups = "drop")
  
  subtype_counts <- meta |>
    group_by(sample_name, subtypes) |>
    summarise(cell_count = n(), .groups = "drop")
  
  subtype_counts |>
    left_join(sample_counts, by = "sample_name") |>
    mutate(proportion = cell_count / total_cells) |>
    left_join(
      meta |> select(sample_name, Response, Timepoint, patients) |> distinct(),
      by = "sample_name"
    ) |>
    filter(!is.na(proportion), !is.na(Response), !is.na(Timepoint)) |>
    mutate(
      Timepoint = factor(Timepoint, levels = keep_timepoints),
      Response = factor(Response, levels = c("Non_Responder", "Responder"))
    )
}
plot_prop_by_timepoint <- function(df_prop, tp_use, show_p = TRUE) {
  
  df_tp <- df_prop |> filter(Timepoint == tp_use)
  
  # Wilcoxon test per subtype: NR vs R
  cap <- df_tp |>
    group_by(subtypes) |>
    summarise(y_max = max(proportion, na.rm = TRUE), .groups = "drop")
  
  stat_tp <- df_tp |>
    group_by(subtypes) |>
    filter(n_distinct(Response) == 2) |>
    wilcox_test(proportion ~ Response, paired = FALSE) |>
    ungroup() |>
    left_join(cap, by = "subtypes") |>
    mutate(
      xmin = "Non_Responder",
      xmax = "Responder",
      y.position = y_max * 1.15,
      label = paste0("p=", signif(p, 2))
    )
  
  cols_resp <- c("Non_Responder" = "#FC8D62", "Responder" = "#66C2A5")
  
  p <- ggplot(df_tp, aes(x = Response, y = proportion, fill = Response)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, color = "black") +
    geom_jitter(width = 0.12, size = 1.1, alpha = 0.75) +
    facet_wrap(~ subtypes, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = cols_resp) +
    theme_bw(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(face = "bold")
    ) +
    labs(title = tp_use, x = "Response", y = "Cell proportion")
  
  if (show_p && nrow(stat_tp) > 0) {
    p <- p +
      ggpubr::stat_pvalue_manual(
        stat_tp,
        label = "label",
        xmin = "xmin", xmax = "xmax",
        y.position = "y.position",
        tip.length = 0.01,
        bracket.size = 0.4,
        size = 3
      )
  }
  
  p
}

# Build df_prop and plot T0/T2
df_prop <- make_prop_df(
  tissue_Myeloid,
  keep_timepoints = c("T0", "T2"),
  exclude_subtypes = c("pDC", "mregDC_LAMP3(low)")
)

p_prop_T0 <- plot_prop_by_timepoint(df_prop, "T0", show_p = TRUE)
p_prop_T2 <- plot_prop_by_timepoint(df_prop, "T2", show_p = TRUE)

p_prop <- p_prop_T0 / p_prop_T2
p_prop

