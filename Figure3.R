# Cell type specific signature from : ref :: Bi, K., Kashima, S., Camp, S. Y., Meli, K., Saad, E., Titchen, B. M., ... & Van Allen, E. M. (2025). Myeloid cells mediate interferon-driven resistance to immunotherapy in advanced renal cell carcinoma. Immunity, 58(11), 2814-2829.
# code ref https://github.com/kevinbi2599/ccRCC_IFN_Resistance_Immunity2025

library(dplyr)
library(tidyr)
library(purrr)
library(fgsea)
library(ggplot2)
library(stringr)
library(tidytext)
library(scales)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
T0_DE <- filter(T0_R_NR,pct.1 > 0.1 & pct.2 > 0.1)
T2_DE <- filter(T2_R_NR,pct.1 > 0.1 & pct.2 > 0.1)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(forcats)
library(scales)
library(forcats)


deg_list <- list(
  T0_R_NR = T0_DE,
  T2_R_NR = T2_DE
) |>
  imap(\(df, nm) {
    df |>
      rename(gene = ...1) |>
      mutate(
        contrast = nm,
        timepoint = if_else(grepl("^T0", nm), "T0", "T2")
      )
  })

split_deg_direction <- function(df, up_cutoff = 0.5, down_cutoff = -0.5) {
  bind_rows(
    df |>
      filter(avg_log2FC >= up_cutoff, p_val_adj < 0.01) |>
      mutate(direction = "R"),
    df |>
      filter(avg_log2FC <= down_cutoff, p_val_adj < 0.01) |>
      mutate(direction = "NR")
  ) |>
    distinct(gene, direction, .keep_all = TRUE)
}

deg_split <- deg_list |>
  imap(\(df, nm) split_deg_direction(df)) |>
  list_rbind()

# ---- 2) Reactome enrichment 함수 ----
ratio_to_numeric <- function(x) {
  parts <- strsplit(x, "/", fixed = TRUE)
  vapply(parts, \(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))
}

run_enrich_reactome <- function(genes_symbol) {
  genes_symbol <- genes_symbol |>
    unique() |>
    toupper()
  
  if (length(genes_symbol) == 0) return(NULL)
  
  mapped <- suppressWarnings(
    clusterProfiler::bitr(
      genes_symbol,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )
  
  if (is.null(mapped) || nrow(mapped) == 0) return(NULL)
  
  entrez_ids <- mapped$ENTREZID |>
    unique() |>
    na.omit()
  
  if (length(entrez_ids) == 0) return(NULL)
  
  er <- ReactomePA::enrichPathway(
    gene = entrez_ids,
    organism = "human",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  if (is.null(er) || nrow(er@result) == 0) return(NULL)
  
  as_tibble(er@result) |>
    transmute(
      term = Description,
      gene_ratio = ratio_to_numeric(GeneRatio),
      count = Count,
      p_adjust = p.adjust
    )
}

enrich_reactome <- deg_split |>
  group_by(timepoint, direction) |>
  summarise(genes = list(unique(gene)), .groups = "drop") |>
  mutate(res = map(genes, run_enrich_reactome)) |>
  select(-genes) |>
  unnest(res)

top_n <- 5
saveRDS(enrich_reactome,'enrich_reactome_Tissue_T.rds')

reactome_top <- enrich_reactome %>%
  filter(p_adjust < 0.05) %>%
  # filter(gene_ratio  >= 0.1) %>%
  group_by(timepoint, direction) %>%
  arrange(p_adjust, desc(gene_ratio), desc(count)) %>%
  slice_head(n = top_n) %>%
  ungroup()

reactome_top <- reactome_top |>
  mutate(
    panel = interaction(timepoint, direction, drop = TRUE),
    term_panel = paste(term, panel, sep = "___"),
    term_panel = fct_reorder(term_panel, gene_ratio)
  )

# ---- 6) Plot: row = timepoint, col = NR/R ----
low_col <- "#2166AC"   # deep red
high_col <- "#B2182B"  # deep blue
saveRDS(reactome_top,'CD8_DE_reactome.rds')

library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(scales)

# reactome_top: columns = timepoint, direction(NR/R), term, gene_ratio, count, p_adjust
reactome_plot_df <- reactome_top |>
  mutate(
    gene_ratio_signed = if_else(direction == "NR", -gene_ratio, gene_ratio),
    term_wrapped = str_wrap(term, width = 45)
  ) |>
  group_by(timepoint) |>
  # timepoint별로 y축 정렬 (좌/우까지 포함한 signed ratio 기준)
  mutate(term_wrapped = fct_reorder(term_wrapped, gene_ratio_signed)) |>
  ungroup()

low_col <-  "#2166AC"  # red
high_col <- "#B2182B"  # blue

ggplot(reactome_plot_df) +
  geom_vline(xintercept = 0, color = "grey35", linewidth = 0.4) +
  geom_point(
    aes(x = gene_ratio_signed, y = term_wrapped, color = p_adjust, size = count),
    alpha = 0.95
  ) +
  facet_wrap(~ timepoint, ncol = 1, scales = "free_y") +
  scale_x_continuous(labels = \(x) abs(x), breaks = pretty_breaks(5)) +
  scale_color_gradient(low = low_col, high = high_col, trans = "reverse") +
  scale_size(range = c(2.5, 7)) +
  labs(
    x = "Gene ratio (NR ← 0 → R)",
    y = NULL,
    title = "Reactome enrichment of DEGs",
    color = "Adj. p",
    size = "Count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )


library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(forcats)

ratio_to_numeric <- function(x) {
  parts <- strsplit(x, "/", fixed = TRUE)
  vapply(parts, \(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))
}

run_enrich_go_bp <- function(genes_symbol) {
  genes_symbol <- genes_symbol |>
    unique() |>
    toupper()
  
  if (length(genes_symbol) == 0) return(NULL)
  
  ego <- clusterProfiler::enrichGO(
    gene = genes_symbol,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05
  )
  
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  
  as_tibble(ego@result) |>
    transmute(
      term = Description,
      gene_ratio = ratio_to_numeric(GeneRatio),
      count = Count,
      p_adjust = p.adjust
    )
}

enrich_go <- deg_split |>
  group_by(timepoint, direction) |>
  summarise(genes = list(unique(gene)), .groups = "drop") |>
  mutate(res = map(genes, run_enrich_go_bp)) |>
  select(-genes) |>
  unnest(res)

top_n <- 5

go_top <- enrich_go |>
  filter(p_adjust < 0.05 ) |>
  group_by(timepoint, direction) |>
  arrange(p_adjust, desc(count)) |>
  slice_head(n = top_n) |>
  ungroup()

go_top <- go_top |>
  mutate(
    panel = interaction(timepoint, direction, drop = TRUE),
    term_panel = paste(term, panel, sep = "___"),
    term_panel = fct_reorder(term_panel, gene_ratio)
  )

# ---- Plot: row=timepoint, col=NR/R ----
low_col <- "#2166AC"   # deep red
high_col <- "#B2182B"  # deep blue

library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(scales)

go_plot_df <- go_top |>
  mutate(
    gene_ratio_signed = if_else(direction == "NR", -gene_ratio, gene_ratio),
    term_wrapped = str_wrap(term, width = 45)
  ) |>
  group_by(timepoint) |>
  mutate(term_wrapped = fct_reorder(term_wrapped, gene_ratio_signed)) |>
  ungroup()

ggplot(go_plot_df) +
  geom_vline(xintercept = 0, color = "grey35", linewidth = 0.4) +
  geom_point(
    aes(
      x = gene_ratio_signed,
      y = term_wrapped,
      size = count,
      color = p_adjust
    ),
    alpha = 0.95
  ) +
  facet_wrap(~ timepoint, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    labels = \(x) abs(x),
    breaks = pretty_breaks(5)
  ) +
  scale_color_gradient(
    low = low_col,
    high = high_col,
    trans = "reverse"
  ) +
  scale_size(range = c(2.5, 7)) +
  labs(
    x = "Gene ratio (NR ← 0 → R)",
    y = NULL,
    title = "GO (BP) enrichment of DEGs",
    color = "Adj. p",
    size = "Count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )


# regulons_asGeneSet_CD8_Regulon : named list (TF -> target genes)
# T0_DE : character vector (DEG genes)

deg <- unique(T0_DE)

res_overlap <- lapply(names(regulons_asGeneSet_CD8_Regulon), function(tf){
  targets <- unique(regulons_asGeneSet_CD8_Regulon[[tf]])
  hit <- intersect(deg, targets)
  
  data.frame(
    TF = tf,
    n_overlap = length(hit),
    overlap_genes = paste(hit, collapse = ", "),
    stringsAsFactors = FALSE
  )
})

res_overlap <- do.call(rbind, res_overlap)
res_overlap2 <- subset(res_overlap, n_overlap > 0)
res_overlap2 <- res_overlap2[order(res_overlap2$n_overlap, decreasing = TRUE), ]

# inputs:
# regulons_asGeneSet_CD8_Regulon : named list (TF -> target genes)
# T0_DE : character vector

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  inter <- length(intersect(a, b))
  uni   <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

deg <- unique(T0_DE)

tf_names <- names(regulons_asGeneSet_CD8_Regulon)
jac_vec <- vapply(tf_names, function(tf){
  jaccard(deg, regulons_asGeneSet_CD8_Regulon[[tf]])
}, numeric(1))

jac_mat <- matrix(jac_vec, ncol = 1)
rownames(jac_mat) <- tf_names
colnames(jac_mat) <- "T0_DE"

top_n <- 50
keep <- order(jac_vec, decreasing = TRUE)[seq_len(min(top_n, length(jac_vec)))]
jac_mat_top <- jac_mat[keep, , drop = FALSE]

# heatmap (pheatmap)
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
pheatmap::pheatmap(jac_mat_top, cluster_rows = TRUE, cluster_cols = FALSE)


# install.packages("ComplexHeatmap")
# install.packages("circlize")

library(ComplexHeatmap)
library(circlize)

# jac_mat_top: TF x 1 (또는 TF x timepoints) matrix
# rownames(jac_mat_top) = TF, colnames(jac_mat_top) = "T0_DE" 등

ht <- Heatmap(
  jac_mat_top,
  name = "Jaccard index",                      
  column_title = "Jaccard index (DEG vs regulon)", 
  row_title = "TF regulon",
  cluster_rows = TRUE,
  cluster_columns = FALSE,                       
  heatmap_legend_param = list(title = "Jaccard index")
)

draw(ht)




suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(lme4)
  library(lmerTest)
  library(UCell)
})

# ----------------------------
# 0) INPUTS (based on your object)
# ----------------------------
seu <- RCC_IO

celltype_col <- "cell_type"   # broad cell type
patient_col  <- "patients"    # patient ID
covars <- c("nFeature_RNA", "histology")

# ATF3 candidate target pool (SCENIC regulon gene set)
atf3_targets <- regulons_asGeneSet_CD8_Regulon$ATF3_extended

# Target broad cell type for CTS derivation (T/NK as target)
target_celltype <- "T/NK cell"

exclude_celltypes <- character(0)

q_cut <- 0.05      # inducible BH q
pct2_cut <- 0.20   # comparator pct < 20%
p_cut <- 0.05      # Bonferroni p

# ----------------------------
# 1) Define the putative pool
# ----------------------------
goi <- intersect(rownames(seu), atf3_targets)
message("ATF3 target pool (goi) size: ", length(goi))

# ----------------------------
# 2) Build an ATF3 activity score
#    - If SCENIC AUC is already present in meta.data, using that is preferable.
#    - Here, we use a gene set-based UCell score as a proxy for ATF3 activity.
# ----------------------------
seu <- AddModuleScore_UCell(seu, features = list(ATF3 = goi), name = "ATF3_UCell")

# Identify the exact meta.data column name created by UCell
score_col <- grep("^ATF3_UCell", colnames(seu@meta.data), value = TRUE)[1]
if (is.na(score_col) || length(score_col) == 0) {
  stop("Could not find ATF3_UCell score column in meta.data. Check AddModuleScore_UCell output.")
}
message("Using ATF3 activity score column: ", score_col)

# ----------------------------
# 3) ATF3-inducible genes (mixed model in T/NK)
# ----------------------------
seu_tnk <- subset(seu, subset = .data[[celltype_col]] == target_celltype)

seu_tnk$ATF3_z <- as.numeric(scale(seu_tnk@meta.data[[score_col]]))

expr <- GetAssayData(seu_tnk, slot = "data")[goi, , drop = FALSE]
meta <- seu_tnk@meta.data %>%
  mutate(cell = rownames(.)) %>%
  transmute(
    cell = cell,
    atf3 = ATF3_z,
    patient = .data[[patient_col]],
    nFeature_RNA = .data[["nFeature_RNA"]],
    histology = .data[["histology"]]
  )

# Use only covariates that actually exist in meta
covars_use <- covars[covars %in% colnames(meta)]

build_formula <- function() {
  rhs <- c("atf3", covars_use, sprintf("(1|%s)", "patient"))
  as.formula(paste0("y ~ ", paste(rhs, collapse = " + ")))
}
form <- build_formula()

fit_one_gene <- function(g) {
  df <- meta
  df$y <- as.numeric(expr[g, df$cell])

  # Skip genes with no variance
  if (sd(df$y, na.rm = TRUE) == 0 || all(is.na(df$y))) {
    return(tibble(gene = g, beta_atf3 = NA_real_, p_atf3 = NA_real_))
  }

  m <- tryCatch(
    lmer(form, data = df),
    error = function(e) NULL
  )
  if (is.null(m)) {
    return(tibble(gene = g, beta_atf3 = NA_real_, p_atf3 = NA_real_))
  }

  co <- summary(m)$coef
  if (!("atf3" %in% rownames(co))) {
    return(tibble(gene = g, beta_atf3 = NA_real_, p_atf3 = NA_real_))
  }

  tibble(
    gene = g,
    beta_atf3 = co["atf3", "Estimate"],
    p_atf3 = co["atf3", "Pr(>|t|)"]
  )
}

message("Fitting mixed models in T/NK for inducibility... (", length(goi), " genes)")
ind_tbl <- map_dfr(goi, fit_one_gene) %>%
  mutate(q_atf3 = p.adjust(p_atf3, method = "BH"))

atf3_inducible <- ind_tbl %>%
  filter(!is.na(beta_atf3), beta_atf3 > 0, !is.na(q_atf3), q_atf3 < q_cut) %>%
  pull(gene)

message("ATF3-inducible genes in T/NK: ", length(atf3_inducible))

# ----------------------------
# 4) T/NK CTS genes (one-vs-each, strict)
#    - target: "T/NK cell"
#    - compare: all other cell types
#    - each comparator must pass cutoffs
#    - LR test + latent vars to adjust for patient/batch/timepoint
# ----------------------------
Idents(seu) <- seu@meta.data[[celltype_col]]

all_types <- setdiff(levels(factor(seu@meta.data[[celltype_col]])), exclude_celltypes)
comp_types <- setdiff(all_types, target_celltype)

latent_vars <- c(patient_col,covars, intersect(c("Timepoint"), colnames(seu@meta.data)))
latent_vars <- unique(latent_vars)

run_one_vs_each <- function(comp) {
  cells_1 <- WhichCells(seu, idents = target_celltype)
  cells_2 <- WhichCells(seu, idents = comp)

  sub <- subset(seu, cells = c(cells_1, cells_2))
  sub$grp <- ifelse(sub@meta.data[[celltype_col]] == target_celltype, "TNK", comp)
  Idents(sub) <- "grp"

  mk <- FindMarkers(
    sub,
    ident.1 = "TNK",
    ident.2 = comp,
    test.use = "LR",
    latent.vars = latent_vars,
    features = goi,
    min.pct = 0,
    logfc.threshold = 0,
    only.pos = TRUE
  ) %>%
    rownames_to_column("gene") %>%
    mutate(comp = comp)

  # Explicitly enforce Bonferroni correction over the goi set only
  mk <- mk %>%
    mutate(p_bonf_goi = p.adjust(p_val, method = "bonferroni", n = length(goi)))

  mk
}

message("Running one-vs-each CTS tests for target = ", target_celltype)
cts_tbl <- map_dfr(comp_types, run_one_vs_each)

cts_pass <- cts_tbl %>%
  mutate(
    pass = avg_log2FC > 0 &
      (pct.1 - pct.2) > 0 &
      pct.2 < pct2_cut &
      p_bonf_goi < p_cut
  )

tnk_cts <- cts_pass %>%
  group_by(gene) %>%
  summarise(
    n_pass = sum(pass, na.rm = TRUE),
    n_test = n_distinct(comp),
    .groups = "drop"
  ) %>%
  filter(n_pass == n_test) %>%
  pull(gene)

message("T/NK CTS genes (strict, pass all comparators): ", length(tnk_cts))

# ----------------------------
# 5) Final CTS ATF3 signature = intersection
# ----------------------------
atf3_tnk_cts_signature <- Reduce(intersect, list(goi, atf3_inducible, tnk_cts))
message("FINAL ATF3 T/NK CTS signature size: ", length(atf3_tnk_cts_signature))

# ----------------------------
# 6) Score the final signature in all cells
#    - UCell score computed on normalized data
# ----------------------------
seu <- AddModuleScore_UCell(
  seu,
  features = list(ATF3_TNK_CTS = atf3_tnk_cts_signature),
  name = "ATF3_TNK_CTS_UCell"
)

# Identify the exact meta.data column name created for the final signature
final_score_col <- grep("^ATF3_TNK_CTS_UCell", colnames(seu@meta.data), value = TRUE)[1]
if (is.na(final_score_col) || length(final_score_col) == 0) {
  stop("Could not find ATF3_TNK_CTS_UCell score column in meta.data. Check AddModuleScore_UCell output.")
}
message("Using final signature score column: ", final_score_col)

# Inspect results
seu@meta.data %>%
  select(all_of(c(celltype_col, patient_col, score_col, final_score_col))) %>%
  head() %>%
  print()

# Save gene sets / tables if needed
out <- list(
  goi = goi,
  inducibility_table = ind_tbl,
  atf3_inducible = atf3_inducible,
  cts_table = cts_tbl,
  tnk_cts = tnk_cts,
  final_signature = atf3_tnk_cts_signature
)

# Example: save to disk
saveRDS(out, file = "ATF3_Expanded_re_TNK_CTS_signature_objects.rds")
# write.table(atf3_tnk_cts_signature, "ATF3_TNK_CTS_signature_genes.txt",
#             quote = FALSE, row.names = FALSE, col.names = FALSE)

RCC_IO <- AddModuleScore_UCell(
  RCC_IO,
  features = list(
    ATF3_Extended_orig = regulons_asGeneSet_CD8_Regulon$ATF3_extended,
    ATF3_Extended_CTS = out$final_signature
  )
)


#Survival analysis general 
GEX <- read_excel("JAVELIN_Renal_101/41591_2020_1044_MOESM3_ESM.xlsx", 
                                                 sheet = "S13_Gene_expression_TPM", skip = 1)
Javelin_Clinical_info <- read_excel("JAVELIN_Renal_101/41591_2020_1044_MOESM3_ESM.xlsx", 
                                           sheet = "S11_Clinical_data", skip = 1)


library(dplyr)
GEX <- data.frame(GEX)
Javelin_Clinical_info <- data.frame(Javelin_Clinical_info)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(GSVA)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(RColorBrewer)
})

## ===============================================================
## 0) Helpers
## ===============================================================

norm_chr <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\s+", " ")
  str_trim(x)
}

clean_geneset <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x)]
  x[nzchar(x)]
}

# If geneID column exists -> median collapse duplicates
# Else assume rownames are already gene symbols
prep_expr_matrix_gene_median <- function(gex_df, gene_id_col = "geneID") {
  gex <- as.data.frame(gex_df, check.names = FALSE)
  
  if (gene_id_col %in% names(gex)) {
    gex <- gex |>
      mutate(geneID = norm_chr(.data[[gene_id_col]])) |>
      filter(!is.na(geneID), geneID != "")
    
    num_cols <- setdiff(names(gex), gene_id_col)
    num_cols <- num_cols[vapply(gex[num_cols], is.numeric, logical(1))]
    if (length(num_cols) < 3) stop("numeric sample columns too few after excluding geneID")
    
    gex_med <- gex |>
      group_by(geneID) |>
      summarise(across(all_of(num_cols), ~ stats::median(.x, na.rm = TRUE)), .groups = "drop")
    
    mat <- as.matrix(gex_med[, num_cols, drop = FALSE])
    rownames(mat) <- gex_med$geneID
    storage.mode(mat) <- "numeric"
    stopifnot(!anyDuplicated(rownames(mat)))
    return(mat)
  }
  
  # no geneID column
  mat <- as.matrix(gex)
  storage.mode(mat) <- "numeric"
  if (is.null(rownames(mat)) || any(!nzchar(rownames(mat)))) {
    stop("GEX has no geneID column and rownames are missing/empty. Provide geneID or set rownames first.")
  }
  mat
}

cox_label <- function(df2, time_col, event_col, group_col = ".group") {
  df3 <- df2 |>
    mutate(
      .time  = as.numeric(.data[[time_col]]),
      .event = as.numeric(.data[[event_col]]),
      .group = factor(.data[[group_col]], levels = c("Low", "High"))
    ) |>
    filter(!is.na(.time), !is.na(.event), !is.na(.group)) |>
    filter(.event %in% c(0, 1))
  
  if (n_distinct(df3$.group) < 2) return("HR=NA")
  
  fit <- survival::coxph(survival::Surv(.time, .event) ~ .group, data = df3)
  s <- summary(fit)
  
  hr <- s$coefficients[1, "exp(coef)"]
  p  <- s$coefficients[1, "Pr(>|z|)"]
  ci <- s$conf.int[1, c("lower .95", "upper .95")]
  
  paste0(
    "HR=", sprintf("%.2f", hr),
    " (95% CI ", sprintf("%.2f", ci[1]), "–", sprintf("%.2f", ci[2]), ")",
    "\nCox p=", format.pval(p, digits = 2, eps = 1e-4)
  )
}

km_plot_cox <- function(df, time_col, event_col, group_col, title_text,
                        show_risktable = TRUE, palette = "Set2") {
  pal <- RColorBrewer::brewer.pal(8, palette)[c(2, 4)]
  names(pal) <- c("Low", "High")
  
  df2 <- df |>
    mutate(
      .time  = as.numeric(.data[[time_col]]),
      .event = as.numeric(.data[[event_col]]),
      .group = factor(.data[[group_col]], levels = c("Low", "High"))
    ) |>
    filter(!is.na(.time), !is.na(.event), !is.na(.group)) |>
    filter(.event %in% c(0, 1))
  
  fit <- survival::survfit(survival::Surv(.time, .event) ~ .group, data = df2)
  hr_lab <- cox_label(df2, time_col = ".time", event_col = ".event", group_col = ".group")
  
  g <- survminer::ggsurvplot(
    fit,
    data = df2,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = show_risktable,
    risk.table.height = 0.26,
    legend.title = NULL,
    legend.labs = c("Low", "High"),
    palette = unname(pal),
    title = title_text,
    xlab = "Time",
    ylab = "Progression-free probability",
    censor.shape = 124,
    censor.size  = 2.6,
    ggtheme = theme_classic(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "top",
        axis.title = element_text(face = "bold")
      ),
    tables.theme = theme_classic(base_size = 10)
  )
  
  g$plot <- g$plot +
    annotate("text", x = Inf, y = Inf, label = hr_lab,
             hjust = 1.05, vjust = 1.15, size = 3.1, fontface = "bold")
  
  g
}

run_gsva_wide <- function(expr_mat, signatures, method = "gsva", prefix = "GSVA_") {
  gset <- lapply(signatures, clean_geneset)
  gs <- GSVA::gsva(expr_mat, gset, method = method)
  
  as.data.frame(t(gs), check.names = FALSE) |>
    rownames_to_column("sample") |>
    rename_with(~ paste0(prefix, .x), -sample)
}

## ===============================================================
## 1) Main: TRT01P strata survival (PFS)
## ===============================================================

run_javelin_pfs_by_trt <- function(
    gex_df,
    clin_df,
    signatures,
    id_col = "ID",
    trt_col = "TRT01P",
    time_col = "PFS_P",
    censor_col = "PFS_P_CNSR",
    gsva_method = "gsva",
    show_risktable = TRUE,
    palette = "Set2",
    cohort_label = "JAVELIN",
    min_n = 20
) {
  expr_mat <- prep_expr_matrix_gene_median(gex_df, gene_id_col = "geneID")
  
  clin <- as.data.frame(clin_df, check.names = FALSE) |>
    transmute(
      sample = norm_chr(.data[[id_col]]),
      TRT01P = factor(norm_chr(.data[[trt_col]])),
      PFS_time = suppressWarnings(as.numeric(.data[[time_col]])),
      cens = suppressWarnings(as.numeric(.data[[censor_col]]))
    ) |>
    mutate(
      # event=1 means progression/death event
      PFS_event = dplyr::case_when(
        cens %in% c(0, 1) ~ 1 - cens,
        TRUE ~ NA_real_
      )
    ) |>
    filter(!is.na(sample), sample != "",
           !is.na(TRT01P),
           !is.na(PFS_time),
           PFS_event %in% c(0, 1)) |>
    distinct(sample, .keep_all = TRUE)
  
  # overall match
  common_all <- intersect(colnames(expr_mat), clin$sample)
  message("Total matched samples with usable PFS: ", length(common_all))
  if (length(common_all) < 10) stop("Too few matched samples overall.")
  
  clin <- clin |>
    filter(sample %in% common_all)
  
  out <- list()
  stats <- list()
  
  for (trt in levels(clin$TRT01P)) {
    clin_t <- clin |>
      filter(TRT01P == trt)
    
    common <- intersect(colnames(expr_mat), clin_t$sample)
    message("[", trt, "] matched: ", length(common))
    
    if (length(common) < min_n) {
      message("[", trt, "] < min_n (", min_n, "), skip")
      next
    }
    
    clin_t2 <- clin_t |>
      filter(sample %in% common) |>
      arrange(match(sample, colnames(expr_mat)))
    
    expr2 <- expr_mat[, clin_t2$sample, drop = FALSE]
    stopifnot(all(colnames(expr2) == clin_t2$sample))
    
    scores <- run_gsva_wide(expr2, signatures, method = gsva_method, prefix = "GSVA_")
    dat <- scores |>
      left_join(clin_t2 |> select(sample, PFS_time, PFS_event, TRT01P), by = "sample")
    
    sigs <- names(signatures)
    plots <- setNames(vector("list", length(sigs)), sigs)
    
    for (s in sigs) {
      col_s <- paste0("GSVA_", s)
      
      d_s <- dat |>
        filter(!is.na(.data[[col_s]])) |>
        mutate(
          Group = ifelse(.data[[col_s]] >= median(.data[[col_s]], na.rm = TRUE), "High", "Low"),
          Group = factor(Group, levels = c("Low", "High"))
        )
      
      plots[[s]] <- km_plot_cox(
        d_s,
        time_col = "PFS_time",
        event_col = "PFS_event",
        group_col = "Group",
        title_text = paste0("[", cohort_label, " | ", trt, "] PFS | ", s, " (median split)"),
        show_risktable = show_risktable,
        palette = palette
      )
    }
    
    out[[as.character(trt)]] <- list(data = dat, plots = plots)
    stats[[as.character(trt)]] <- tibble(TRT01P = as.character(trt), n = nrow(dat))
  }
  
  list(results = out, summary = bind_rows(stats))
}

save_trt_surv_pdf <- function(res, file, width = 12, height = 8) {
  pdf(file, width = width, height = height, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  
  for (trt in names(res$results)) {
    for (s in names(res$results[[trt]]$plots)) {
      print(res$results[[trt]]$plots[[s]])
    }
  }
}

## ===============================================================
## 2) RUN (사용 예시)
## ===============================================================

ATF3 <- list(ATF3_TCS = ATF3_Expanded_re_TNK_CTS_signature_objects$final_signature)
rownames(GEX) <- GEX$HUGO
GEX <- GEX[,-1]
res_jav <- run_javelin_pfs_by_trt(
  gex_df = GEX,
  clin_df = Javelin_Clinical_info,
  signatures = ATF3,
  id_col = "ID",
  trt_col = "TRT01P",
  time_col = "PFS_P",
  censor_col = "PFS_P_CNSR",
  gsva_method = "gsva",
  palette = "Set2",
  cohort_label = "JAVELIN",
  min_n = 20
)

print(res_jav$summary)
res_jav$results$`Avelumab+Axitinib`
res_jav$results$Sunitinib

save_trt_surv_pdf(
  res_jav,
  file = "JAVELIN_GSVA_PFS_by_TRT01P_ATF3.pdf",
  width = 6,
  height = 8
)




suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(readxl)
  library(readr)
  library(GSVA)
  library(survival)
  library(ggplot2)
  library(forcats)
  library(Biobase)
})

# ===============================================================
# Main ATF3 continuous Cox workflow across cohorts
# - Continuous Cox model per cohort/subgroup
# - Separate forest plots for PFS and OS
# - Cohort-level GSVA cache to avoid recomputation
# - JAVELIN: arm-specific only
# - UC_Genome: Immunotherapy yes / Chemotherapy yes + Immunotherapy not yes
# - HUGABOOM: PFS
# ===============================================================

# ===============================================================
# 0) Signature and output directory
# ===============================================================

signatures <- list(
  ATF3_final = ATF3_Expanded_re_TNK_CTS_signature_objects$final_signature
)

out_dir <- "./ATF3_continuous_cox_main"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cache_dir <- file.path(out_dir, "gsva_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# ===============================================================
# 1) General helpers
# ===============================================================

norm_chr <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\s+", " ")
  str_trim(x)
}

is_yes <- function(x) {
  x2 <- str_to_upper(norm_chr(x))
  !is.na(x2) & x2 == "YES"
}

clean_geneset <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x)]
  x[nzchar(x)]
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

clean_dimnames <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

ensure_matrix_dimnames <- function(mat, cohort_label) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  
  rn <- clean_dimnames(rownames(mat))
  cn <- clean_dimnames(colnames(mat))
  
  rownames(mat) <- rn
  colnames(mat) <- cn
  
  if (is.null(rn) || any(is.na(rn))) {
    stop(cohort_label, ": expression matrix must have non-empty gene-symbol rownames.")
  }
  
  if (is.null(cn) || any(is.na(cn))) {
    stop(cohort_label, ": expression matrix must have non-empty sample colnames.")
  }
  
  if (anyDuplicated(rn)) {
    stop(cohort_label, ": duplicated gene-symbol rownames detected.")
  }
  
  if (anyDuplicated(cn)) {
    stop(cohort_label, ": duplicated sample colnames detected.")
  }
  
  mat
}

extract_uc_sample4 <- function(x) {
  x <- as.character(x)
  m <- regexpr("217-[0-9]{4}", x)
  hit <- ifelse(m > 0, regmatches(x, m), NA_character_)
  ifelse(!is.na(hit), sub("^217-", "", hit), NA_character_)
}

normalize_unc108_rs <- function(x) {
  x <- norm_chr(x)
  x <- gsub("\\.", "-", x)
  x
}

cache_path <- function(cohort_name) {
  file.path(cache_dir, paste0("gsva_", cohort_name, ".rds"))
}

read_or_build_gsva_data <- function(cache_file, builder_fun, use_cache = TRUE, force_recompute = FALSE) {
  if (use_cache && !force_recompute && file.exists(cache_file)) {
    message("Loading cached GSVA data: ", cache_file)
    return(readRDS(cache_file))
  }
  
  dat <- builder_fun()
  saveRDS(dat, cache_file)
  message("Saved GSVA data: ", cache_file)
  
  dat
}

prep_expr_matrix_gene_median <- function(gex_df, gene_id_col = NULL) {
  gex <- as.data.frame(gex_df, check.names = FALSE)
  
  if (ncol(gex) < 2) {
    stop("Expression data must contain at least one gene column and sample columns.")
  }
  
  if (is.null(gene_id_col)) {
    if ("geneID" %in% names(gex)) {
      gene_id_col <- "geneID"
    } else {
      non_numeric_cols <- names(gex)[
        vapply(gex, function(x) !is.numeric(x), logical(1))
      ]
      
      if (length(non_numeric_cols) >= 1) {
        gene_id_col <- non_numeric_cols[1]
      } else {
        gene_id_col <- names(gex)[1]
      }
    }
  }
  
  if (!gene_id_col %in% names(gex)) {
    stop("`gene_id_col` was not found in expression data: ", gene_id_col)
  }
  
  gene_raw <- norm_chr(gex[[gene_id_col]])
  keep_gene <- !is.na(gene_raw) & gene_raw != ""
  gex <- gex[keep_gene, , drop = FALSE]
  gene_raw <- gene_raw[keep_gene]
  
  sample_cols <- setdiff(names(gex), gene_id_col)
  if (length(sample_cols) < 3) {
    stop("Too few sample columns after excluding the gene column.")
  }
  
  gex_num <- gex[, sample_cols, drop = FALSE]
  gex_num[] <- lapply(gex_num, safe_numeric)
  
  num_cols <- names(gex_num)[vapply(gex_num, is.numeric, logical(1))]
  if (length(num_cols) < 3) {
    stop("Too few numeric sample columns after excluding the gene column.")
  }
  
  gex_num <- gex_num[, num_cols, drop = FALSE]
  gex_num$gene_id <- gene_raw
  
  gex_med <- gex_num |>
    group_by(gene_id) |>
    summarise(
      across(everything(), ~ stats::median(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  mat <- as.matrix(gex_med[, setdiff(names(gex_med), "gene_id"), drop = FALSE])
  rownames(mat) <- gex_med$gene_id
  storage.mode(mat) <- "numeric"
  
  if (anyDuplicated(rownames(mat))) {
    stop("Duplicated gene IDs remained after median collapsing.")
  }
  
  mat
}

run_gsva_wide <- function(expr_mat, signatures, method = "gsva", prefix = "GSVA_") {
  stopifnot(is.list(signatures), length(signatures) >= 1)
  stopifnot(!is.null(names(signatures)), all(names(signatures) != ""))
  
  expr_mat <- ensure_matrix_dimnames(expr_mat, cohort_label = "GSVA_input")
  
  row_sd <- apply(expr_mat, 1, stats::sd, na.rm = TRUE)
  expr_mat <- expr_mat[!is.na(row_sd) & row_sd > 0, , drop = FALSE]
  expr_mat <- ensure_matrix_dimnames(expr_mat, cohort_label = "GSVA_input_after_filter")
  
  if (nrow(expr_mat) < 2) {
    stop("Less than two variable genes remain after filtering constant genes.")
  }
  
  gset <- lapply(signatures, function(x) {
    genes <- clean_geneset(x)
    intersect(genes, rownames(expr_mat))
  })
  
  gset_size <- vapply(gset, length, integer(1))
  if (any(gset_size < 1)) {
    stop(
      "At least one signature has no overlapping genes in the expression matrix: ",
      paste(names(gset)[gset_size < 1], collapse = ", ")
    )
  }
  
  gs <- GSVA::gsva(expr_mat, gset, method = method)
  
  as.data.frame(t(gs), check.names = FALSE) |>
    rownames_to_column("sample") |>
    rename_with(~ paste0(prefix, .x), -sample)
}

fit_continuous_cox <- function(data, time_col, event_col, score_col, cohort, subgroup = "All") {
  dat <- data |>
    transmute(
      time = safe_numeric(.data[[time_col]]),
      event = safe_numeric(.data[[event_col]]),
      score = safe_numeric(.data[[score_col]])
    ) |>
    filter(!is.na(time), !is.na(event), !is.na(score), event %in% c(0, 1))
  
  n_total <- nrow(dat)
  n_event <- sum(dat$event == 1, na.rm = TRUE)
  
  if (n_total < 5 || stats::sd(dat$score, na.rm = TRUE) == 0) {
    return(tibble(
      cohort = cohort,
      subgroup = subgroup,
      n = n_total,
      n_event = n_event,
      hr = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      beta = NA_real_,
      se = NA_real_
    ))
  }
  
  fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ score, data = dat),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble(
      cohort = cohort,
      subgroup = subgroup,
      n = n_total,
      n_event = n_event,
      hr = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      beta = NA_real_,
      se = NA_real_
    ))
  }
  
  sm <- summary(fit)
  
  tibble(
    cohort = cohort,
    subgroup = subgroup,
    n = n_total,
    n_event = n_event,
    hr = unname(sm$coefficients[1, "exp(coef)"]),
    conf_low = unname(sm$conf.int[1, "lower .95"]),
    conf_high = unname(sm$conf.int[1, "upper .95"]),
    p_value = unname(sm$coefficients[1, "Pr(>|z|)"]),
    beta = unname(sm$coefficients[1, "coef"]),
    se = unname(sm$coefficients[1, "se(coef)"])
  )
}

make_forest_plot <- function(cox_tbl, file, title_text, endpoint_filter) {
  plot_df <- cox_tbl |>
    filter(endpoint == endpoint_filter) |>
    mutate(
      display_label = case_when(
        cohort == "JAVELIN" ~ paste0("JAVELIN: ", subgroup),
        cohort == "UC_Genome" & subgroup == "Immunotherapy_Yes" ~
          "UC-GENOME: Immunotherapy yes",
        cohort == "UC_Genome" & subgroup == "Chemotherapy_Yes_Immunotherapy_NotYes" ~
          "UC-GENOME: Chemotherapy yes / Immunotherapy not yes",
        cohort == "TCGA_KIRC" ~ "TCGA-KIRC",
        cohort == "UNC108" ~ "UNC108",
        cohort == "HUGABOOM" & subgroup == "PFS" ~ "HUGABOOM: PFS",
        cohort == "HUGABOOM" & subgroup == "OS" ~ "HUGABOOM: OS",
        TRUE ~ paste(cohort, subgroup, sep = " | ")
      ),
      hr_ci_text = ifelse(
        is.na(hr),
        "HR=NA",
        paste0(
          "HR ", sprintf("%.2f", hr),
          " (", sprintf("%.2f", conf_low), "-", sprintf("%.2f", conf_high), ")"
        )
      ),
      p_text = ifelse(
        is.na(p_value),
        "p=NA",
        paste0("p=", format.pval(p_value, digits = 2, eps = 1e-4))
      ),
      n_text = paste0("n=", n, ", events=", n_event),
      label_text = paste(hr_ci_text, p_text, n_text, sep = " | "),
      display_label = forcats::fct_rev(factor(display_label, levels = display_label))
    )
  
  x_min <- min(plot_df$conf_low, na.rm = TRUE)
  x_max <- max(plot_df$conf_high, na.rm = TRUE)
  text_x <- x_max * 1.8
  
  p <- ggplot(plot_df, aes(x = hr, y = display_label)) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
    geom_errorbarh(
      aes(xmin = conf_low, xmax = conf_high),
      height = 0.18,
      na.rm = TRUE,
      color = "black"
    ) +
    geom_point(size = 2.8, na.rm = TRUE) +
    geom_text(
      aes(x = text_x, label = label_text),
      hjust = 0,
      size = 3.2,
      na.rm = TRUE
    ) +
    scale_x_log10(limits = c(x_min * 0.8, text_x * 1.7)) +
    labs(
      title = title_text,
      x = "Hazard ratio per unit increase in ATF3 GSVA score",
      y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9),
      plot.margin = margin(5.5, 220, 5.5, 5.5)
    )
  
  ggsave(file, p, width = 13, height = max(4.5, 0.75 * nrow(plot_df)))
}
# ===============================================================
# 2) JAVELIN: PFS continuous Cox, arm-specific only
# ===============================================================

load_javelin_data <- function() {
  gex <- readxl::read_excel(
    "JAVELIN_Renal_101/41591_2020_1044_MOESM3_ESM.xlsx",
    sheet = "S13_Gene_expression_TPM",
    skip = 1
  )
  
  clin <- readxl::read_excel(
    "JAVELIN_Renal_101/41591_2020_1044_MOESM3_ESM.xlsx",
    sheet = "S11_Clinical_data",
    skip = 1
  )
  
  list(gex = as.data.frame(gex), clin = as.data.frame(clin))
}
analyze_javelin <- function(
    signatures,
    gsva_method = "gsva",
    use_cache = TRUE,
    force_recompute = FALSE
) {
  dat <- read_or_build_gsva_data(
    cache_file = cache_path("JAVELIN"),
    use_cache = use_cache,
    force_recompute = force_recompute,
    builder_fun = function() {
      dat0 <- load_javelin_data()
      expr_mat <- prep_expr_matrix_gene_median(dat0$gex)
      
      clin <- dat0$clin |>
        transmute(
          sample = norm_chr(ID),
          arm = norm_chr(TRT01P),
          PFS_time = safe_numeric(PFS_P),
          cens = safe_numeric(PFS_P_CNSR)
        ) |>
        mutate(
          PFS_event = case_when(
            cens %in% c(0, 1) ~ 1 - cens,
            TRUE ~ NA_real_
          )
        ) |>
        filter(!is.na(sample), sample != "", !is.na(arm), arm != "")
      
      common <- intersect(colnames(expr_mat), clin$sample)
      
      clin2 <- clin |>
        filter(sample %in% common) |>
        distinct(sample, .keep_all = TRUE) |>
        arrange(match(sample, common))
      
      expr2 <- expr_mat[, clin2$sample, drop = FALSE]
      expr2 <- ensure_matrix_dimnames(expr2, cohort_label = "JAVELIN")
      
      scores <- run_gsva_wide(expr2, signatures, method = gsva_method, prefix = "GSVA_")
      
      clin2 |>
        left_join(scores, by = "sample")
    }
  )
  
  arm_list <- split(dat, dat$arm)
  
  arm_res <- lapply(arm_list, function(df_arm) {
    fit_continuous_cox(
      df_arm,
      time_col = "PFS_time",
      event_col = "PFS_event",
      score_col = "GSVA_ATF3_final",
      cohort = "JAVELIN",
      subgroup = unique(df_arm$arm)[1]
    )
  }) |>
    bind_rows()
  
  arm_res
}
# ===============================================================
# 3) UC_Genome: ICI and chemo/no-ICI
# ===============================================================

prep_expr_matrix_uc <- function(expr_df, gene_col_guess = 1) {
  expr_df <- as.data.frame(expr_df, check.names = FALSE)
  
  cn <- colnames(expr_df)
  if (is.null(cn)) {
    cn <- rep("", ncol(expr_df))
  }
  if (any(!nzchar(cn))) {
    idx <- which(!nzchar(cn))
    cn[idx] <- paste0("V", idx)
    colnames(expr_df) <- cn
  }
  
  gene_col <- colnames(expr_df)[gene_col_guess]
  genes <- as.character(expr_df[[gene_col]])
  
  expr_num <- expr_df[, setdiff(colnames(expr_df), gene_col), drop = FALSE]
  expr_num[] <- lapply(expr_num, safe_numeric)
  
  num_cols <- colnames(expr_num)[vapply(expr_num, is.numeric, logical(1))]
  if (length(num_cols) < 3) {
    stop("Too few numeric expression columns found in UC_Genome.")
  }
  
  expr_num <- expr_num[, num_cols, drop = FALSE]
  row_mean <- rowMeans(expr_num, na.rm = TRUE)
  
  tmp <- expr_num
  tmp$gene <- genes
  tmp$row_mean <- row_mean
  
  tmp2 <- tmp |>
    arrange(desc(row_mean)) |>
    distinct(gene, .keep_all = TRUE)
  
  mat <- tmp2 |>
    select(-gene, -row_mean) |>
    as.matrix()
  
  rownames(mat) <- tmp2$gene
  storage.mode(mat) <- "numeric"
  
  new_ids <- extract_uc_sample4(colnames(mat))
  keep <- !is.na(new_ids)
  mat <- mat[, keep, drop = FALSE]
  colnames(mat) <- new_ids[keep]
  
  if (anyDuplicated(colnames(mat))) {
    stop("Duplicated UC_Genome sample IDs after extracting 217-XXXX.")
  }
  
  mat
}

prep_mapping_from_key <- function(key_df) {
  key_df |>
    as_tibble() |>
    transmute(
      sample = extract_uc_sample4(SEQUENCE_NO),
      patient_id = as.character(`RndSubNo...1`)
    ) |>
    filter(!is.na(sample), sample != "", !is.na(patient_id), patient_id != "") |>
    distinct(sample, .keep_all = TRUE)
}

prep_uc_meta <- function(meta_df) {
  meta_df |>
    slice(-(1:4)) |>
    transmute(
      patient_id = as.character(`#Patient Identifier`),
      survival_status = norm_chr(`Survival Status`),
      OS_time = safe_numeric(`Survival Time`),
      immunotherapy = norm_chr(Immunotherapy),
      chemotherapy = norm_chr(Chemotherapy)
    ) |>
    mutate(
      OS_event = case_when(
        str_to_upper(survival_status) %in% c("DECEASED", "DEAD", "DIED", "1") ~ 1,
        str_to_upper(survival_status) %in% c("LIVING", "ALIVE", "0") ~ 0,
        TRUE ~ NA_real_
      ),
      ici_yes = is_yes(immunotherapy),
      chemo_yes = is_yes(chemotherapy)
    )
}

analyze_uc_genome <- function(
    signatures,
    gsva_method = "gsva",
    use_cache = TRUE,
    force_recompute = FALSE
) {
  dat <- read_or_build_gsva_data(
    cache_file = cache_path("UC_Genome"),
    use_cache = use_cache,
    force_recompute = force_recompute,
    builder_fun = function() {
      gex <- read.table(
        "UC_Genome_UNC/salmon_tpm_gene.matrix.tsv",
        header = TRUE,
        sep = "\t",
        check.names = FALSE,
        quote = "",
        comment.char = ""
      )
      
      meta <- readxl::read_excel("UC_Genome_UNC/cBioPortal_Upload_26MAR2020.xlsx")
      key <- read.table(
        "UC_Genome_UNC/BCAN_ID_KEY.txt",
        header = TRUE,
        sep = "\t",
        check.names = FALSE,
        quote = "",
        comment.char = ""
      )
      
      expr_mat <- prep_expr_matrix_uc(gex)
      mapping <- prep_mapping_from_key(key)
      meta2 <- prep_uc_meta(meta)
      
      sample_meta <- mapping |>
        inner_join(meta2, by = "patient_id") |>
        filter(!is.na(OS_time), OS_event %in% c(0, 1))
      
      common <- intersect(colnames(expr_mat), sample_meta$sample)
      
      sample_meta <- sample_meta |>
        filter(sample %in% common) |>
        distinct(sample, .keep_all = TRUE) |>
        arrange(match(sample, common))
      
      expr2 <- expr_mat[, sample_meta$sample, drop = FALSE]
      expr2 <- ensure_matrix_dimnames(expr2, cohort_label = "UC_Genome")
      
      scores <- run_gsva_wide(expr2, signatures, method = gsva_method, prefix = "GSVA_")
      
      sample_meta |>
        left_join(scores, by = "sample")
    }
  )
  
  res_ici <- dat |>
    filter(ici_yes) |>
    fit_continuous_cox(
      time_col = "OS_time",
      event_col = "OS_event",
      score_col = "GSVA_ATF3_final",
      cohort = "UC_Genome",
      subgroup = "Immunotherapy_Yes"
    )
  
  res_chemo <- dat |>
    filter(chemo_yes, !ici_yes | is.na(ici_yes)) |>
    fit_continuous_cox(
      time_col = "OS_time",
      event_col = "OS_event",
      score_col = "GSVA_ATF3_final",
      cohort = "UC_Genome",
      subgroup = "Chemotherapy_Yes_Immunotherapy_NotYes"
    )
  
  bind_rows(res_ici, res_chemo)
}

# ===============================================================
# 4) TCGA-KIRC: OS continuous Cox
# ===============================================================

prep_tcga_os <- function(meta_df) {
  meta_df <- as.data.frame(meta_df, check.names = FALSE, stringsAsFactors = FALSE)
  meta_df[meta_df == "--"] <- NA
  
  meta_df$demographic.days_to_death <- safe_numeric(meta_df$demographic.days_to_death)
  meta_df$diagnoses.days_to_last_follow_up <- safe_numeric(meta_df$diagnoses.days_to_last_follow_up)
  
  max_or_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA_real_)
    }
    max(x)
  }
  
  first_non_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA_character_)
    }
    as.character(x[1])
  }
  
  meta_df |>
    group_by(cases.submitter_id) |>
    summarise(
      vital_status = first_non_na(demographic.vital_status),
      days_to_death = max_or_na(demographic.days_to_death),
      days_to_last_follow_up = max_or_na(diagnoses.days_to_last_follow_up),
      .groups = "drop"
    ) |>
    mutate(
      OS_event = ifelse(vital_status == "Dead", 1, 0),
      OS_time = ifelse(OS_event == 1, days_to_death, days_to_last_follow_up)
    ) |>
    filter(!is.na(cases.submitter_id), !is.na(OS_time), OS_event %in% c(0, 1))
}

analyze_tcga_kirc <- function(
    signatures,
    gsva_method = "gsva",
    use_cache = TRUE,
    force_recompute = FALSE
) {
  dat <- read_or_build_gsva_data(
    cache_file = cache_path("TCGA_KIRC"),
    use_cache = use_cache,
    force_recompute = force_recompute,
    builder_fun = function() {
      gex <- read.table(
        "TCGA-KIRC/salmon_tpm_gene.matrix.tsv",
        header = TRUE,
        sep = "\t",
        check.names = FALSE,
        quote = "",
        comment.char = "",
        stringsAsFactors = FALSE
      )
      
      meta <- read.table(
        "TCGA-KIRC/meta/clinical.tsv",
        header = TRUE,
        sep = "\t",
        check.names = FALSE,
        quote = "",
        comment.char = "",
        stringsAsFactors = FALSE
      )
      
      gex <- as.data.frame(gex, check.names = FALSE)
      rownames(gex) <- gex[[1]]
      gex <- gex[, -1, drop = FALSE]
      gex[] <- lapply(gex, safe_numeric)
      
      expr_mat <- as.matrix(gex)
      storage.mode(expr_mat) <- "numeric"
      
      gex_cols <- colnames(expr_mat)
      gex_info <- tibble(
        sample = gex_cols,
        patient_id = substr(gex_cols, 1, 12),
        sample_type = vapply(
          strsplit(gex_cols, "-"),
          function(x) substr(x[4], 1, 2),
          character(1)
        )
      ) |>
        filter(sample_type == "01") |>
        distinct(patient_id, .keep_all = TRUE)
      
      os_tbl <- prep_tcga_os(meta)
      common <- intersect(gex_info$patient_id, os_tbl$cases.submitter_id)
      
      gex_info2 <- gex_info |>
        filter(patient_id %in% common) |>
        arrange(match(patient_id, common))
      
      os_tbl2 <- os_tbl |>
        filter(cases.submitter_id %in% common) |>
        arrange(match(cases.submitter_id, common))
      
      expr2 <- expr_mat[, gex_info2$sample, drop = FALSE]
      colnames(expr2) <- clean_dimnames(gex_info2$patient_id)
      expr2 <- ensure_matrix_dimnames(expr2, cohort_label = "TCGA_KIRC")
      
      row_var <- apply(expr2, 1, stats::sd, na.rm = TRUE)
      expr2 <- expr2[!is.na(row_var) & row_var > 0, , drop = FALSE]
      
      scores <- run_gsva_wide(expr2, signatures, method = gsva_method, prefix = "GSVA_")
      
      os_tbl2 |>
        transmute(
          sample = cases.submitter_id,
          OS_time,
          OS_event
        ) |>
        left_join(scores, by = "sample")
    }
  )
  
  fit_continuous_cox(
    dat,
    time_col = "OS_time",
    event_col = "OS_event",
    score_col = "GSVA_ATF3_final",
    cohort = "TCGA_KIRC",
    subgroup = "All"
  )
}

# ===============================================================
# 5) UNC108: OS continuous Cox
# ===============================================================

prep_expr_matrix_unc108 <- function(gex_df) {
  gex_df <- as.data.frame(gex_df, check.names = FALSE)
  rownames(gex_df) <- gex_df[[1]]
  gex_df <- gex_df[, -1, drop = FALSE]
  gex_df[] <- lapply(gex_df, safe_numeric)
  mat <- as.matrix(gex_df)
  storage.mode(mat) <- "numeric"
  mat
}

prep_unc108_meta <- function(meta_df, key_df) {
  key2 <- key_df |>
    transmute(
      sample_baci = `Sample ID`,
      sample_rs = normalize_unc108_rs(`Omniseq_RS_ID (RNAseq)`)
    ) |>
    filter(!is.na(sample_baci), !is.na(sample_rs), sample_baci != "", sample_rs != "")
  
  meta_df |>
    transmute(
      sample_baci = `RNA Data Sample name`,
      OS_time = safe_numeric(`Overall Survival`),
      OS_event = case_when(
        Alive %in% c("No", "Dead", "Deceased") ~ 1,
        Alive %in% c("Yes", "Alive") ~ 0,
        TRUE ~ NA_real_
      )
    ) |>
    left_join(key2, by = "sample_baci") |>
    filter(!is.na(sample_rs), !is.na(OS_time), OS_event %in% c(0, 1))
}

analyze_unc108 <- function(rds_file) {
  dat <- readRDS(rds_file)
  
  check_required_columns(
    dat,
    c("GSVA_ATF3_final"),
    "UNC108 RDS"
  )
  
  res <- list()
  
  if (has_required_columns(dat, c("OS_time", "OS_event"))) {
    res <- append(
      res,
      list(
        run_one_endpoint(
          dat = dat,
          cohort = "UNC108",
          subgroup = "All",
          endpoint = "OS",
          time_col = "OS_time",
          event_col = "OS_event"
        )
      )
    )
  } else {
    message("UNC108: OS_time / OS_event not found, skipping OS.")
  }
  
  if (has_required_columns(dat, c("PFS_time", "PFS_event"))) {
    res <- append(
      res,
      list(
        run_one_endpoint(
          dat = dat,
          cohort = "UNC108",
          subgroup = "All",
          endpoint = "PFS",
          time_col = "PFS_time",
          event_col = "PFS_event"
        )
      )
    )
  } else {
    message("UNC108: PFS_time / PFS_event not found, skipping PFS.")
  }
  
  if (length(res) == 0) {
    stop("UNC108 RDS does not contain usable OS or PFS columns.")
  }
  
  bind_rows(res)
}

# ===============================================================
# 7) HUGABOOM: PFS continuous Cox
# ===============================================================

prep_expr_matrix_hugaboom <- function(expr_df, gene_col = "Description") {
  expr_df <- as.data.frame(expr_df, check.names = FALSE)
  
  sample_cols <- setdiff(colnames(expr_df), gene_col)
  expr_df[sample_cols] <- lapply(expr_df[sample_cols], safe_numeric)
  num_cols <- sample_cols[vapply(expr_df[sample_cols], is.numeric, logical(1))]
  expr_df$row_mean <- rowMeans(expr_df[, num_cols, drop = FALSE], na.rm = TRUE)
  
  expr_unique <- expr_df |>
    arrange(desc(row_mean)) |>
    distinct(.data[[gene_col]], .keep_all = TRUE)
  
  rownames(expr_unique) <- expr_unique[[gene_col]]
  mat <- as.matrix(expr_unique[, num_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat
}

prep_hugaboom_meta <- function(meta_df) {
  meta_df |>
    transmute(
      sample = SUBJECT_ID,
      PFS_time = safe_numeric(PFS),
      PFS_event = safe_numeric(PFS_CNSR),
      OS_time = safe_numeric(OS),
      OS_event = safe_numeric(OS_CNSR)
    ) |>
    filter(
      !is.na(sample),
      sample != "",
      !is.na(PFS_time),
      PFS_event %in% c(0, 1)
    )
}

analyze_hugaboom <- function(rds_file) {
  dat <- readRDS(rds_file)
  
  check_required_columns(
    dat,
    c("GSVA_ATF3_final"),
    "HUGABOOM RDS"
  )
  
  res <- list()
  
  if (has_required_columns(dat, c("PFS_time", "PFS_event"))) {
    res <- append(
      res,
      list(
        run_one_endpoint(
          dat = dat,
          cohort = "HUGABOOM",
          subgroup = "All",
          endpoint = "PFS",
          time_col = "PFS_time",
          event_col = "PFS_event"
        )
      )
    )
  } else {
    message("HUGABOOM: PFS_time / PFS_event not found, skipping PFS.")
  }
  
  if (has_required_columns(dat, c("OS_time", "OS_event"))) {
    res <- append(
      res,
      list(
        run_one_endpoint(
          dat = dat,
          cohort = "HUGABOOM",
          subgroup = "All",
          endpoint = "OS",
          time_col = "OS_time",
          event_col = "OS_event"
        )
      )
    )
  } else {
    message("HUGABOOM: OS_time / OS_event not found, skipping OS.")
  }
  
  if (length(res) == 0) {
    stop("HUGABOOM RDS does not contain usable OS or PFS columns.")
  }
  
  bind_rows(res)
}
# ===============================================================
# 8) Main runner
# ===============================================================

run_main_atf3_continuous_cox <- function(
    signatures,
    out_dir,
    gsva_method = "gsva",
    use_cache = TRUE,
    force_recompute = FALSE
) {
  res_list <- list(
    analyze_javelin(
      signatures,
      gsva_method = gsva_method,
      use_cache = use_cache,
      force_recompute = force_recompute
    ),
    analyze_uc_genome(
      signatures,
      gsva_method = gsva_method,
      use_cache = use_cache,
      force_recompute = force_recompute
    ),
    analyze_tcga_kirc(
      signatures,
      gsva_method = gsva_method,
      use_cache = use_cache,
      force_recompute = force_recompute
    ),
    analyze_unc108(
      signatures,
      gsva_method = gsva_method,
      use_cache = use_cache,
      force_recompute = force_recompute
    ),
    analyze_hugaboom(
      signatures,
      gsva_method = gsva_method,
      use_cache = use_cache,
      force_recompute = force_recompute
    )
  )
  
  cox_tbl <- bind_rows(res_list) |>
    mutate(
      endpoint = case_when(
        cohort == "JAVELIN" ~ "PFS",
        cohort == "HUGABOOM" ~ "PFS",
        TRUE ~ "OS"
      )
    ) |>
    relocate(cohort, subgroup, endpoint)
  
  readr::write_csv(
    cox_tbl,
    file.path(out_dir, "atf3_continuous_cox_results.csv")
  )
  
  saveRDS(
    cox_tbl,
    file.path(out_dir, "atf3_continuous_cox_results.rds")
  )
  
  make_forest_plot(
    cox_tbl,
    file = file.path(out_dir, "atf3_continuous_cox_forest_pfs.pdf"),
    title_text = "ATF3 continuous Cox across PFS cohorts",
    endpoint_filter = "PFS"
  )
  
  make_forest_plot(
    cox_tbl,
    file = file.path(out_dir, "atf3_continuous_cox_forest_os.pdf"),
    title_text = "ATF3 continuous Cox across OS cohorts",
    endpoint_filter = "OS"
  )
  
  invisible(cox_tbl)
}

# ===============================================================
# 9) Run
# ===============================================================

cox_results <- run_main_atf3_continuous_cox(
  signatures = signatures,
  out_dir = out_dir,
  gsva_method = "gsva",
  use_cache = TRUE,
  force_recompute = FALSE
)

print(cox_results)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(survival)
  library(ggplot2)
  library(forcats)
  library(purrr)
})

# ===============================================================
# 0) Output directory
# ===============================================================

out_dir <- "./ATF3_continuous_cox_main_from_rds"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ===============================================================
# 1) Helpers
# ===============================================================

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

fit_continuous_cox <- function(
    data,
    time_col,
    event_col,
    score_col = "GSVA_ATF3_final",
    cohort,
    subgroup = "All",
    endpoint
) {
  dat <- data |>
    transmute(
      time = safe_numeric(.data[[time_col]]),
      event = safe_numeric(.data[[event_col]]),
      score = safe_numeric(.data[[score_col]])
    ) |>
    filter(!is.na(time), !is.na(event), !is.na(score), event %in% c(0, 1))
  
  n_total <- nrow(dat)
  n_event <- sum(dat$event == 1, na.rm = TRUE)
  
  if (n_total < 5 || is.na(stats::sd(dat$score, na.rm = TRUE)) || stats::sd(dat$score, na.rm = TRUE) == 0) {
    return(
      tibble(
        cohort = cohort,
        subgroup = subgroup,
        endpoint = endpoint,
        n = n_total,
        n_event = n_event,
        hr = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        p_value = NA_real_,
        beta = NA_real_,
        se = NA_real_
      )
    )
  }
  
  fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ score, data = dat),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(
      tibble(
        cohort = cohort,
        subgroup = subgroup,
        endpoint = endpoint,
        n = n_total,
        n_event = n_event,
        hr = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        p_value = NA_real_,
        beta = NA_real_,
        se = NA_real_
      )
    )
  }
  
  sm <- summary(fit)
  
  tibble(
    cohort = cohort,
    subgroup = subgroup,
    endpoint = endpoint,
    n = n_total,
    n_event = n_event,
    hr = unname(sm$coefficients[1, "exp(coef)"]),
    conf_low = unname(sm$conf.int[1, "lower .95"]),
    conf_high = unname(sm$conf.int[1, "upper .95"]),
    p_value = unname(sm$coefficients[1, "Pr(>|z|)"]),
    beta = unname(sm$coefficients[1, "coef"]),
    se = unname(sm$coefficients[1, "se(coef)"])
  )
}

check_required_columns <- function(dat, cols, object_name) {
  missing_cols <- setdiff(cols, names(dat))
  
  if (length(missing_cols) > 0) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  invisible(dat)
}

make_forest_plot <- function(cox_tbl, file, title_text, endpoint_filter) {
  plot_df <- cox_tbl |>
    filter(endpoint == endpoint_filter) |>
    mutate(
      display_label = case_when(
        cohort == "JAVELIN" ~ paste0("JAVELIN: ", subgroup),
        cohort == "UC_Genome" & subgroup == "Immunotherapy_Yes" ~
          "UC-GENOME: Immunotherapy yes",
        cohort == "UC_Genome" & subgroup == "Chemotherapy_Yes_Immunotherapy_NotYes" ~
          "UC-GENOME: Chemotherapy yes / Immunotherapy not yes",
        cohort == "TCGA_KIRC" ~ "TCGA-KIRC",
        cohort == "UNC108" & endpoint == "OS" ~ "UNC108: OS",
        cohort == "UNC108" & endpoint == "PFS" ~ "UNC108: PFS",
        cohort == "HUGABOOM" & endpoint == "OS" ~ "HUGABOOM: OS",
        cohort == "HUGABOOM" & endpoint == "PFS" ~ "HUGABOOM: PFS",
        TRUE ~ paste(cohort, subgroup, endpoint, sep = " | ")
      ),
      hr_ci_text = ifelse(
        is.na(hr),
        "HR=NA",
        paste0(
          "HR ", sprintf("%.2f", hr),
          " (", sprintf("%.2f", conf_low), "-", sprintf("%.2f", conf_high), ")"
        )
      ),
      p_text = ifelse(
        is.na(p_value),
        "p=NA",
        paste0("p=", format.pval(p_value, digits = 2, eps = 1e-4))
      ),
      n_text = paste0("n=", n, ", events=", n_event),
      label_text = paste(hr_ci_text, p_text, n_text, sep = " | "),
      display_label = forcats::fct_rev(factor(display_label, levels = display_label))
    )
  
  if (nrow(plot_df) == 0) {
    warning("No rows available for endpoint: ", endpoint_filter)
    return(invisible(NULL))
  }
  
  x_min <- min(plot_df$conf_low, na.rm = TRUE)
  x_max <- max(plot_df$conf_high, na.rm = TRUE)
  text_x <- x_max * 1.8
  
  p <- ggplot(plot_df, aes(x = hr, y = display_label)) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
    geom_errorbarh(
      aes(xmin = conf_low, xmax = conf_high),
      height = 0.18,
      na.rm = TRUE,
      color = "black"
    ) +
    geom_point(size = 2.8, na.rm = TRUE) +
    geom_text(
      aes(x = text_x, label = label_text),
      hjust = 0,
      size = 3.2,
      na.rm = TRUE
    ) +
    scale_x_log10(limits = c(x_min * 0.8, text_x * 1.7)) +
    labs(
      title = title_text,
      x = "Hazard ratio per unit increase in ATF3 GSVA score",
      y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9),
      plot.margin = margin(5.5, 220, 5.5, 5.5)
    )
  
  ggsave(file, p, width = 13, height = max(4.5, 0.75 * nrow(plot_df)))
}

run_one_endpoint <- function(dat, cohort, subgroup, endpoint, time_col, event_col) {
  fit_continuous_cox(
    data = dat,
    time_col = time_col,
    event_col = event_col,
    score_col = "GSVA_ATF3_final",
    cohort = cohort,
    subgroup = subgroup,
    endpoint = endpoint
  )
}

# ===============================================================
# 2) Cohort analyzers from precomputed RDS
# ===============================================================

analyze_javelin <- function(rds_file) {
  dat <- readRDS(rds_file)
  
  check_required_columns(
    dat,
    c("arm", "GSVA_ATF3_final", "PFS_time", "PFS_event"),
    "JAVELIN RDS"
  )
  
  split(dat, dat$arm) |>
    imap(\(df_arm, arm_name) {
      run_one_endpoint(
        dat = df_arm,
        cohort = "JAVELIN",
        subgroup = arm_name,
        endpoint = "PFS",
        time_col = "PFS_time",
        event_col = "PFS_event"
      )
    }) |>
    bind_rows()
}

analyze_uc_genome <- function(rds_file) {
  dat <- readRDS(rds_file)
  
  check_required_columns(
    dat,
    c(
      "GSVA_ATF3_final",
      "OS_time", "OS_event",
      "ici_yes", "chemo_yes"
    ),
    "UC_Genome RDS"
  )
  
  res_ici <- dat |>
    filter(ici_yes) |>
    run_one_endpoint(
      cohort = "UC_Genome",
      subgroup = "Immunotherapy_Yes",
      endpoint = "OS",
      time_col = "OS_time",
      event_col = "OS_event"
    )
  
  res_chemo <- dat |>
    filter(chemo_yes, !ici_yes | is.na(ici_yes)) |>
    run_one_endpoint(
      cohort = "UC_Genome",
      subgroup = "Chemotherapy_Yes_Immunotherapy_NotYes",
      endpoint = "OS",
      time_col = "OS_time",
      event_col = "OS_event"
    )
  
  bind_rows(res_ici, res_chemo)
}

analyze_tcga_kirc <- function(rds_file) {
  dat <- readRDS(rds_file)
  
  check_required_columns(
    dat,
    c("GSVA_ATF3_final", "OS_time", "OS_event"),
    "TCGA_KIRC RDS"
  )
  
  run_one_endpoint(
    dat = dat,
    cohort = "TCGA_KIRC",
    subgroup = "All",
    endpoint = "OS",
    time_col = "OS_time",
    event_col = "OS_event"
  )
}

analyze_unc108 <- function(rds_file) {
  dat <- readRDS(rds_file)
  
  check_required_columns(
    dat,
    c("GSVA_ATF3_final", "OS_time", "OS_event", "PFS_time", "PFS_event"),
    "UNC108 RDS"
  )
  
  bind_rows(
    run_one_endpoint(
      dat = dat,
      cohort = "UNC108",
      subgroup = "All",
      endpoint = "OS",
      time_col = "OS_time",
      event_col = "OS_event"
    ),
    run_one_endpoint(
      dat = dat,
      cohort = "UNC108",
      subgroup = "All",
      endpoint = "PFS",
      time_col = "PFS_time",
      event_col = "PFS_event"
    )
  )
}

analyze_hugaboom <- function(rds_file) {
  dat <- readRDS(rds_file)
  
  check_required_columns(
    dat,
    c("GSVA_ATF3_final", "OS_time", "OS_event", "PFS_time", "PFS_event"),
    "HUGABOOM RDS"
  )
  
  bind_rows(
    run_one_endpoint(
      dat = dat,
      cohort = "HUGABOOM",
      subgroup = "All",
      endpoint = "PFS",
      time_col = "PFS_time",
      event_col = "PFS_event"
    ),
    run_one_endpoint(
      dat = dat,
      cohort = "HUGABOOM",
      subgroup = "All",
      endpoint = "OS",
      time_col = "OS_time",
      event_col = "OS_event"
    )
  )
}

# ===============================================================
# 3) Main runner
# ===============================================================

run_main_atf3_continuous_cox_from_rds <- function(
    javelin_rds,
    uc_genome_rds,
    tcga_kirc_rds,
    unc108_rds,
    hugaboom_rds,
    out_dir
) {
  res_list <- list(
    analyze_javelin(javelin_rds),
    analyze_uc_genome(uc_genome_rds),
    analyze_tcga_kirc(tcga_kirc_rds),
    analyze_unc108(unc108_rds),
    analyze_hugaboom(hugaboom_rds)
  )
  
  cox_tbl <- bind_rows(res_list) |>
    relocate(cohort, subgroup, endpoint)
  
  readr::write_csv(
    cox_tbl,
    file.path(out_dir, "atf3_continuous_cox_results.csv")
  )
  
  saveRDS(
    cox_tbl,
    file.path(out_dir, "atf3_continuous_cox_results.rds")
  )
  
  make_forest_plot(
    cox_tbl,
    file = file.path(out_dir, "atf3_continuous_cox_forest_pfs.pdf"),
    title_text = "ATF3 continuous Cox across PFS cohorts",
    endpoint_filter = "PFS"
  )
  
  make_forest_plot(
    cox_tbl,
    file = file.path(out_dir, "atf3_continuous_cox_forest_os.pdf"),
    title_text = "ATF3 continuous Cox across OS cohorts",
    endpoint_filter = "OS"
  )
  
  cox_tbl
}

# ===============================================================
# 4) File paths
# ===============================================================
# Cache directory path for precomputed GSVA RDS files.
path = './ATF3_continuous_cox_main/gsva_cache/'

javelin_rds <- "./ATF3_continuous_cox_main/gsva_cache/gsva_JAVELIN.rds"
uc_genome_rds <- "./ATF3_continuous_cox_main/gsva_cache/gsva_UC_Genome.rds"
tcga_kirc_rds <- "./ATF3_continuous_cox_main/gsva_cache/gsva_TCGA_KIRC.rds"
unc108_rds <- "./ATF3_continuous_cox_main/gsva_cache/gsva_UNC108.rds"
hugaboom_rds <- "./ATF3_continuous_cox_main/gsva_cache/gsva_hugaboom.rds"

# ===============================================================
# 5) Run
# ===============================================================
# Check whether all requested columns are present in the input data frame.
has_required_columns <- function(dat, cols) {
  all(cols %in% names(dat))
}
cox_results <- run_main_atf3_continuous_cox_from_rds(
  javelin_rds = javelin_rds,
  uc_genome_rds = uc_genome_rds,
  tcga_kirc_rds = tcga_kirc_rds,
  unc108_rds = unc108_rds,
  hugaboom_rds = hugaboom_rds,
  out_dir = out_dir
)

print(cox_results)





library(dplyr)
library(forestplot)
library(grid)

# Rename subgroup labels for the final figure output.
cox_results$subgroup <- c("Avel+Axi",'Suni','ICI','Chemo_only','ALL','ALL','ALL','ALL','ALL')


# Use the combined result table directly for the final forest plots.
cox2 <- cox_results

# PFS
pfs <- cox2 %>%
  filter(endpoint == "PFS") %>%
  mutate(
    study = ifelse(subgroup == "All", cohort, paste0(cohort, " | ", subgroup)),
    events_n = paste0(n_event, "/", n),
    hr_ci = sprintf("%.2f (%.2f-%.2f)", hr, conf_low, conf_high),
    p_txt = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  )

pfs_text <- cbind(
  c("Cohort / subgroup", pfs$study),
  c("Events/N", pfs$events_n),
  c("", rep("", nrow(pfs))),
  c("HR (95% CI)", pfs$hr_ci),
  c("P value", pfs$p_txt)
)

grDevices::cairo_pdf("forest_PFS.pdf", width = 10, height = 4.8, family = "Helvetica")
grid.newpage()
forestplot(
  labeltext = pfs_text,
  mean = c(NA, pfs$hr),
  lower = c(NA, pfs$conf_low),
  upper = c(NA, pfs$conf_high),
  is.summary = c(TRUE, rep(FALSE, nrow(pfs))),
  graph.pos = 3,
  align = c("l", "c", "c", "c", "c"),
  zero = 1,
  xlog = TRUE,
  clip = c(0.1, 4),
  xticks = c(0.25, 0.5, 1, 2, 4),
  boxsize = 0.2,
  lwd.ci = 2,
  lwd.zero = 1.2,
  ci.vertices = TRUE,
  ci.vertices.height = 0.1,
  col = fpColors(box = "black", line = "black", zero = "gray50"),
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.9),
    ticks = gpar(cex = 0.85),
    xlab = gpar(cex = 0.9),
    title = gpar(cex = 1.1, fontface = "bold")
  ),
  hrzl_lines = list(
    "1" = gpar(lwd = 1.2),
    "2" = gpar(lwd = 0.8, col = "gray70")
  ),
  colgap = unit(5, "mm"),
  lineheight = unit(8, "mm"),
  xlab = "Hazard ratio",
  title = "PFS"
)
dev.off()

# OS
os <- cox2 %>%
  filter(endpoint == "OS") %>%
  mutate(
    study = ifelse(subgroup == "All", cohort, paste0(cohort, " | ", subgroup)),
    events_n = paste0(n_event, "/", n),
    hr_ci = sprintf("%.2f (%.2f-%.2f)", hr, conf_low, conf_high),
    p_txt = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  )

os_text <- cbind(
  c("Cohort / subgroup", os$study),
  c("Events/N", os$events_n),
  c("", rep("", nrow(os))),
  c("HR (95% CI)", os$hr_ci),
  c("P value", os$p_txt)
)

grDevices::cairo_pdf("forest_OS.pdf", width = 10, height = 5.2, family = "Helvetica")
grid.newpage()
forestplot(
  labeltext = os_text,
  mean = c(NA, os$hr),
  lower = c(NA, os$conf_low),
  upper = c(NA, os$conf_high),
  is.summary = c(TRUE, rep(FALSE, nrow(os))),
  graph.pos = 3,
  align = c("l", "c", "c", "c", "c"),
  zero = 1,
  xlog = TRUE,
  clip = c(0.1, 4),
  xticks = c(0.25, 0.5, 1, 2, 4),
  boxsize = 0.2,
  lwd.ci = 2,
  lwd.zero = 1.2,
  ci.vertices = TRUE,
  ci.vertices.height = 0.1,
  col = fpColors(box = "black", line = "black", zero = "gray50"),
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.9),
    ticks = gpar(cex = 0.85),
    xlab = gpar(cex = 0.9),
    title = gpar(cex = 1.1, fontface = "bold")
  ),
  hrzl_lines = list(
    "1" = gpar(lwd = 1.2),
    "2" = gpar(lwd = 0.8, col = "gray70")
  ),
  colgap = unit(5, "mm"),
  lineheight = unit(8, "mm"),
  xlab = "Hazard ratio",
  title = "OS"
)
dev.off()



