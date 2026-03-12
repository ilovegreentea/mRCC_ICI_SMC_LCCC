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

# ---- 패널별 y축 정렬 ----
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

# DEG가 1개라도 걸리는 TF만 보고 싶으면
res_overlap2 <- subset(res_overlap, n_overlap > 0)

# overlap 많은 순으로
res_overlap2 <- res_overlap2[order(res_overlap2$n_overlap, decreasing = TRUE), ]

res_overlap2 %>% head()
 



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

# TF x 1 matrix (heatmap용)
jac_mat <- matrix(jac_vec, ncol = 1)
rownames(jac_mat) <- tf_names
colnames(jac_mat) <- "T0_DE"

# 보기 좋게: 상위 TF만 (예: top 50)
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
  name = "Jaccard index",                         # <- colorbar 라벨
  column_title = "Jaccard index (DEG vs regulon)", # <- 타이틀
  row_title = "TF regulon",
  cluster_rows = TRUE,
  cluster_columns = FALSE,                         # 여러 컬럼이면 TRUE로
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



