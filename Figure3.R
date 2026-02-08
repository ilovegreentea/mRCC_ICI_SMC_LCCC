# Cell type specific signature from : ref :: Bi, K., Kashima, S., Camp, S. Y., Meli, K., Saad, E., Titchen, B. M., ... & Van Allen, E. M. (2025). Myeloid cells mediate interferon-driven resistance to immunotherapy in advanced renal cell carcinoma. Immunity, 58(11), 2814-2829.
# code ref https://github.com/kevinbi2599/ccRCC_IFN_Resistance_Immunity2025

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

latent_vars <- c(patient_col, intersect(c("Timepoint"), colnames(seu@meta.data)))
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




