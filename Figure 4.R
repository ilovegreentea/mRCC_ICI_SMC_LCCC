#### Clonal bubble plot and hyper-extended standard :: ref :: Kjær, A., Kristjánsdóttir, N., Juul, R. I., Nordentoft, I., Birkenkamp-Demtröder, K., Ahrenfeldt, J., ... & Dyrskjøt, L. (2025). Low T cell diversity associates with poor outcome in bladder cancer: A comprehensive longitudinal analysis of the T cell receptor repertoire. Cell Reports Medicine, 6(5).
#### https://zenodo.org/records/14823578
####

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
})

## =========================
## Settings
## =========================
N_TOP <- 10                 # Responders with T0: Top N pre-existing + expanded
N_TOP_MISSING_T0 <- 10      # Responders without T0: Top N largest clones at T2
MIN_TCR_CELLS <- 500        # Minimum number of TCR+ cells per sample_name (= sample)
TIME_T0 <- "T0"
TIME_T2 <- "T2"
OUTDIR <- "plots/tumor_only"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## Set TRUE to avoid Unicode warnings (αβΔ–)
ASCII_HEADERS <- TRUE

## =========================
## 1) Prepare meta
## =========================
meta <- T_cell@meta.data %>% rownames_to_column("barcode")

clone_id_col <- if ("CTnt" %in% colnames(meta)) "CTnt" else "CTstrict"
stopifnot(clone_id_col %in% colnames(meta))
stopifnot(all(c("sample_name","patients","Timepoint","Response") %in% colnames(meta)))

meta_tcr %>% summarise(
  n_cells = n(),
  n_CTstrict = sum(!is.na(CTstrict) & CTstrict != ""),
  n_CTnt     = sum(!is.na(CTnt) & CTnt != ""),
  n_CTaa     = sum(!is.na(CTaa) & CTaa != ""),
  n_CTgene   = sum(!is.na(CTgene) & CTgene != "")
)

has_CTaa   <- "CTaa"   %in% colnames(meta)
has_CTgene <- "CTgene" %in% colnames(meta)
has_TNK    <- "T_NK_cell" %in% colnames(meta)

if (!has_TNK) stop("T_NK_cell column is missing in meta.data.")
if (!has_CTaa) message("Warning: CTaa is missing; CDR3α/β will be shown as '(NA)'.")
if (!has_CTgene) message("Warning: CTgene is missing; TRAV/J and TRBV/J will be shown as '(NA)'.")

## =========================
## 2) Keep TCR+ cells and define sample_id (sample_name)
##    (No tumor filter is applied here since everything is tumor)
## =========================
meta_tcr <- meta %>%
  filter(!is.na(.data[[clone_id_col]]), .data[[clone_id_col]] != "") %>%
  filter(!is.na(sample_name), sample_name != "") %>%
  filter(Timepoint %in% c(TIME_T0, TIME_T2)) %>%
  mutate(sample_id = sample_name)

## =========================
## 3) Apply sample-level TCR+ cell cutoff
## =========================
keep_samples <- meta_tcr %>%
  count(sample_id, name = "n_tcr_cells") %>%
  filter(n_tcr_cells >= MIN_TCR_CELLS) %>%
  pull(sample_id)

meta_tcr <- meta_tcr %>% filter(sample_id %in% keep_samples)

## =========================
## 4) Compute within-sample clonotype frequency (paper-style)
##    n_cells_per_clone / total_TCR_cells_in_sample
## =========================
clone_sample <- meta_tcr %>%
  group_by(sample_id, patients, Response, Timepoint,
           clone_id = .data[[clone_id_col]]) %>%
  summarise(
    n_cells = n(),
    T_NK_cell_mode = names(sort(table(T_NK_cell), decreasing = TRUE))[1],
    CTaa_mode = if (has_CTaa) names(sort(table(CTaa), decreasing = TRUE))[1] else NA_character_,
    CTgene_mode = if (has_CTgene) names(sort(table(CTgene), decreasing = TRUE))[1] else NA_character_,
    .groups = "drop"
  ) %>%
  group_by(sample_id) %>%
  mutate(freq = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  mutate(expanded = ifelse(freq > 0.002, "Hyper-Expanded", "Non-Expanded"))  # 0.2%

## =========================
## 5) Check whether each patient has T0/T2
## =========================
pt_tp <- meta_tcr %>%
  distinct(patients, Response, Timepoint) %>%
  mutate(has_tp = TRUE) %>%
  pivot_wider(names_from = Timepoint, values_from = has_tp, values_fill = FALSE) %>%
  mutate(
    T0_available = .data[[TIME_T0]],
    T2_available = .data[[TIME_T2]]
  )

## =========================
## 6) Patient-level wide table (fill missing timepoints with 0)
## =========================
cl_wide <- clone_sample %>%
  select(patients, Response, Timepoint, clone_id, freq, n_cells,
         T_NK_cell_mode, CTaa_mode, CTgene_mode) %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = c(freq, n_cells),
    values_fill = 0
  ) %>%
  left_join(pt_tp %>% select(patients, Response, T0_available, T2_available),
            by = c("patients", "Response"))

freq_t0_col <- paste0("freq_", TIME_T0)
freq_t2_col <- paste0("freq_", TIME_T2)
stopifnot(all(c(freq_t0_col, freq_t2_col) %in% colnames(cl_wide)))

cl_wide2 <- cl_wide %>%
  mutate(
    delta  = .data[[freq_t2_col]] - .data[[freq_t0_col]],
    hyper_T2 = .data[[freq_t2_col]] > 0.002
  )

## =========================
## 7) Select "important" clonotypes
##    A) Responders with T0: pre-existing (T0>0) and expanded (delta>0) Top N
##    B) Responders without T0: Top N largest clones at T2
## =========================
important_A <- cl_wide2 %>%
  filter(Response == "Responder", T0_available, T2_available) %>%
  filter(.data[[freq_t0_col]] > 0) %>%
  filter(delta > 0) %>%
  arrange(desc(delta), desc(.data[[freq_t2_col]])) %>%
  slice_head(n = N_TOP)

important_B <- cl_wide2 %>%
  filter(Response == "Responder", !T0_available, T2_available) %>%
  arrange(desc(.data[[freq_t2_col]])) %>%
  slice_head(n = N_TOP_MISSING_T0)

important <- bind_rows(important_A, important_B) %>%
  mutate(
    TRA_aa = if (has_CTaa) str_split_fixed(CTaa_mode, "_", 2)[,1] else NA_character_,
    TRB_aa = if (has_CTaa) str_split_fixed(CTaa_mode, "_", 2)[,2] else NA_character_,
    T0_available = ifelse(T0_available, "Yes", "No")
  )

write_tsv(important, file.path(OUTDIR, "important_clonotypes_Responder.tsv"))

## =========================
## 8) T_NK_cell composition at T2 for selected clonotypes
##    - Top 3 states + the rest as "Other"
## =========================
meta_sel <- meta_tcr %>%
  mutate(clone_id = .data[[clone_id_col]]) %>%
  semi_join(important %>% distinct(patients, clone_id),
            by = c("patients", "clone_id"))

tnk_prop_T2 <- meta_sel %>%
  filter(Timepoint == TIME_T2) %>%
  count(patients, clone_id, T_NK_cell, name = "n") %>%
  group_by(patients, clone_id) %>%
  mutate(prop = n / sum(n)) %>%
  arrange(patients, clone_id, desc(prop)) %>%
  mutate(rank = row_number()) %>%
  summarise(
    `T_NK_cell composition (T2)` = {
      top <- cur_data() %>%
        filter(rank <= 3) %>%
        mutate(txt = paste0(T_NK_cell, " ", sprintf("%.0f%%", 100 * prop)))
      
      other <- cur_data() %>%
        filter(rank > 3) %>%
        summarise(other_prop = sum(prop), .groups = "drop") %>%
        pull(other_prop)
      
      out <- top$txt
      if (length(other) == 1 && is.finite(other) && other > 0) {
        out <- c(out, paste0("Other ", sprintf("%.0f%%", 100 * other)))
      }
      paste(out, collapse = "; ")
    },
    .groups = "drop"
  )

## =========================
## 9) Figure table data (no log2FC/selection columns)
## =========================
table_df <- important %>%
  transmute(
    Patient = patients,
    CloneID = clone_id,
    
    `CDR3α (aa)` = ifelse(is.na(TRA_aa) | TRA_aa == "", "(NA)", TRA_aa),
    `CDR3β (aa)` = ifelse(is.na(TRB_aa) | TRB_aa == "", "(NA)", TRB_aa),
    
    `TRAV/J (from CTgene)` = if (has_CTgene) {
      str_extract(CTgene_mode, "^TRAV[^_]+") %>% replace_na("(NA)")
    } else "(NA)",
    
    `TRBV/J (from CTgene)` = if (has_CTgene) {
      str_extract(CTgene_mode, "TRBV.*$") %>% replace_na("(NA)")
    } else "(NA)",
    
    `T0 (%)` = 100 * .data[[freq_t0_col]],
    `T2 (%)` = 100 * .data[[freq_t2_col]],
    `ΔT2–T0 (%)` = 100 * delta,
    `Hyper @T2` = ifelse(hyper_T2, "Yes", "No"),
    `T0 sample available` = T0_available
  ) %>%
  left_join(tnk_prop_T2, by = c("Patient" = "patients", "CloneID" = "clone_id")) %>%
  mutate(
    `T_NK_cell composition (T2)` = replace_na(`T_NK_cell composition (T2)`, "(NA)"),
    Row = paste0(Patient, " | ", `CDR3β (aa)`, " | T0:", `T0 sample available`)
  ) %>%
  arrange(desc(`ΔT2–T0 (%)`), desc(`T2 (%)`))

## =========================
## 10) (Optional) Convert headers to ASCII to avoid Unicode warnings
## =========================
if (ASCII_HEADERS) {
  table_df <- table_df %>%
    rename(
      `CDR3a (aa)` = `CDR3α (aa)`,
      `CDR3b (aa)` = `CDR3β (aa)`,
      `dT2-T0 (%)` = `ΔT2–T0 (%)`
    )
  dcol <- "dT2-T0 (%)"
  col_levels <- c(
    "CDR3a (aa)",
    "CDR3b (aa)",
    "TRAV/J (from CTgene)",
    "TRBV/J (from CTgene)",
    "T0 (%)",
    "T2 (%)",
    dcol,
    "Hyper @T2",
    "T0 sample available",
    "T_NK_cell composition (T2)"
  )
} else {
  dcol <- "ΔT2–T0 (%)"
  col_levels <- c(
    "CDR3α (aa)",
    "CDR3β (aa)",
    "TRAV/J (from CTgene)",
    "TRBV/J (from CTgene)",
    "T0 (%)",
    "T2 (%)",
    dcol,
    "Hyper @T2",
    "T0 sample available",
    "T_NK_cell composition (T2)"
  )
}

## =========================
## 11) Table plot with per-column widths (geom_rect)
##     Numeric columns (T0/T2/Δ) are narrow; string columns are wider
## =========================
row_levels <- table_df$Row

long <- table_df %>%
  select(-Patient, -CloneID) %>%
  pivot_longer(
    cols = -Row,
    names_to = "col",
    values_to = "val",
    values_transform = list(val = as.character)
  ) %>%
  mutate(
    col = factor(col, levels = col_levels),
    val_num = suppressWarnings(as.numeric(val)),
    label = case_when(
      col %in% c("T0 (%)", "T2 (%)", dcol) ~ sprintf("%.2f", val_num),
      TRUE ~ val
    ),
    heat = ifelse(col == dcol, val_num, NA_real_)
  )

## Column width specification (adjust if needed)
col_widths <- c(
  "CDR3a (aa)" = 3.4,
  "CDR3b (aa)" = 3.4,
  "TRAV/J (from CTgene)" = 4,
  "TRBV/J (from CTgene)" = 3,
  "T_NK_cell composition (T2)" = 5,
  "T0 (%)" = 0.9,
  "T2 (%)" = 0.9,
  "Hyper @T2" = 0.9,
  "T0 sample available" = 0.9
)
col_widths[dcol] <- 1.0

## Assign a default width for any column missing from col_widths
missing_cols <- setdiff(col_levels, names(col_widths))
if (length(missing_cols) > 0) col_widths[missing_cols] <- 1.6

layout_df <- tibble(
  col = col_levels,
  w = as.numeric(col_widths[col_levels])
) %>%
  mutate(
    xmin = lag(cumsum(w), default = 0),
    xmax = cumsum(w),
    x = (xmin + xmax) / 2
  )

row_df <- tibble(Row = row_levels) %>%
  mutate(y = rev(seq_along(Row)))  # Place the first row at the top

plot_df <- long %>%
  left_join(layout_df, by = "col") %>%
  left_join(row_df, by = "Row") %>%
  mutate(
    ymin = y - 0.5,
    ymax = y + 0.5
  )

p_table <- ggplot(plot_df) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = heat),
    color = "grey85", linewidth = 0.35
  ) +
  geom_text(
    aes(x = x, y = y, label = label),
    size = 3, lineheight = 0.95
  ) +
  scale_fill_gradient2(
    low = "#2c7bb6", mid = "white", high = "#d7191c",
    midpoint = 0, na.value = "white", name = dcol
  ) +
  scale_x_continuous(
    breaks = layout_df$x,
    labels = layout_df$col,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = row_df$y,
    labels = row_df$Row,
    expand = c(0, 0)
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Responders: post-treatment expansion of tumor clonotypes",
    subtitle = paste0(
      "Clone frequencies computed within each sample; samples with < ",
      MIN_TCR_CELLS,
      " TCR+ cells removed. T_NK_cell composition shown for T2 (top3 + Other)."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9, angle = 35, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_table)

ggsave(file.path(OUTDIR, "important_clonotypes_table_with_TNKcomposition.pdf"),
       p_table, width = 20, height = 5.2)
ggsave(file.path(OUTDIR, "important_clonotypes_table_with_TNKcomposition.png"),
       p_table, width = 20, height = 5.2, dpi = 300)

## =========================
## 12) (Optional) Sanity check: selected clones T0 vs T2 (%)
## =========================
p_check <- table_df %>%
  select(Row, `T0 (%)`, `T2 (%)`) %>%
  pivot_longer(cols = c(`T0 (%)`, `T2 (%)`), names_to = "tp", values_to = "pct") %>%
  mutate(pct = as.numeric(pct)) %>%
  ggplot(aes(x = tp, y = pct, group = Row)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1.8) +
  coord_flip() +
  theme_minimal(base_size = 11) +
  labs(title = "Selected clonotypes: T0 vs T2 frequency (%)", x = NULL, y = "%")

print(p_check)
ggsave(file.path(OUTDIR, "important_clonotypes_T0_T2_check.pdf"), p_check, width = 10, height = 5)

table_df %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  print()

table_df %>%
  summarise(across(everything(), ~sum(as.character(.) %in% c("NA","<NA>","None")))) %>%
  print()

important %>%
  filter(!is.na(CTgene_mode), CTgene_mode != "") %>%
  mutate(trbv = str_extract(CTgene_mode, "TRBV.*$")) %>%
  filter(is.na(trbv)) %>%
  distinct(CTgene_mode) %>%
  print(n = Inf)

T_cell@meta.data[grepl("CASSTTDYNEQFF", T_cell@meta.data$CTaa), ]

filter(T_cell@meta.data, )
DefaultAssay(T_cell) <- "RNA"
FeaturePlot(T_cell, features = "TRAV13-2")
FeaturePlot(T_cell, features = "TRBV25-1")
Idents(T_cell) <- "T_NK_cell"
DotPlot(T_cell, c("TRAV13-2", "TRBV25-1"), group.by = "T_R")

library(scRepertoire)
scRep_example <- highlightClones(
  T_cell,
  cloneCall = "aa",
  sequence = c(T_cell@meta.data$CTaa)
)

Seurat::DimPlot(scRep_example, split.by = "highlight") +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  ggplot2::theme(
    plot.title = element_blank(),
    legend.position = "bottom"
  )

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

# =========================================================
# 1) Build the T0 clonotype set from tissue and label T2 cells as Maintained/New
#    -> Aggregate by Response to form a 2x2 contingency table (cell-weighted)
# =========================================================
make_fisher_table_T2 <- function(
    seu,
    clonotype_col = "CTaa",
    patient_col = "patients",
    response_col = "Response",
    time_col = "Timepoint",
    t0 = "T0",
    t2 = "T2"
) {
  md <- seu@meta.data %>%
    as.data.frame() %>%
    filter(
      .data[[time_col]] %in% c(t0, t2),
      !is.na(.data[[clonotype_col]]),
      .data[[clonotype_col]] != ""
    )

  # Patients with both T0 and T2 present + minimal clonotype-positive cell count criteria
  keep_pat <- md %>%
    count(.data[[patient_col]], .data[[time_col]], name = "n_cells") %>%
    group_by(.data[[patient_col]]) %>%
    summarise(
      has_both = all(c(t0, t2) %in% .data[[time_col]]),
      min_cells = min(n_cells),
      .groups = "drop"
    ) %>%
    pull(.data[[patient_col]])

  md <- md %>% filter(.data[[patient_col]] %in% keep_pat)

  # Per-patient T0 clonotype set
  t0_sets <- md %>%
    filter(.data[[time_col]] == t0) %>%
    group_by(.data[[patient_col]]) %>%
    summarise(t0_set = list(unique(.data[[clonotype_col]])), .groups = "drop")

  # Label each T2 cell as maintained/new
  t2_labeled <- md %>%
    filter(.data[[time_col]] == t2) %>%
    left_join(t0_sets, by = patient_col) %>%
    mutate(status = ifelse(.data[[clonotype_col]] %in% unlist(t0_set), "Maintained", "New"))

  # Per-patient summary (optional; can be used for dot plots etc.)
  per_patient <- t2_labeled %>%
    group_by(.data[[patient_col]], .data[[response_col]], status) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    pivot_wider(names_from = status, values_from = n_cells, values_fill = 0) %>%
    rename(patients = !!patient_col, Response = !!response_col)

  # Aggregate by Response -> 2x2 contingency table
  grp <- per_patient %>%
    group_by(Response) %>%
    summarise(
      Maintained = sum(Maintained),
      New = sum(New),
      .groups = "drop"
    )

  tab <- as.matrix(grp[, c("Maintained", "New")])
  rownames(tab) <- grp$Response

  list(per_patient = per_patient, group_table = tab)
}

# =========================================================
# 2) Fix row/column order (Non_Responder first)
# =========================================================
order_table <- function(tab,
                        row_order = c("Non_Responder", "Responder"),
                        col_order = c("Maintained", "New")) {
  # rows
  if (all(row_order %in% rownames(tab))) tab <- tab[row_order, , drop = FALSE]
  # columns
  if (all(col_order %in% colnames(tab))) tab <- tab[, col_order, drop = FALSE]
  tab
}

# =========================================================
# 3) Proportion stacked bar + annotate p/OR/CI
# =========================================================
plot_proportion_with_stats <- function(tab, title = "T2 clonotype status in tissue") {
  tab <- order_table(tab)

  ft <- fisher.test(tab, alternative = "two.sided")

  df <- as.data.frame(tab) %>%
    tibble::rownames_to_column("Response") %>%
    pivot_longer(cols = colnames(tab), names_to = "Status", values_to = "Count") %>%
    group_by(Response) %>%
    mutate(Total = sum(Count), Prop = Count / Total) %>%
    ungroup() %>%
    mutate(
      Response = factor(Response, levels = rownames(tab)),
      Status   = factor(Status, levels = colnames(tab)),
      label_in = paste0(Count, "\n(", percent(Prop, accuracy = 0.1), ")")
    )

  ann <- paste0(
    "Fisher's exact test (two-sided)\n",
    "p = ", format.pval(ft$p.value, digits = 2, eps = 1e-16), "\n",
    "OR = ", signif(unname(ft$estimate), 3),
    " (95% CI ", signif(ft$conf.int[1], 3), "–", signif(ft$conf.int[2], 3), ")\n",
    "Row totals: ", paste0(rownames(tab), "=", rowSums(tab), collapse = ", ")
  )

  p <- ggplot(df, aes(x = Response, y = Prop, fill = Status)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.25) +
    geom_text(
      aes(label = label_in),
      position = position_stack(vjust = 0.5),
      size = 3.2,
      lineheight = 0.95
    ) +
    scale_y_continuous(
      limits = c(0, 1.18),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = percent_format(accuracy = 1)
    ) +
    labs(
      title = title,
      subtitle = "Maintained vs New in T2 (cell-weighted)",
      x = "",
      y = "Proportion of T2 TCR+ cells"
    ) +
    annotate(
      "label",
      x = 1.5, y = 1.165,  # Place annotation so the p-value is clearly visible
      label = ann,
      size = 3.2,
      label.size = 0.25,
      label.r = unit(0.15, "lines"),
      hjust = 0.5, vjust = 1
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(margin = margin(b = 6)),
      axis.text.x = element_text(face = "bold"),
      plot.margin = margin(10, 16, 8, 8)
    )

  list(plot = p, fisher = ft, table = tab)
}

# =========================================================
# 5) RUN (execute only below)
# =========================================================
out <- make_fisher_table_T2(T_cell, clonotype_col = "CTaa")

tab <- order_table(out$group_table, row_order = c("Non_Responder", "Responder"))
res <- plot_proportion_with_stats(tab, title = "T2 clonotype status in tissue (CTaa)")

print(res$table)    # Check 2x2
print(res$fisher)   # Fisher results
print(res$plot)     # Proportion plot

# Save
ggsave("Figure_Fisher_T2_Maintained_vs_New_proportion.pdf", res$plot, width = 7.2, height = 4.9)
ggsave("Figure_Fisher_T2_Maintained_vs_New_proportion.png", res$plot, width = 7.2, height = 4.9, dpi = 600)

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
})

## =========================
## INPUT (must already exist)
##  - important: includes patients, clone_id (=CTnt), and CTaa_mode (or TRA_aa/TRB_aa)
##  - PBMC_T: Seurat object (PBMC)
## =========================

OUTDIR <- "plots/pbmc_clones"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## Label options
LABEL_MIN_PROP <- 0.05   # Label only >=5%; set to 0 to label all
LABEL_DIGITS   <- 1
TIMEPOINT_LEVELS <- c("T0", "T1", "T2", "T3")

## ---------------------------------------------------------
## 0) Validate keys: clone_id should be CTnt
## ---------------------------------------------------------
stopifnot(all(c("patients","clone_id") %in% colnames(important)))

## If patients have a trailing underscore (e.g., RCC5_), IDs may not match PBMC
## Set TRUE to normalize patient IDs on both sides
NORMALIZE_PATIENT_ID <- TRUE
norm_patient <- function(x) {
  x <- as.character(x)
  if (NORMALIZE_PATIENT_ID) x <- str_replace(x, "_$", "")
  x
}

important_key <- important %>%
  transmute(
    Patient = norm_patient(patients),
    CloneID = as.character(clone_id),        # <-- CTnt
    CTaa_label = if ("CTaa_mode" %in% colnames(important)) as.character(CTaa_mode) else NA_character_,
    TRA_aa = if ("TRA_aa" %in% colnames(important)) as.character(TRA_aa) else NA_character_,
    TRB_aa = if ("TRB_aa" %in% colnames(important)) as.character(TRB_aa) else NA_character_
  ) %>%
  distinct()

## ---------------------------------------------------------
## 1) Prepare PBMC meta
## ---------------------------------------------------------
pbmc_meta <- PBMC_T@meta.data %>% rownames_to_column("barcode")

need_cols <- c("patients","Timepoint","Response","Merged","CTnt","CTaa")
miss <- setdiff(need_cols, colnames(pbmc_meta))
if (length(miss) > 0) stop("PBMC_T@meta.data is missing required columns: ", paste(miss, collapse = ", "))

pbmc_tcr <- pbmc_meta %>%
  mutate(
    Patient = norm_patient(patients),
    Timepoint = as.character(Timepoint),
    CellType  = as.character(Merged),
    CloneID_CTnt = as.character(CTnt),
    CTaa_pbmc = as.character(CTaa)
  ) %>%
  filter(!is.na(CloneID_CTnt), CloneID_CTnt != "") %>%
  filter(!is.na(CellType), CellType != "") %>%
  filter(!is.na(Timepoint), Timepoint != "")

## ---------------------------------------------------------
## 2) Identify tumor-selected clones in PBMC (Patient + CTnt)
## ---------------------------------------------------------
pbmc_sel <- pbmc_tcr %>%
  transmute(
    barcode,
    Patient,
    CloneID = CloneID_CTnt,     # Match by CTnt
    Timepoint,
    Response = as.character(Response),
    CellType,
    CTaa_pbmc
  ) %>%
  semi_join(important_key %>% select(Patient, CloneID), by = c("Patient","CloneID"))

found_key <- pbmc_sel %>% distinct(Patient, CloneID)
missing_key <- anti_join(
  important_key %>% distinct(Patient, CloneID),
  found_key,
  by = c("Patient","CloneID")
)

message("Number of clones found in PBMC: ", nrow(found_key))
message("Number of clones NOT found in PBMC: ", nrow(missing_key))

write_tsv(found_key,   file.path(OUTDIR, "found_clones_in_PBMC.tsv"))
write_tsv(missing_key, file.path(OUTDIR, "missing_clones_in_PBMC.tsv"))

if (nrow(pbmc_sel) == 0) stop("No clones matched in PBMC. Re-check patient IDs / CTnt column / duplicates.")

tp_use <- intersect(TIMEPOINT_LEVELS, sort(unique(pbmc_sel$Timepoint)))
pbmc_sel <- pbmc_sel %>% mutate(Timepoint = factor(Timepoint, levels = tp_use))

## ---------------------------------------------------------
## 3) Build clone labels for figures (display CTaa)
##    - Prefer CTaa_mode from important if available
## ---------------------------------------------------------
clone_labels <- important_key %>%
  mutate(
    CloneLabel = case_when(
      !is.na(CTaa_label) & CTaa_label != "" ~ paste0("CTaa: ", CTaa_label),
      !is.na(TRA_aa) & !is.na(TRB_aa) & TRA_aa != "" & TRB_aa != "" ~ paste0("CTaa: ", TRA_aa, "_", TRB_aa),
      TRUE ~ paste0("CTnt: ", str_sub(CloneID, 1, 18), "…")
    )
  ) %>%
  select(Patient, CloneID, CloneLabel) %>%
  distinct()

pbmc_sel <- pbmc_sel %>%
  left_join(clone_labels, by = c("Patient","CloneID"))

## ---------------------------------------------------------
## 4) Figure A: clone-level composition in PBMC (which subtypes host each clone?)
##    - x: Timepoint, y: within-clone proportion
##    - fill: Merged (CellType)
##    - label: n (xx%)
## ---------------------------------------------------------
pbmc_clone_comp <- pbmc_sel %>%
  count(Patient, CloneID, CloneLabel, Timepoint, CellType, name = "n") %>%
  group_by(Patient, CloneID, CloneLabel, Timepoint) %>%
  mutate(
    total = sum(n),
    prop = n / total,
    label = ifelse(
      prop >= LABEL_MIN_PROP,
      paste0(n, " (", sprintf(paste0("%.", LABEL_DIGITS, "f"), 100 * prop), "%)"),
      ""
    )
  ) %>%
  ungroup()

p_clone <- ggplot(pbmc_clone_comp, aes(x = Timepoint, y = prop, fill = CellType)) +
  geom_col(width = 0.85) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 2.8, color = "black", lineheight = 0.9) +
  facet_grid(CloneLabel ~ Patient, scales = "free_y", space = "free_y") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(size = 7),
    legend.position = "bottom"
  ) +
  labs(
    title = "PBMC phenotypes of tumor-selected clonotypes (matched by CTnt, labeled by CTaa)",
    subtitle = paste0("Labels: n (%, within clonotype/timepoint). Label threshold: ", 100 * LABEL_MIN_PROP, "%"),
    x = "Timepoint",
    y = "Within-clonotype proportion",
    fill = "Merged"
  )

ggsave(file.path(OUTDIR, "PBMC_clone_level_stacked_by_Merged_labeled.pdf"),
       p_clone, width = 12, height = 9)
ggsave(file.path(OUTDIR, "PBMC_clone_level_stacked_by_Merged_labeled.png"),
       p_clone, width = 12, height = 9, dpi = 300)

print(p_clone)

## ---------------------------------------------------------
## 5) Figure B: pooled composition (aggregate selected clones per patient/timepoint)
## ---------------------------------------------------------
pbmc_patient_comp <- pbmc_sel %>%
  count(Patient, Timepoint, CellType, name = "n") %>%
  group_by(Patient, Timepoint) %>%
  mutate(
    total = sum(n),
    prop = n / total,
    label = ifelse(
      prop >= LABEL_MIN_PROP,
      paste0(n, " (", sprintf(paste0("%.", LABEL_DIGITS, "f"), 100 * prop), "%)"),
      ""
    )
  ) %>%
  ungroup()

p_pooled <- ggplot(pbmc_patient_comp, aes(x = Timepoint, y = prop, fill = CellType)) +
  geom_col(width = 0.85) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 3.0, color = "black", lineheight = 0.9) +
  facet_wrap(~ Patient, nrow = 1) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    title = "PBMC phenotypes of tumor-selected clonotypes (pooled across clonotypes)",
    subtitle = paste0("Labels: n (%, within patient/timepoint). Label threshold: ", 100 * LABEL_MIN_PROP, "%"),
    x = "Timepoint",
    y = "Proportion",
    fill = "Merged"
  )

ggsave(file.path(OUTDIR, "PBMC_patient_timepoint_pooled_stacked_by_Merged_labeled.pdf"),
       p_pooled, width = 12, height = 4.2)
ggsave(file.path(OUTDIR, "PBMC_patient_timepoint_pooled_stacked_by_Merged_labeled.png"),
       p_pooled, width = 12, height = 4.2, dpi = 300)

print(p_pooled)

## ---------------------------------------------------------
## 6) Presence table by timepoint
## ---------------------------------------------------------
presence_tbl <- pbmc_sel %>%
  distinct(Patient, CloneID, Timepoint) %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = Timepoint, values_from = present, values_fill = 0L)

write_tsv(presence_tbl, file.path(OUTDIR, "PBMC_clone_presence_by_timepoint.tsv"))
presence_tbl

pbmc_sel

# =========================
# Settings
# =========================
outdir <- "plots/pbmc_origin_in_tumor"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tcr_col <- "CTstrict"
tumor_state_col <- "T_NK_cell"
pbmc_tp_col <- "Timepoint"
tumor_tp_col <- "Timepoint"

use_tps <- c("T0", "T2")

pbmc_labels <- list(pEx = "CD8_pEx", Temra = "CD8_Temra")
pbmc_t1_prolif_label <- "Proliferative T cell"

stopifnot(all(c(tcr_col, "Merged", pbmc_tp_col) %in% colnames(PBMC_T@meta.data)))
stopifnot(all(c(tcr_col, tumor_state_col, "Response", tumor_tp_col) %in% colnames(T_cell@meta.data)))

tnk_palette <- c(
  "CD56-high_NK" = "#9BC8BE",
  "CD56-low_NK"  = "#F5F0B9",
  "CD8_STR"      = "#B9B4CD",
  "tEx"          = "#D27D73",
  "ILC"          = "#8CA5C3",
  "iEx"          = "#E1AF6E",
  "Naive T cell" = "#B4CD78",
  "pEx"          = "#EBC8DC",
  "Proliferative T cell" = "#D7D7D7",
  "Temra"        = "#A57DAF",
  "Treg"         = "#CDE1C3"
)

normalize_state <- function(x) {
  dplyr::case_when(
    x == "Exhausted T cell(tEx)" ~ "tEx",
    x == "Intermediate exhausted T cell(iEx)" ~ "iEx",
    x == "Progenitor exhausted T cell(pEx)" ~ "pEx",
    TRUE ~ as.character(x)
  )
}

chisq_p <- function(tab) {
  tab <- as.matrix(tab)
  if (nrow(tab) < 2 || ncol(tab) < 2) return(list(p = NA_real_, warn = NA_character_))
  w <- NULL
  p <- withCallingHandlers(
    chisq.test(tab)$p.value,
    warning = function(e) { w <<- conditionMessage(e); invokeRestart("muffleWarning") }
  )
  list(p = unname(p), warn = w)
}

# =========================================================
# (1) Same-timepoint mapping: PBMC(T0/T2) -> Tumor(T0/T2)
# =========================================================
pbmc_md <- PBMC_T@meta.data %>%
  mutate(
    pbmc_tp = as.character(.data[[pbmc_tp_col]]),
    clone = .data[[tcr_col]]
  ) %>%
  filter(pbmc_tp %in% use_tps, !is.na(clone), clone != "")

tumor_md <- T_cell@meta.data %>%
  mutate(
    tumor_tp = as.character(.data[[tumor_tp_col]]),
    clone = .data[[tcr_col]],
    tumor_state = normalize_state(.data[[tumor_state_col]]),
    Response = as.character(Response)
  ) %>%
  filter(tumor_tp %in% use_tps, !is.na(clone), clone != "") %>%
  filter(!is.na(tumor_state), tumor_state != "", !is.na(Response), Response != "")

make_origin_map <- function(tp) {
  pex <- pbmc_md %>%
    filter(pbmc_tp == tp, Merged == pbmc_labels$pEx) %>%
    distinct(clone) %>% pull(clone)

  temra <- pbmc_md %>%
    filter(pbmc_tp == tp, Merged == pbmc_labels$Temra) %>%
    distinct(clone) %>% pull(clone)

  both <- intersect(pex, temra)

  tibble(
    tumor_tp = tp,
    clone = unique(c(pex, temra)),
    origin = case_when(
      clone %in% both ~ "Both",
      clone %in% pex ~ "pEx",
      clone %in% temra ~ "Temra",
      TRUE ~ NA_character_
    )
  ) %>%
    filter(origin %in% c("pEx", "Temra"))
}

origin_map <- bind_rows(lapply(use_tps, make_origin_map))

tumor_matched <- tumor_md %>%
  left_join(origin_map, by = c("tumor_tp", "clone")) %>%
  filter(origin %in% c("pEx", "Temra"))

if (nrow(tumor_matched) == 0) stop("No pEx/Temra tumor cells after same-timepoint mapping.")

write.table(
  tumor_matched %>% select(clone, origin, tumor_tp, Response, tumor_state),
  file.path(outdir, "tumor_cells_sameTP_matched.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -------------------------
# TEST #1: Per-timepoint origin × state (chi-square) + BH adjust across timepoints
# -------------------------
chisq_overall_tp <- tumor_matched %>%
  count(tumor_tp, origin, tumor_state, name = "n") %>%
  group_by(tumor_tp) %>%
  group_modify(~{
    tab <- xtabs(n ~ origin + tumor_state, data = .x)
    res <- chisq_p(tab)
    tibble(p_value = res$p, warn = res$warn)
  }) %>%
  ungroup() %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    tp_label = paste0(
      tumor_tp, " (origin × state) BH p=",
      format.pval(p_adj_BH, digits = 4, eps = 1e-300)
    )
  )

write.table(
  chisq_overall_tp,
  file.path(outdir, "chisq_by_timepoint_origin_vs_state_BH.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

plot_overall_tp <- tumor_matched %>%
  count(tumor_tp, origin, tumor_state, name = "n") %>%
  group_by(tumor_tp, origin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

strip_lab <- setNames(chisq_overall_tp$tp_label, chisq_overall_tp$tumor_tp)

pA <- ggplot(plot_overall_tp, aes(origin, prop, fill = tumor_state)) +
  geom_col(width = 0.8) +
  facet_wrap(~ tumor_tp, nrow = 1, labeller = labeller(tumor_tp = strip_lab)) +
  scale_fill_manual(values = tnk_palette, drop = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    title = "Same-timepoint mapping (T0/T2): tumor-state by PBMC origin, split by tumor timepoint",
    x = "PBMC origin (within same TP)",
    y = "Proportion (within origin)",
    fill = "Tumor state"
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "PLOT_A_byTP_origin_vs_state_chisq_BH.pdf"),
       pA, width = 10, height = 4.5)
ggsave(file.path(outdir, "PLOT_A_byTP_origin_vs_state_chisq_BH.png"),
       pA, width = 10, height = 4.5, dpi = 300)
pA

# -------------------------
# TEST #2: Per (origin, TP) Response × state (chi-square) + BH within origin
# -------------------------
chisq_by_origin_tp <- tumor_matched %>%
  count(origin, tumor_tp, Response, tumor_state, name = "n") %>%
  group_by(origin, tumor_tp) %>%
  group_modify(~{
    tab <- xtabs(n ~ Response + tumor_state, data = .x)
    res <- chisq_p(tab)
    tibble(p_value = res$p, warn = res$warn)
  }) %>%
  ungroup() %>%
  group_by(origin) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(label = paste0(
    "BH p=",
    format.pval(p_adj_BH, digits = 4, eps = 1e-300)
  ))

write.table(
  chisq_by_origin_tp,
  file.path(outdir, "chisq_by_origin_tp_Response_vs_state.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

plot_B <- tumor_matched %>%
  count(origin, tumor_tp, Response, tumor_state, name = "n") %>%
  group_by(origin, tumor_tp, Response) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    Response = factor(Response, levels = c("Non_Responder", "Responder")),
    origin = factor(origin, levels = c("pEx", "Temra"))
  )

pB <- ggplot(plot_B, aes(Response, prop, fill = tumor_state)) +
  geom_col(width = 0.8) +
  facet_grid(origin ~ tumor_tp) +
  geom_text(
    data = chisq_by_origin_tp,
    aes(x = 1.5, y = 1.06, label = label),
    inherit.aes = FALSE,
    size = 3.1
  ) +
  scale_fill_manual(values = tnk_palette, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  theme_bw(base_size = 12) +
  labs(
    title = "Same-timepoint mapping (T0/T2): response shifts within each origin & tumor timepoint",
    subtitle = "Per (origin, tumor TP): Chi-square on Response×tumor_state; BH adjusted within origin.",
    x = "Response",
    y = "Proportion within (origin, TP, response)",
    fill = "Tumor state"
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "PLOT_B_origin_tp_xResponse_chisq.pdf"), pB, width = 11, height = 6.5)
ggsave(file.path(outdir, "PLOT_B_origin_tp_xResponse_chisq.png"), pB, width = 11, height = 6.5, dpi = 300)

pA
pB

## -----------------------------
## 1) meta
## -----------------------------
T_cell@meta.data$T_NK_cell %>% table()
T_cell <- T_cell_TCR

## Normalize tumor-state labels to short keys
T_cell@meta.data$T_NK_cell <- ifelse(T_cell@meta.data$T_NK_cell == "Exhausted T cell(tEx)", "tEx", T_cell@meta.data$T_NK_cell)
T_cell@meta.data$T_NK_cell <- ifelse(T_cell@meta.data$T_NK_cell == "Intermediate exhausted T cell(iEx)", "iEx", T_cell@meta.data$T_NK_cell)
T_cell@meta.data$T_NK_cell <- ifelse(T_cell@meta.data$T_NK_cell == "Progenitor exhausted T cell(pEx)", "pEx", T_cell@meta.data$T_NK_cell)

meta <- T_cell@meta.data %>% filter(Timepoint == "T2") %>% rownames_to_column("barcode")

clone_id_col <- if ("CTnt" %in% colnames(meta)) "CTnt" else "CTstrict"
stopifnot(clone_id_col %in% colnames(meta))
stopifnot(all(c("T_NK_cell","sample_name") %in% colnames(meta)))

meta <- meta %>% mutate(sample_id = sample_name)

## Keep only TCR+ cells
meta_tcr <- meta %>%
  filter(!is.na(.data[[clone_id_col]]), .data[[clone_id_col]] != "") %>%
  filter(!is.na(sample_id), sample_id != "") %>%
  filter(!is.na(T_NK_cell), T_NK_cell != "")

## -----------------------------
## 2) Sample depth cutoff (>= 500 TCR+ cells)
## -----------------------------
good_samples <- meta_tcr %>%
  dplyr::count(sample_id, name = "TCR_cells") %>%
  filter(TCR_cells >= min_sample_tcr_cells) %>%
  pull(sample_id)

meta_tcr <- meta_tcr %>% filter(sample_id %in% good_samples)

## -----------------------------
## 3) Paper-style: within-sample clonotype frequency
## -----------------------------
clone_sample <- meta_tcr %>%
  group_by(sample_id, clone_id = .data[[clone_id_col]]) %>%
  summarise(
    n_cells = n(),
    T_NK_cell_mode = names(sort(table(T_NK_cell), decreasing = TRUE))[1],
    .groups = "drop"
  ) %>%
  group_by(sample_id) %>%
  mutate(freq_in_sample = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  mutate(
    expanded = ifelse(
      freq_in_sample > hyper_freq_cut & n_cells >= min_clone_cells_for_hyper,
      "Hyper-Expanded", "Non-Expanded"
    )
  )

## -----------------------------
## 4) Re-normalize for visualization (hyper vs non)
##    - "Plot-only" step; matches the paper bubble logic
## -----------------------------
clone_group <- clone_sample %>%
  group_by(expanded, clone_id) %>%
  summarise(
    n_cells_group = sum(n_cells),
    T_NK_cell_mode = names(sort(table(T_NK_cell_mode), decreasing = TRUE))[1],
    .groups = "drop"
  ) %>%
  group_by(expanded) %>%
  mutate(Proportion = n_cells_group / sum(n_cells_group)) %>%
  ungroup() %>%
  mutate(group_key = expanded)

## -----------------------------
## 5) Helper to generate bubble polygons
## -----------------------------
make_bubble_plots <- function(df, prop_cut = 0.90, topn = 1e6,
                              area_scale = 1e9, npoints = 20) {
  function(idx) {

    tmp <- df %>%
      filter(group_key == idx) %>%
      arrange(desc(Proportion)) %>%
      mutate(cs = cumsum(Proportion) - Proportion) %>%
      filter(cs < prop_cut) %>%
      slice_head(n = topn) %>%
      select(group_key, clone_id, Proportion)

    if (nrow(tmp) == 0) return(tibble())

    pack <- circleProgressiveLayout(tmp$Proportion * area_scale, sizetype = "area")

    circleLayoutVertices(pack, npoints = npoints) %>%
      mutate(
        clone_id = rep(tmp$clone_id, each = npoints + 1),
        Proportion = rep(tmp$Proportion, each = npoints + 1)
      )
  }
}

make_bubble <- make_bubble_plots(
  clone_group %>% select(group_key, clone_id, Proportion),
  prop_cut = prop_cut_for_bubbles,
  area_scale = area_scale,
  npoints = npoints
)

bubble_plot_df <- clone_group %>%
  distinct(group_key) %>%
  mutate(plot_df = map(group_key, make_bubble)) %>%
  unnest(plot_df) %>%
  left_join(
    clone_group %>% select(group_key, clone_id, T_NK_cell_mode),
    by = c("group_key", "clone_id")
  ) %>%
  mutate(
    expanded = factor(group_key, levels = c("Hyper-Expanded", "Non-Expanded"))
  )

## -----------------------------
## 6) Final plot
## -----------------------------
p_bubble <- ggplot(bubble_plot_df) +
  geom_polygon(aes(x, y, group = id, fill = T_NK_cell_mode),
               color = "white", alpha = 1, linewidth = 0.01) +
  scale_fill_manual(values = tnk_palette, drop = FALSE) +
  facet_wrap(~expanded, ncol = 2) +
  coord_fixed(expand = FALSE) +
  theme_void() +
  labs(fill = NULL) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.8, "lines"),
    legend.text = element_text(size = 7),
    strip.text = element_text(size = 9)
  )

print(p_bubble)

ggsave(file.path(OUTDIR, "hyper_vs_non_TNK_samplefreq_cut0.002_depth500.pdf"),
       p_bubble, width = 120, height = 70, units = "mm")
ggsave(file.path(OUTDIR, "hyper_vs_non_TNK_samplefreq_cut0.002_depth500.png"),
       p_bubble, width = 120, height = 70, units = "mm", dpi = PNG_dpi)

## -----------------------------
## 7) QC output (verify sample-level hyper fraction)
## -----------------------------
qc_sample <- clone_sample %>%
  group_by(sample_id) %>%
  summarise(
    TCR_cells = sum(n_cells),
    n_clones = n(),
    frac_hyper = mean(freq_in_sample > hyper_freq_cut & n_cells >= min_clone_cells_for_hyper),
    .groups = "drop"
  ) %>%
  arrange(TCR_cells)

print(qc_sample, n = 50)







suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(packcircles)
})

PNG_dpi <- 300
OUTDIR <- "plots/bubble"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## Fixed colors
tnk_palette <- c(
  "CD56-high_NK" = "#9BC8BE",
  "CD56-low_NK"  = "#F5F0B9",
  "CD8_STR"      = "#B9B4CD",
  "Exhausted T cell(tEx)"   = "#D27D73",
  "ILC"          = "#8CA5C3",
  "Intermediate exhausted T cell(iEx)" = "#E1AF6E",
  "Naive T cell" = "#B4CD78",
  "Progenitor exhausted T cell(pEx)"   = "#EBC8DC",
  "Proliferative T cell" = "#D7D7D7",
  "Temra"        = "#A57DAF",
  "Treg"         = "#CDE1C3"
)

## -----------------------------
## Thresholds
## -----------------------------
min_sample_tcr_cells <- 500      # Per-sample TCR+ cell cutoff
hyper_freq_cut <- 0.002          # 0.2% (strictly >)
min_clone_cells_for_hyper <- 2   # Set to 1 if desired

prop_cut_for_bubbles <- 0.90     # Draw up to 90% cumulative proportion per panel
area_scale <- 1e9
npoints <- 20

## -----------------------------
## 1) meta
## -----------------------------
meta <- T_cell@meta.data %>% rownames_to_column("barcode")

clone_id_col <- if ("CTnt" %in% colnames(meta)) "CTnt" else "CTstrict"
stopifnot(clone_id_col %in% colnames(meta))
stopifnot(all(c("T_NK_cell","T_R","sample_name") %in% colnames(meta)))

meta <- meta %>% mutate(sample_id = sample_name)

meta_tcr <- meta %>%
  filter(!is.na(.data[[clone_id_col]]), .data[[clone_id_col]] != "") %>%
  filter(!is.na(sample_id), sample_id != "") %>%
  filter(!is.na(T_NK_cell), T_NK_cell != "") %>%
  filter(!is.na(T_R), T_R != "")

## -----------------------------
## 2) Sample depth cutoff (>= 500 TCR+ cells)
## -----------------------------
good_samples <- meta_tcr %>%
  count(sample_id, name = "TCR_cells") %>%
  filter(TCR_cells >= min_sample_tcr_cells) %>%
  pull(sample_id)

meta_tcr <- meta_tcr %>% filter(sample_id %in% good_samples)

## -----------------------------
## 3) Paper-style: within-sample clonotype frequency + hyper call
## -----------------------------
clone_sample <- meta_tcr %>%
  group_by(sample_id, clone_id = .data[[clone_id_col]]) %>%
  summarise(
    n_cells = n(),
    T_NK_cell_mode = names(sort(table(T_NK_cell), decreasing = TRUE))[1],
    T_R_mode       = names(sort(table(T_R),      decreasing = TRUE))[1],
    .groups = "drop"
  ) %>%
  group_by(sample_id) %>%
  mutate(freq_in_sample = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  mutate(
    expanded = ifelse(
      freq_in_sample > hyper_freq_cut & n_cells >= min_clone_cells_for_hyper,
      "Hyper-Expanded", "Non-Expanded"
    )
  )

## -----------------------------
## 4) Aggregate clones within each panel (T_R_mode × expanded)
##    and re-normalize within panel (plot-only)
## -----------------------------
clone_panel <- clone_sample %>%
  group_by(T_R_mode, expanded, clone_id) %>%
  summarise(
    n_cells_panel = sum(n_cells),
    T_NK_cell_mode = names(sort(table(T_NK_cell_mode), decreasing = TRUE))[1],
    .groups = "drop"
  ) %>%
  group_by(T_R_mode, expanded) %>%
  mutate(
    Proportion = n_cells_panel / sum(n_cells_panel),
    group_key = paste(T_R_mode, expanded, sep = "__")
  ) %>%
  ungroup()

## Fix factor levels (force all panels to show)
tr_levels <- meta_tcr %>%
  distinct(T_R) %>% pull(T_R) %>% unique()

clone_panel <- clone_panel %>%
  mutate(
    T_R_mode = factor(T_R_mode, levels = tr_levels),
    expanded = factor(expanded, levels = c("Hyper-Expanded", "Non-Expanded"))
  )

## -----------------------------
## 5) Helper to generate bubble polygons
## -----------------------------
make_bubble_plots <- function(df, prop_cut = 0.90, topn = 1e6,
                              area_scale = 1e9, npoints = 20) {
  function(idx) {

    tmp <- df %>%
      filter(group_key == idx) %>%
      arrange(desc(Proportion)) %>%
      mutate(cs = cumsum(Proportion) - Proportion) %>%
      filter(cs < prop_cut) %>%
      slice_head(n = topn) %>%
      select(group_key, clone_id, Proportion)

    if (nrow(tmp) == 0) return(tibble())

    pack <- circleProgressiveLayout(tmp$Proportion * area_scale, sizetype = "area")

    circleLayoutVertices(pack, npoints = npoints) %>%
      mutate(
        clone_id = rep(tmp$clone_id, each = npoints + 1),
        Proportion = rep(tmp$Proportion, each = npoints + 1)
      )
  }
}

make_bubble <- make_bubble_plots(
  clone_panel %>% select(group_key, clone_id, Proportion),
  prop_cut = prop_cut_for_bubbles,
  area_scale = area_scale,
  npoints = npoints
)

bubble_plot_df <- clone_panel %>%
  distinct(T_R_mode, expanded, group_key) %>%
  mutate(plot_df = map(group_key, make_bubble)) %>%
  unnest(plot_df) %>%
  left_join(
    clone_panel %>% select(group_key, clone_id, T_NK_cell_mode, T_R_mode, expanded),
    by = c("group_key", "clone_id", "T_R_mode", "expanded")
  )

## -----------------------------
## 6) Final plot: rows = T_R_mode, cols = expanded
## -----------------------------
p_bubble_TR <- ggplot(bubble_plot_df) +
  geom_polygon(aes(x, y, group = id, fill = T_NK_cell_mode),
               color = "white", alpha = 1, linewidth = 0.01) +
  scale_fill_manual(values = tnk_palette, drop = FALSE) +
  facet_grid(rows = vars(T_R_mode), cols = vars(expanded), drop = FALSE) +
  coord_fixed(expand = FALSE) +
  theme_void() +
  labs(fill = NULL) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.8, "lines"),
    legend.text = element_text(size = 7),
    strip.text = element_text(size = 8)
  )

print(p_bubble_TR)

ggsave(file.path(OUTDIR, "hyper_vs_non_TNK_by_TR_paper_method_depth500.pdf"),
       p_bubble_TR, width = 160, height = 180, units = "mm")
ggsave(file.path(OUTDIR, "hyper_vs_non_TNK_by_TR_paper_method_depth500.png"),
       p_bubble_TR, width = 160, height = 180, units = "mm", dpi = PNG_dpi)

## -----------------------------
## 7) QC (optional)
## -----------------------------
qc_sample <- clone_sample %>%
  group_by(sample_id) %>%
  summarise(
    TCR_cells = sum(n_cells),
    n_clones  = n(),
    frac_hyper = mean(freq_in_sample > hyper_freq_cut & n_cells >= min_clone_cells_for_hyper),
    .groups = "drop"
  ) %>%
  arrange(TCR_cells)

print(qc_sample, n = 50)








suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(dittoSeq)
})

## -----------------------
## Thresholds
## -----------------------
min_sample_tcr_cells <- 500
hyper_freq_cut <- 0.002          # 0.2% (strictly >)
min_clone_cells_for_hyper <- 2   # Set to 1 if desired

## -----------------------
## Extract meta
## -----------------------
meta <- T_cell@meta.data %>% rownames_to_column("barcode")

clone_id_col <- if ("CTnt" %in% colnames(meta)) "CTnt" else "CTstrict"
stopifnot(all(c("sample_name","T_NK_cell") %in% colnames(meta)))
stopifnot(clone_id_col %in% colnames(meta))

meta <- meta %>% mutate(sample_id = sample_name)

## Keep only TCR+ cells (for labeling)
meta_tcr <- meta %>%
  filter(!is.na(.data[[clone_id_col]]), .data[[clone_id_col]] != "") %>%
  filter(!is.na(sample_id), sample_id != "")

## Sample depth filter (>= 500)
good_samples <- meta_tcr %>%
  count(sample_id, name = "TCR_cells") %>%
  filter(TCR_cells >= min_sample_tcr_cells) %>%
  pull(sample_id)

meta_tcr <- meta_tcr %>% filter(sample_id %in% good_samples)

## -----------------------
## Paper-style hyper call (within sample)
## -----------------------
clone_sample <- meta_tcr %>%
  group_by(sample_id, clone_id = .data[[clone_id_col]]) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(freq_in_sample = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  mutate(
    expanded = ifelse(
      freq_in_sample > hyper_freq_cut & n_cells >= min_clone_cells_for_hyper,
      "Hyper-Expanded", "Non-Expanded"
    )
  )

## Build cell-level expanded labels
cell_expanded <- meta_tcr %>%
  transmute(
    barcode,
    sample_id,
    clone_id = .data[[clone_id_col]]
  ) %>%
  left_join(clone_sample %>% select(sample_id, clone_id, expanded),
            by = c("sample_id","clone_id")) %>%
  select(barcode, expanded)

## -----------------------
## Attach to Seurat meta.data
## -----------------------
T_cell$expanded <- NA_character_
T_cell$expanded[cell_expanded$barcode] <- cell_expanded$expanded
T_cell$expanded <- factor(T_cell$expanded, levels = c("Non-Expanded","Hyper-Expanded"))

## -----------------------
## Downstream plotting
## -----------------------
NA_T <- subset(T_cell, expanded != "NA")
NA_T_2 <- subset(T_cell, expanded == "Hyper-Expanded")

dittoSeq::dittoBarPlot(
  NA_T_2,
  var = "T_NK_cell",
  group.by = "T_R",
  color.panel = tnk_palette
)





library(dplyr)
library(ggplot2)
library(scales)   

overlap_obj <- clonalOverlap(
  T_PBMC_T,
  cloneCall = "strict",
  method = "jaccard",
  group.by = "T_NK_cell"
)

df_overlap <- overlap_obj@data

# Tissue(Var1 = T_) × PBMC(Var2 = B_) 
df_tb <- df_overlap %>%
  filter(
    grepl("^Tissue_", Var2),
    grepl("^PBMC_", Var1)
  )

p_overlap_pbmc_tissue <- ggplot(df_tb,
                                aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "grey70") +
  scale_fill_gradient(
    low  = "white",
    high = "red",
    na.value = "white",
    limits = c(min(df_tb$value, na.rm = TRUE),
               max(df_tb$value, na.rm = TRUE))
  ) +
  theme_bw() +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "Clonal overlap (Jaccard): Tissue (Y) vs PBMC (X)",
    x     = "PBMC T/NK subsets",
    y     = "Tissue T/NK subsets",
    fill  = "Overlap jaccard Score"
  )

p_overlap_pbmc_tissue
