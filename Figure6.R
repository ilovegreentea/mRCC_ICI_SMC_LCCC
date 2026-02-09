# ---- Library paths and package loading ----
.libPaths("~/R/rhel9/4.4.0")

library(Seurat)
library(SingleR)
library(dplyr)
library(SingleCellExperiment)

print("PACKAGE LOADED")

# ---- Load Seurat objects ----
# Combined Xenium query object
combined <- readRDS("/work/users/k/y/kyutae/re_IO/script/combined.rds")

# Reference scRNA-seq object with known subtypes
RCC_IO <- readRDS("/work/users/k/y/kyutae/Shared/RCC_IO_tissue_1030.rds")

# ---- Metadata preprocessing ----
# Extract batch information from cell barcodes (first 4 characters)
combined@meta.data$batch <- substr(rownames(combined@meta.data), 1, 4)

# Set Xenium assay as default for downstream processing
DefaultAssay(combined) <- "Xenium"

# ---- Split object by sample ----
# Each sample will be processed independently before integration
combined.list <- SplitObject(combined, split.by = "sample_name")

# Increase allowed memory for future-based parallelization
options(future.globals.maxSize = 80000 * 1024^2)

# ---- Per-sample normalization and feature selection ----
combined.list <- lapply(combined.list, NormalizeData)

combined.list <- lapply(
  combined.list,
  FindVariableFeatures,
  selection.method = "vst",
  nfeatures = 500
)

# ---- Sketching (cell subsampling) ----
# Uniformly sample up to 2,000 cells per sample
# This reduces computation while preserving structure
combined.list <- lapply(
  combined.list,
  SketchData,
  ncells = 2000,
  method = "Uniform",
  sketched.assay = "sketch"
)

# ---- Re-normalize and find variable features on sketched data ----
combined.list <- lapply(combined.list, function(obj) {
  DefaultAssay(obj) <- "sketch"
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 500
  )
  obj
})

# ---- Merge all samples back together ----
combined <- merge(
  x = combined.list[[1]],
  y = combined.list[-1]
)

# ---- Scaling and dimensional reduction ----
combined <- ScaleData(combined)

# PCA on sketched assay
combined <- RunPCA(combined, npcs = 50)

# Batch correction using Harmony (by sample)
combined <- RunHarmony(
  combined,
  group.by.vars = "sample_name"
)

# ---- Clustering and visualization ----
combined <- RunUMAP(
  combined,
  reduction = "harmony",
  dims = 1:30
)

combined <- FindNeighbors(
  combined,
  reduction = "harmony",
  dims = 1:30
)

combined <- FindClusters(
  combined,
  resolution = 0.1
)

# ---- Prepare reference and query assays ----
# Set RNA assay for reference object
DefaultAssay(RCC_IO) <- "RNA"

# Use sketched assay for query (Xenium)
DefaultAssay(combined) <- "sketch"

# Ensure assays are of class "Assay"
RCC_IO[["RNA"]] <- as(RCC_IO[["RNA"]], Class = "Assay")
combined[["sketch"]] <- as(combined[["sketch"]], Class = "Assay")

# ---- Build SingleCellExperiment objects ----
# Extract normalized expression matrices
ref_data <- LayerData(
  RCC_IO,
  assay = "RNA",
  layer = "data"
)

ref_sce <- SingleCellExperiment(
  assays = list(logcounts = ref_data)
)

# Store cell-type labels in reference SCE
ref_sce$label <- RCC_IO$subtypes

query_data <- LayerData(
  combined,
  assay = "sketch",
  layer = "data"
)

query_sce <- SingleCellExperiment(
  assays = list(logcounts = query_data)
)

# ---- Gene filtering ----
# Keep only genes shared between reference and query
common_genes <- intersect(
  rownames(ref_sce),
  rownames(query_sce)
)

# Stop if gene overlap is too small (e.g. Xenium panel too limited)
if (length(common_genes) < 200) {
  stop(
    "Too few shared genes: ",
    length(common_genes),
    ". Consider expanding the gene panel or using a different reference."
  )
}

ref_sce_f <- ref_sce[common_genes, , drop = FALSE]
query_sce_f <- query_sce[common_genes, , drop = FALSE]

# Remove genes rarely expressed in the query
keep <- rowMeans(
  assay(query_sce_f, "logcounts") > 0
) >= 0.01

ref_sce_f <- ref_sce_f[keep, , drop = FALSE]
query_sce_f <- query_sce_f[keep, , drop = FALSE]

# ---- Run SingleR annotation ----
singleR_pred <- SingleR(
  test = query_sce_f,
  ref = ref_sce_f,
  labels = ref_sce_f$label,
  assay.type.test = "logcounts",
  assay.type.ref  = "logcounts",
  genes = "de",           # Use DE genes for matching
  de.method = "wilcox",   # Differential expression test
  de.n = 200,             # Number of marker genes
  aggr.ref = TRUE,        # Aggregate reference by label
  quantile = 0.8,         # Similarity quantile
  fine.tune = TRUE,       # Fine-tune labels
  tune.thresh = 0.05,     # Threshold for fine-tuning
  prune = TRUE,           # Remove ambiguous assignments
  BPPARAM = BiocParallel::SerialParam()
)

# ---- Store SingleR results in Seurat metadata ----
sketch_cells <- colnames(query_sce_f)

combined$SingleR_labels_subtypes <- NA_character_
combined$SingleR_labels_subtypes[sketch_cells] <- singleR_pred$labels

# ---- Project annotations to full (non-sketched) data ----
combined <- ProjectData(
  object = combined,
  assay = "Xenium",
  full.reduction = "harmony.full",
  sketched.assay = "sketch",
  sketched.reduction = "harmony",
  umap.model = "umap",
  dims = 1:30,
  refdata = list(
    SingleR_labels_subtypes_full = "SingleR_labels_subtypes"
  )
)

# ---- Save final annotated object ----
saveRDS(
  combined,
  "/work/users/k/y/kyutae/Shared/combied_subtypes_all_1030.rds"
)




library(dplyr)
library(tibble)
library(ggplot2)
library(MASS)

# ============================================================
# SETTINGS
# ============================================================
obj <- combined
fov_name <- "post_ICI"
DefaultFOV(obj) <- fov_name

core_var <- "sample_name"
facet_var <- "Response"

# squidpy z-score matrix (wide): row = ...1 (source), col = target
source_type <- "C1QC-TAM"
target_type <- "Exhausted T cell(tEx)"

# KDE + contour
grid_n <- 220
bandwidth <- c(40, 40)
contour_probs <- c(0.95, 0.97, 0.99, 0.995, 0.998)

# packed layout
tile_margin <- 150
ncol_grid <- NULL

# colors
type_cols <- c(
  "C1QC-TAM" = "#0072B2",
  "Exhausted T cell(tEx)" = "#D55E00"
)

# ============================================================
# 0) Pick representative POST cores:
#    2 cores closest to group MEAN z-score per Response
# ============================================================
pair_z <- all_zscores %>%
  filter(...1 == source_type) %>%
  transmute(sample_name, z_score = .data[[target_type]])

resp_map <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  distinct(sample_name, Response)

pair_z2 <- pair_z %>%
  left_join(resp_map, by = "sample_name") %>%
  filter(grepl("^POST_", sample_name)) %>%
  filter(!is.na(Response), !is.na(z_score)) %>%
  mutate(Response = factor(Response, levels = c("Non-responder", "Responder")))

rep_post_median2_each <- pair_z2 %>%
  group_by(Response) %>%
  mutate(
    mean_z = mean(z_score, na.rm = TRUE),
    abs_from_mean = abs(z_score - mean_z)
  ) %>%
  arrange(abs_from_mean) %>%
  slice_head(n = 6) %>%
  ungroup() %>%
  dplyr::select(sample_name, Response, z_score, mean_z, abs_from_mean) %>%
  arrange(Response, abs_from_mean)

rep_post_median2_each


sel_cores <- rep_post_median2_each$sample_name
# sel_cores <- c("POST_T4",'POST_T6','POST_T15','POST_T14')
message("Selected representative cores (mean-based): ", paste(sel_cores, collapse = ", "))

# ============================================================
# 1) Coordinates + metadata (only selected cores, only two cell types)
# ============================================================
xy <- GetTissueCoordinates(obj, fov = fov_name)
if (is.matrix(xy)) xy <- as.data.frame(xy)

if (!"cell" %in% names(xy)) {
  xy <- xy %>% rownames_to_column("cell")
} else {
  rownames(xy) <- NULL
  xy <- as_tibble(xy)
}

if (!all(c("x", "y") %in% names(xy))) {
  if (all(c("imagecol", "imagerow") %in% names(xy))) {
    xy <- xy %>% rename(x = imagecol, y = imagerow)
  } else if (all(c("col", "row") %in% names(xy))) {
    xy <- xy %>% rename(x = col, y = row)
  } else {
    stop("Could not find x/y columns in xy: ", paste(names(xy), collapse = ", "))
  }
}
df <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  dplyr::select(cell, all_of(group_var), all_of(core_var), all_of(facet_var)) %>%
  dplyr::inner_join(xy %>% dplyr::select(cell, x, y), by = "cell") %>%  # <- 여기!
  dplyr::filter(.data[[core_var]] %in% sel_cores) %>%
  dplyr::filter(.data[[group_var]] %in% two_types) %>%
  dplyr::mutate(Response = factor(.data[[facet_var]], levels = c("Non-responder", "Responder")))

# ============================================================
# 2) Core-wise local coords + ranking within each Response panel
# ============================================================
core_info <- df %>%
  group_by(Response, .data[[core_var]]) %>%
  summarise(
    cx = median(x),
    cy = median(y),
    r  = quantile(sqrt((x - median(x))^2 + (y - median(y))^2), 0.98),
    n  = n(),
    .groups = "drop"
  ) %>%
  arrange(Response, desc(r), desc(n)) %>%
  group_by(Response) %>%
  mutate(core_rank = row_number()) %>%
  ungroup()

df2 <- df %>%
  dplyr::inner_join(
    core_info %>% dplyr::select(Response, all_of(core_var), cx, cy, r, core_rank),
    by = c("Response", core_var)
  ) %>%
  dplyr::mutate(
    x_local = x - cx,
    y_local = y - cy
  )

# ============================================================
# 3) Pack cores within each Response panel
# ============================================================
pack_by_panel <- function(dsub, info_sub) {
  k <- nrow(info_sub)
  ncol_use <- if (is.null(ncol_grid)) ceiling(sqrt(k)) else ncol_grid
  
  rmax <- max(info_sub$r, na.rm = TRUE)
  spacing <- 2 * rmax + tile_margin
  
  info_sub2 <- info_sub %>%
    mutate(
      grid_col = (core_rank - 1) %% ncol_use,
      grid_row = (core_rank - 1) %/% ncol_use
    )
  
  dsub %>%
    left_join(info_sub2 %>% dplyr::select(all_of(core_var), grid_col, grid_row), by = core_var) %>%
    mutate(
      x_pack = x_local + grid_col * spacing,
      y_pack = y_local + grid_row * spacing
    )
}

d_pack <- df2 %>%
  group_split(Response) %>%
  lapply(function(dsub) {
    resp <- unique(dsub$Response)
    info_sub <- core_info %>% filter(Response == resp)
    pack_by_panel(dsub, info_sub)
  }) %>%
  bind_rows()

core_labels <- d_pack %>%
  group_by(Response, sample_name) %>%
  summarise(
    x_lab = median(x_pack, na.rm = TRUE),
    y_lab = max(y_pack, na.rm = TRUE) + 120,  # 코어 위쪽에 살짝 띄움
    .groups = "drop"
  )


# ============================================================
# 4) Overlap KDE per Response panel (same bandwidth), global breaks
# ============================================================
make_overlap_grid <- function(dsub, n = 220) {
  n1 <- sum(dsub[[group_var]] == two_types[1])
  n2 <- sum(dsub[[group_var]] == two_types[2])
  if (n1 < 50 || n2 < 50) return(NULL)
  
  xlim <- range(dsub$x_pack, na.rm = TRUE)
  ylim <- range(dsub$y_pack, na.rm = TRUE)
  
  k1 <- MASS::kde2d(
    dsub$x_pack[dsub[[group_var]] == two_types[1]],
    dsub$y_pack[dsub[[group_var]] == two_types[1]],
    n = n, h = bandwidth, lims = c(xlim, ylim)
  )
  k2 <- MASS::kde2d(
    dsub$x_pack[dsub[[group_var]] == two_types[2]],
    dsub$y_pack[dsub[[group_var]] == two_types[2]],
    n = n, h = bandwidth, lims = c(xlim, ylim)
  )
  
  expand.grid(x = k1$x, y = k1$y) %>%
    mutate(score = log1p(pmin(as.vector(k1$z), as.vector(k2$z))))
}

grid_pack <- d_pack %>%
  group_split(Response) %>%
  lapply(function(dsub) {
    g <- make_overlap_grid(dsub, n = grid_n)
    if (is.null(g)) return(NULL)
    g$Response <- unique(dsub$Response)
    g
  }) %>%
  bind_rows()

z <- grid_pack$score[is.finite(grid_pack$score) & grid_pack$score > 0]
lev <- quantile(z, probs = contour_probs, na.rm = TRUE)
lev <- unique(as.numeric(lev[is.finite(lev)]))
lev <- lev[lev > 0]

# ============================================================
# 5) Plot
# ============================================================
p <- ggplot() +
  geom_point(
    data = d_pack,
    aes(x_pack, y_pack, color = .data[[group_var]]),
    size = 0.2, alpha = 0.4
  ) +
  geom_text(                      
    data = core_labels,
    aes(x = x_lab, y = y_lab, label = sample_name),
    size = 3,
    fontface = "bold"
  ) +
  geom_contour(
    data = grid_pack,
    aes(x, y, z = score),
    breaks = lev,
    color = "red",
    linewidth = 0.1,alpha = 0.8
  ) +
  coord_fixed(expand = FALSE) +
  facet_grid(.~Response) +
  scale_color_manual(values = type_cols) +
  theme_minimal(base_size = 12) +
  theme(
    panel.spacing = grid::unit(0.3, "lines"),
    legend.position = "right"
  ) +
  labs(
    title = paste0("Representative POST cores (2 mean-adjacent cores per Response): ",
                   two_types[1], " & ", two_types[2]),
    subtitle = paste0(
      "Cores retiled for visualization; KDE bandwidth=",
      paste(bandwidth, collapse = ", "),
      "; contour breaks are global (same across panels).",
      caption = "\nRepresentative POST cores: four cores per Response selected closest to the group mean neighborhood enrichment z-score (squidpy)."
      
    ),
    color = "type"
  )

print(p)


library(RANN)
library(dplyr)
library(ggplot2)

radius_um <- 15          # <- 원하는 반경(µm)
grid_step <- 10         # <- grid 해상도(µm). 5~20 사이 조절(작을수록 세밀/느림)

# Response 패널별로 proximity grid 계산
make_proximity_grid <- function(dsub, step = 10, r = 40) {
  # 전체 패널 범위에서 grid
  xg <- seq(min(dsub$x_pack), max(dsub$x_pack), by = step)
  yg <- seq(min(dsub$y_pack), max(dsub$y_pack), by = step)
  grid <- expand.grid(x = xg, y = yg)
  
  # 두 타입 좌표
  tam <- dsub %>% filter(.data[[group_var]] == two_types[1]) %>% dplyr::select(x_pack, y_pack)
  tex <- dsub %>% filter(.data[[group_var]] == two_types[2]) %>% dplyr::select(x_pack, y_pack)
  
  if (nrow(tam) < 5 || nrow(tex) < 5) return(NULL)
  
  # grid point -> nearest distance (k=1)
  nn_tam <- RANN::nn2(as.matrix(tam), as.matrix(grid), k = 1)
  nn_tex <- RANN::nn2(as.matrix(tex), as.matrix(grid), k = 1)
  
  grid %>%
    mutate(
      d_tam = nn_tam$nn.dists[, 1],
      d_tex = nn_tex$nn.dists[, 1],
      within = as.integer(d_tam <= r & d_tex <= r)
    )
}

prox_grid <- d_pack %>%
  group_split(Response) %>%
  lapply(function(dsub) {
    g <- make_proximity_grid(dsub, step = grid_step, r = radius_um)
    if (is.null(g)) return(NULL)
    g$Response <- unique(dsub$Response)
    g
  }) %>%
  bind_rows()

sel_text <- paste(sel_cores, collapse = ", ")

p40 <- ggplot() +
  geom_point(
    data = d_pack,
    aes(x_pack, y_pack, color = .data[[group_var]]),
    size = 0.18, alpha = 0.28
  ) +
  geom_contour(
    data = prox_grid,
    aes(x, y, z = within),
    breaks = 0.5,
    color = "red",
    linewidth = 0.15
  ) +
  geom_text(
    data = core_labels,
    aes(x = x_lab, y = y_lab, label = sample_name),
    size = 3,
    fontface = "bold"
  ) +
  coord_fixed(expand = FALSE) +
  facet_grid(Response ~. ) +
  scale_color_manual(values = type_cols) +
  theme_minimal(base_size = 12) +
  theme(
    panel.spacing = grid::unit(0.3, "lines"),
    legend.position = "right",
    plot.caption = element_text(hjust = 0, size = 9)
  ) +
  labs(
    title = paste0("15µm proximity contour (", two_types[1], " & ", two_types[2], ")"),
    subtitle = paste0("Contour: locations within ", radius_um, "µm of both cell types (grid step=", grid_step, "µm)."),
    caption = "Representative POST cores: six cores per Response selected closest to the group mean neighborhood enrichment z-score (squidpy)."
    ,
    color = "type"
  )
p40


