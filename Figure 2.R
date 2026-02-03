# Figure 2. Temporal PBMC analysis



#Figure 2 PBMC processing

library(Seurat)
library(harmony)
library(dplyr)
library(stringr)
library(Seurat)
library(dplyr)
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(dplyr)

PBMC_IO <- readRDS("/work/users/k/y/kyutae/final4_PBMC.rds")
PBMC_IO <- JoinLayers(PBMC_IO)
PBMC[['percent.mt']] <- PercentageFeatureSet(PBMC, pattern = '^MT-')
DefaultAssay(PBMC_IO) <- "RNA"  
sce <- scDblFinder(LayerData(PBMC_IO, assay = "RNA", layer = "counts"),  samples= PBMC_IO@meta.data$orig.ident)
PBMC_IO$scDblFinder.score <- sce$scDblFinder.score
PBMC_IO$scDblFinder.class <- sce$scDblFinder.class 
VlnPlot(PBMC_IO,c("nFeature_RNA",'percent.mt','scDblFinder.score'),pt.size = 0)

PBMC_IO <- subset(PBMC_IO,subset = nFeature_RNA > 300 &
                 nFeature_RNA < 6000 &
                 percent.mt < 10  & 
                 scDblFinder.class == 'singlet')
VlnPlot(PBMC_IO,c("nFeature_RNA",'percent.mt','scDblFinder.score'),pt.size = 0)

PBMC.list <- SplitObject(PBMC_IO, split.by = "orig.ident")
PBMC.list <- lapply(PBMC.list, NormalizeData)
PBMC.list <- lapply(PBMC.list, FindVariableFeatures, selection.method = "vst", nfeatures = 2000)
features <- SelectIntegrationFeatures(PBMC.list, nfeatures = 3000)  
anchors <- FindIntegrationAnchors(
  object.list = PBMC.list, 
  anchor.features = features, 
  reduction = "cca",  
  dims = 1:15
)
PBMC_IO <- IntegrateData(anchorset = anchors, dims = 1:15)
PBMC_IO <- ScaleData(PBMC_IO)
PBMC_IO <- RunPCA(PBMC_IO, npcs = 50)
ElbowPlot(PBMC,ndims = 50)
PBMC_IO <- FindNeighbors(PBMC_IO, dims = 1:15, graph.name = 'mygraph')
PBMC_IO <- FindClusters(PBMC_IO, resolution = 1.2, graph.name = 'mygraph',random.seed = 42)
PBMC_IO <- RunUMAP(PBMC_IO, dims = 1:15,dim = 1:15,seed.use = 42)
saveRDS(PBMC_IO,"/work/users/k/y/kyutae/IO_PBMC_splited.rds")
DefaultAssay(PBMC_IO) <- "RNA"  
PBMC_IO <- JoinLayers(PBMC_IO)


#Supplemental figure 1B
DotPlot(PBMC_IO,group.by = 'seurat_clusters', features = c('CD3D','CD3E','CD8A','CD8B','CD4','CD40LG','NCAM1','GNLY','MKI67','TOP2A','CD14','CD163','LYZ','S100A12','S100A9','S100A8','FCGR3A','IGKC','IGHM','FCER1A','CD1E','CD1D','CD79A','MS4A1','BANK1','SELP','LILRA4','TCF4','TPSAB1','KIT'))+ RotatedAxis() + scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") 
#Supplemental figure 1A
PBMC_IO@meta.data <- PBMC_IO@meta.data %>%
  mutate(cell_type = case_when(
    seurat_clusters == "0" ~ "C_Monocyte", 
    seurat_clusters == "1" ~ "CD4+ T cell", 
    seurat_clusters == "2" ~ "NK cell", 
    seurat_clusters == "3" ~ "C_Monocyte", 
    seurat_clusters == "4" ~ "CD8+ T cell", 
    seurat_clusters == "5" ~ "CD4+ T cell",
    seurat_clusters == "6" ~ "NK cell",
    seurat_clusters == "7" ~ "NK cell",
    seurat_clusters == "8" ~ "NC_Monocyte",
    seurat_clusters == "9" ~ "CD8+ T cell",
    seurat_clusters == "10" ~ "B cell",
    seurat_clusters == '11' ~ 'CD8+ T cell',
    seurat_clusters == '12' ~ 'C_Monocyte',
    seurat_clusters == '13' ~ 'CD8+ T cell',
    seurat_clusters == '14' ~ 'CD4+ T cell',
    seurat_clusters == '15' ~ 'Proliferative_T',
    seurat_clusters == '16' ~ 'C_Monocyte',
    seurat_clusters == '17' ~ 'Platelet',
    seurat_clusters == '18' ~ 'Dendritic cell',
    seurat_clusters == '19' ~ 'Plasma cell',
    seurat_clusters == '20' ~ 'Mast cell',
    seurat_clusters == '21' ~ 'C_Monocyte',
    TRUE ~ NA_character_
  ))

plot <- DimPlot(
  object = PBMC_IO,
  group.by = "cell_type",
  reduction = "umap",
  label = F,  
  label.size = 6, 
  pt.size = 1.2 
  ,raster = F)

library(RColorBrewer )

set3_colors <- c(brewer.pal(n = 12, name = "Set3"), 
                 '#D2691E',  
                 '#B39DDB',  
                 '#FFAB91',  
                 '#4DB6AC',
                 '#F48FB1', 
                 '#C5E1A5',
                 '#FFD54F'  
)


plot + 
  scale_color_manual(values = set3_colors) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  
    axis.title = element_text(size = 14, face = "bold"),  
    axis.text = element_text(size = 12),  
    legend.title = element_text(size = 14, face = "bold"),  
    legend.text = element_text(size = 12), 
    legend.key.size = unit(1.5, "lines"),  
    legend.position = "right"  
  ) +
  labs(
    title = "Tissue-T/NK cells",  
    x = "UMAP1",  
    y = "UMAP2",  
    color = "Cell types" 
  ) 




#PBMC_global cell type proportional analysis
sample_counts <- PBMC_IO@meta.data %>%
  filter(Timepoint != 'T3') %>% 
  group_by(sample_name) %>%
  summarise(total_cells = n(), .groups = "drop")

cell_type_counts <- PBMC_IO@meta.data %>%
  group_by(sample_name, cell_type) %>%
  summarise(cell_count = n(), .groups = "drop")

cell_type_proportions <- cell_type_counts %>%
  left_join(sample_counts, by = "sample_name") %>%
  mutate(proportion = cell_count / total_cells) %>%
  left_join(
    PBMC_IO@meta.data %>%
      select(sample_name, Timepoint, Response, patients) %>%
      distinct(),
    by = "sample_name"
  )

df_all <- cell_type_proportions %>%
  filter(
    !is.na(Timepoint),
    !is.na(proportion),
    !is.na(Response)
  ) %>%
  droplevels()

pd <- position_dodge(width = 0.6)

library(dplyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(RColorBrewer)
library(patchwork)   

df_time_resp <- df_all %>%
  filter(!is.na(Response),
         !is.na(Timepoint),
         !is.na(proportion)) %>%
  droplevels()

df_time_resp <- df_time_resp %>%
  mutate(
    Timepoint = factor(Timepoint),
    Response  = factor(Response),
    cell_type = factor(cell_type)
  )

x_ref <- levels(df_time_resp$Timepoint)[1]  
panel_cap_resp <- df_time_resp %>%
  group_by(cell_type, Response) %>%
  summarise(
    y_min = min(proportion, na.rm = TRUE),
    y_max = max(proportion, na.rm = TRUE),
    .groups = "drop"
  )

pw_resp <- df_time_resp %>%
  group_by(cell_type, Response) %>%
  filter(n_distinct(Timepoint) >= 2) %>%
  wilcox_test(
    proportion ~ Timepoint,
    p.adjust.method = "BH"
  ) %>%
  ungroup() %>%
  # filter(!( (group1 == "T0" & group2 == "T2") | (group1 == "T2" & group2 == "T0") )) %>%
  add_x_position(x = "Timepoint") %>%
  left_join(panel_cap_resp, by = c("cell_type", "Response"))

alpha_cutoff <- 1  

pw_resp_plot <- pw_resp %>%
  filter(!is.na(p.adj), p.adj < alpha_cutoff) %>%  
  group_by(cell_type, Response) %>%
  arrange(p.adj, .by_group = TRUE) %>%            
  mutate(
    y_range = y_max - y_min,
    y_step  = ifelse(y_range == 0, 0.02, 0.10 * y_range),  
    y.position = y_max + y_step * row_number(),
    label = paste0("adj p=", rstatix::p_format(p.adj, accuracy = 0.001))
  ) %>%
  ungroup()

time_levels_resp <- levels(df_time_resp$Timepoint)
time_colors_resp <- brewer.pal(n = max(3, length(time_levels_resp)), "Set2")[seq_along(time_levels_resp)]
names(time_colors_resp) <- time_levels_resp

make_resp_plot <- function(resp_level, show_x = TRUE) {
  
  df_sub <- df_time_resp %>% filter(Response == resp_level)
  pw_sub <- pw_resp_plot  %>% filter(Response == resp_level)
  
  panel_blank <- pw_sub %>%
    group_by(cell_type) %>%
    summarise(
      y_blank = ifelse(all(is.na(y.position)), NA_real_, max(y.position, na.rm = TRUE) * 1.03),
      .groups = "drop"
    ) %>%
    mutate(Timepoint = x_ref)
  
  ggplot(df_sub, aes(x = Timepoint, y = proportion, fill = Timepoint)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, color = "black") +
    geom_jitter(
      aes(color = Timepoint),
      position = position_jitter(width = 0.1, height = 0),
      size = 1.2, alpha = 0.8, show.legend = FALSE
    ) +
    geom_blank(data = panel_blank, aes(x = Timepoint, y = y_blank), inherit.aes = FALSE) +
    ggpubr::stat_pvalue_manual(
      data = pw_sub,
      label = "label",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.5,
      size = 3.2,
      na.rm = TRUE
    ) +
    facet_wrap(~ cell_type, scales = "free_y", nrow = 1) +  
    scale_fill_manual(values = time_colors_resp) +
    scale_color_manual(values = time_colors_resp) +
    theme_bw() +
    theme(
      text = element_text(size = 12),
      strip.text = element_text(size = 10),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 13, face = "bold"),
      legend.position = "none"
    ) +
    labs(
      title = resp_level,
      x = if (show_x) "Timepoint" else NULL,
      y = "Proportion",
      fill = "Timepoint"
    ) +
    {if(!show_x) theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) else theme()}
}

p_non  <- make_resp_plot("Non_Responder", show_x = FALSE)
p_resp <- make_resp_plot("Responder",     show_x = TRUE)
#Figure 2A
p_time_by_response_2row <- (p_non / p_resp) +
  plot_layout(heights = c(1, 1)) &
  theme(legend.position = "bottom")

p_time_by_response_2row


PBMC_T <- subset(PBMC_IO, cell_type %in% c('CD4+ T cell','CD8+ T cell','NK cell','Proliferative_T'))
DefaultAssay(PBMC_T) <- "integrated"
PBMC_B@assays$integrated
PBMC_T <- RunPCA(PBMC_T, features = VariableFeatures(object= PBMC_T),npcs = 50)
print("RunPCA DONE")
ElbowPlot(PBMC_T,ndims = 50)
ElbowPlot(PBMC_T,ndims = 50)
PBMC_T <- FindNeighbors(PBMC_T, dims= 1:10,reduction = 'pca', graph.name = 'mygraph')
PBMC_T <- FindClusters(PBMC_T, resolution = 1, graph.name = 'mygraph',random.seed = 42)
PBMC_T <- RunUMAP(PBMC_T,dims = 1:10,seed.use = 42,reduction = 'pca')
DotPlot(PBMC_T, group.by = 'seurat_clusters', 
        features = c('CD3D','CD3E','CD4','CD40LG','CD8A','CD8B','NCAM1','MKI67','TOP2A','CD44','GZMB','PRF1','GZMK','IFNG','TRDV2','TRGV9','TRGV10','SLC4A10','TRAV1-2','KLRD1','FGFBP2','CX3CR1','KLRG1','HBB','FOXP3','CTLA4','SELL','IL7R','GPR183','CD69','SLAMF6','PDCD1','TCF7','TOX','HAVCR2')) + 
  ggtitle("PBMC T/NK cell Markers") + 
  RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
    geom_vline(xintercept = c(9.5, 16.5), linetype = "dashed", color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

PBMC_T@meta.data <- PBMC_T@meta.data %>%
  mutate(T_NK_cell = case_when(
    seurat_clusters == "0" ~ 'NKT',
    seurat_clusters == "1" ~ "CD4_Naive",
    seurat_clusters == "2" ~ "CD56-low_NK",
    seurat_clusters == "3" ~ "CD8_Temra" ,
    seurat_clusters == "4" ~ "CD4_Naive",
    seurat_clusters == "5" ~ "CD8_EFF",
    seurat_clusters == "6" ~ "CD56-high_NK",
    seurat_clusters == "7" ~ "CD4_Treg",
    seurat_clusters == "8" ~ "CD4_CM",
    seurat_clusters == "9" ~ "Tgd",
    seurat_clusters == "10" ~ "CCD4_Naive",
    seurat_clusters == "11" ~ "CD8_Naive",
    seurat_clusters == "12" ~ "CD56-low_NK",
    seurat_clusters == "13" ~ "CD56-low_NK",
    seurat_clusters == "14" ~ "CD8_MAIT",
    seurat_clusters == "15" ~ "Proliferative",
    seurat_clusters == "16" ~ "CD8_pEx",
    seurat_clusters == "17" ~ "Erythrocyte",
    TRUE ~ NA_character_
  ))


#PBMC T/NK cell proportional analysis
sample_counts <- PBMC_T@meta.data %>%
  filter(Timepoint != 'T3') %>% 
  group_by(sample_name) %>%
  summarise(total_cells = n(), .groups = "drop")

cell_type_counts <- PBMC_T@meta.data %>%
  group_by(sample_name, T_NK_cell) %>%
  summarise(cell_count = n(), .groups = "drop")

cell_type_proportions <- cell_type_counts %>%
  left_join(sample_counts, by = "sample_name") %>%
  mutate(proportion = cell_count / total_cells) %>%
  left_join(
    PBMC_IO@meta.data %>%
      select(sample_name, Timepoint, Response, patients) %>%
      distinct(),
    by = "sample_name"
  )

df_all <- cell_type_proportions %>%
  filter(
    !is.na(Timepoint),
    !is.na(proportion),
    !is.na(Response)
  ) %>%
  droplevels()

pd <- position_dodge(width = 0.6)

library(dplyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(RColorBrewer)
library(patchwork)   

df_time_resp <- df_all %>%
  filter(!is.na(Response),
         !is.na(Timepoint),
         !is.na(proportion)) %>%
  droplevels()

df_time_resp <- df_time_resp %>%
  mutate(
    Timepoint = factor(Timepoint),
    Response  = factor(Response),
    T_NK_cell = factor(T_NK_cell)
  )

x_ref <- levels(df_time_resp$Timepoint)[1]  
panel_cap_resp <- df_time_resp %>%
  group_by(T_NK_cell, Response) %>%
  summarise(
    y_min = min(proportion, na.rm = TRUE),
    y_max = max(proportion, na.rm = TRUE),
    .groups = "drop"
  )

pw_resp <- df_time_resp %>%
  group_by(T_NK_cell, Response) %>%
  filter(n_distinct(Timepoint) >= 2) %>%
  wilcox_test(
    proportion ~ Timepoint,
    p.adjust.method = "BH"
  ) %>%
  ungroup() %>%
  # filter(!( (group1 == "T0" & group2 == "T2") | (group1 == "T2" & group2 == "T0") )) %>%
  add_x_position(x = "Timepoint") %>%
  left_join(panel_cap_resp, by = c("cell_type", "Response"))

alpha_cutoff <- 1  

pw_resp_plot <- pw_resp %>%
  filter(!is.na(p.adj), p.adj < alpha_cutoff) %>%  
  group_by(T_NK_cell, Response) %>%
  arrange(p.adj, .by_group = TRUE) %>%            
  mutate(
    y_range = y_max - y_min,
    y_step  = ifelse(y_range == 0, 0.02, 0.10 * y_range),  
    y.position = y_max + y_step * row_number(),
    label = paste0("adj p=", rstatix::p_format(p.adj, accuracy = 0.001))
  ) %>%
  ungroup()

time_levels_resp <- levels(df_time_resp$Timepoint)
time_colors_resp <- brewer.pal(n = max(3, length(time_levels_resp)), "Set2")[seq_along(time_levels_resp)]
names(time_colors_resp) <- time_levels_resp

make_resp_plot <- function(resp_level, show_x = TRUE) {
  
  df_sub <- df_time_resp %>% filter(Response == resp_level)
  pw_sub <- pw_resp_plot  %>% filter(Response == resp_level)
  
  panel_blank <- pw_sub %>%
    group_by(cell_type) %>%
    summarise(
      y_blank = ifelse(all(is.na(y.position)), NA_real_, max(y.position, na.rm = TRUE) * 1.03),
      .groups = "drop"
    ) %>%
    mutate(Timepoint = x_ref)
  
  ggplot(df_sub, aes(x = Timepoint, y = proportion, fill = Timepoint)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, color = "black") +
    geom_jitter(
      aes(color = Timepoint),
      position = position_jitter(width = 0.1, height = 0),
      size = 1.2, alpha = 0.8, show.legend = FALSE
    ) +
    geom_blank(data = panel_blank, aes(x = Timepoint, y = y_blank), inherit.aes = FALSE) +
    ggpubr::stat_pvalue_manual(
      data = pw_sub,
      label = "label",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.5,
      size = 3.2,
      na.rm = TRUE
    ) +
    facet_wrap(~ T_NK_cell, scales = "free_y", nrow = 1) + 
    scale_fill_manual(values = time_colors_resp) +
    scale_color_manual(values = time_colors_resp) +
    theme_bw() +
    theme(
      text = element_text(size = 12),
      strip.text = element_text(size = 10),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 13, face = "bold"),
      legend.position = "none"
    ) +
    labs(
      title = resp_level,
      x = if (show_x) "Timepoint" else NULL,
      y = "Proportion",
      fill = "Timepoint"
    ) +
    {if(!show_x) theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) else theme()}
}

p_non  <- make_resp_plot("Non_Responder", show_x = FALSE)
p_resp <- make_resp_plot("Responder",     show_x = TRUE)
#Figure 2C
p_time_by_response_2row <- (p_non / p_resp) +
  plot_layout(heights = c(1, 1)) &
  theme(legend.position = "bottom")

p_time_by_response_2row

PBMC_T@meta.data$T_R <- paste0(PBMC_T@meta.data$Timepoint,'_',PBMC_T@meta.data$Response)
#Figure2D
dittoSeq::dittoBarPlot(
  object = PBMC_T,
  var = "T_NK_cell",
  group.by = "T_R",
  ylab = "Proportion",  
  color.panel = set3_colors,
  theme = theme_minimal(base_size = 14) )+ 
  # scale_fill_manual(values = set3_colors,
  #                   breaks = levels(T_cell_T1@meta.data$TCF7_PDCD1))+
  theme(
    axis.text = element_text(size = 12, family = "serif"),
    axis.title = element_text(size = 11, face = "bold", family = "serif"),
    legend.title = element_text(size = 12, family = "serif"),
    legend.text = element_text(size = 10, family = "serif"),
    plot.title = element_text(size = 12, family = "serif")
  ) 


C15 <- subset(PBMC_T, T_NK_cell == 'Proliferative_T')
DefaultAssay(C15) <- 'integrated'
C15 <- RunUMAP(C15, dims = 1:10)
C15 <- FindNeighbors(C15, dims = 1:10)
C15 <- FindClusters(C15, resolution = 1)
DefaultAssay(C15) <- 'RNA'
DotPlot(C15,group.by = 'seurat_clusters',c('NCAM1','CD3D','CD8A','CD4'),scale = T)

C15@meta.data <- C15@meta.data %>%
  mutate(Pro_types = case_when(
    seurat_clusters == "0" ~ "CD8", 
    seurat_clusters == "1" ~ "NK", 
    seurat_clusters == "2" ~ "CD8", 
    seurat_clusters == "3" ~ "NK", 
    seurat_clusters == "4" ~ "Treg", 
    seurat_clusters == "5" ~ "Treg",
    seurat_clusters == "6" ~ "NK",
    seurat_clusters == "7" ~ "NK",
    seurat_clusters == "8" ~ "CD8",
    seurat_clusters == "9" ~ "CD4",
    seurat_clusters == "10" ~ "CD8",
    seurat_clusters == '11' ~ 'CD4CD8',
    seurat_clusters == '12' ~ 'CD4CD8',
    TRUE ~ NA_character_
  ))

  #Figure 2E
dittoSeq::dittoBarPlot(C15,var = 'Pro_types',group.by = 'T_R')

C15_T0 <- subset(PBMC_T,T_NK_cell == 'Proliferative' & Timepoint == 'T0')
C15_T1 <- subset(PBMC_T,T_NK_cell == 'Proliferative' & Timepoint == 'T1')
C15_T2 <- subset(PBMC_T,T_NK_cell == 'Proliferative' & Timepoint == 'T2')

C15_DE <- FindMarkers(C15_T1,ident.1 = 'Responder',ident.2 = 'Non_Responder',group.by = 'Response')
C15_DE_2 <- FindMarkers(C15_T2,ident.1 = 'Responder',ident.2 = 'Non_Responder',group.by = 'Response')
C15_DE_3 <- FindMarkers(C15_T2,ident.1 = 'Responder',ident.2 = 'Non_Responder',group.by = 'Response')

C15_DE <- filter(C15_DE,pct.1 > 0.1 & pct.2 > 0.1)
C15_DE_2 <- filter(C15_DE_2,pct.1 > 0.1 & pct.2 > 0.1)
C15_DE_3 <- filter(C15_DE_3,pct.1 > 0.1 & pct.2 > 0.1)

EnhancedVolcano(C15_DE,x = 'avg_log2FC',y = 'p_val_adj',lab = rownames(C15_DE),pCutoff = 0.01,FCcutoff = 0.5) + 
  EnhancedVolcano(C15_DE_2,x = 'avg_log2FC',y = 'p_val_adj',lab = rownames(C15_DE_2),pCutoff = 0.01,FCcutoff = 0.5,selectLab = c("NKG7",'LAG3','FGFBP2','PRF1','CCL4','GZMB','HAVCR2','CST7','GZMA','PRDM1','CX3CR1','NFKBIA','IFI6','MAP3K8','IL2RB','MX1','GZMH','CD52','LEF1','TCF7','IL7R'),max.overlaps = Inf,drawConnectors = T) +
   EnhancedVolcano(C15_DE_3,x = 'avg_log2FC',y = 'p_val_adj',lab = rownames(C15_DE_3),pCutoff = 0.01,FCcutoff = 0.5,selectLab = c('CX3CR1','HLA-B','CXCR4','CTLA4','FABP5'),max.overlaps = Inf,drawConnectors = T) 





library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(rlang)

PBMC_T <- AddModuleScore(
  PBMC_T,
  features = as.list(T_cell), #T cell atlas signature
  name     = names(T_cell)  
)

meta0 <- PBMC_T@meta.data

## =========================
## 1) 메타데이터 정리
## =========================

df <- meta0 %>%
  mutate(
    # T_R: 예) "T0_Responder", "T1_PD_at_T4" 이런 형식이라고 가정
    TimePoint = sub("_.*", "", T_R),      # "T0","T1","T2","T3",...
    subtype   = as.character(Merged),
    
    # 🔴 ResponseGroup: Responder vs Non_Responder만 사용, PD_at_T4는 제외
    ResponseGroup = case_when(
      Response == "Responder"      ~ "Responder",
      Response == "Non_Responder"  ~ "Non_Responder",
      TRUE                         ~ NA_character_       # PD_at_T4 포함해서 NA 처리
    )
  ) %>%
  # T0, T1, T2만 사용 (필요하면 T3 추가 가능)
  filter(TimePoint %in% c("T0", "T1", "T2")) %>%
  # CD8 / Tgd subset만 보고 싶으면 (원래 코드 그대로)
  filter(grepl("CD8|Tgd", subtype)) %>%
  # PD_at_T4 등 ResponseGroup이 NA인 애들은 비교에서 제외
  filter(!is.na(ResponseGroup))

# 🔴 시그니처 컬럼: meta에서 위치로 가져온 뒤, 실제 있는 것만 사용
sig_cols <- colnames(meta0)[39:57]
sig_cols <- sig_cols[sig_cols %in% colnames(df)]
if (length(sig_cols) == 0) {
  stop("sig_cols 범위(39:57)에 시그니처 컬럼이 없습니다. 위치를 다시 확인하세요.")
}

## =========================
## 2) 한 feature × timepoint × subtype의 log2FC 계산 함수
## =========================

calc_log2fc_one <- function(dat, feature, tp, subtype) {
  sub_df <- dat %>%
    filter(TimePoint == tp, subtype == !!subtype) %>%
    mutate(val = .data[[feature]]) %>%
    select(ResponseGroup, val) %>%
    filter(!is.na(val), !is.na(ResponseGroup))
  
  # 해당 조건에 셀 없거나, 그룹이 하나 뿐이면 NA
  if (nrow(sub_df) == 0 || dplyr::n_distinct(sub_df$ResponseGroup) < 2) {
    return(tibble::tibble(
      feature   = feature,
      timepoint = tp,
      subtype   = subtype,
      log2FC    = NA_real_,
      p         = NA_real_
    ))
  }
  
  # 🔑 음수/0 방지: PBMC_T 코드에서처럼 shift + eps
  eps   <- 1e-6
  shift <- -min(sub_df$val, na.rm = TRUE)
  sub_df <- mutate(sub_df, val_pos = val + shift + eps)
  
  mean_resp <- mean(sub_df$val_pos[sub_df$ResponseGroup == "Responder"],     na.rm = TRUE)
  mean_non  <- mean(sub_df$val_pos[sub_df$ResponseGroup == "Non_Responder"], na.rm = TRUE)
  
  log2fc <- log2(mean_resp / mean_non)
  
  pval <- tryCatch(
    stats::wilcox.test(val ~ ResponseGroup, data = sub_df)$p.value,
    error = function(e) NA_real_
  )
  
  tibble::tibble(
    feature   = feature,
    timepoint = tp,
    subtype   = subtype,
    log2FC    = log2fc,
    p         = pval
  )
}

## =========================
## 3) 전체 feature × timepoint × subtype에 대해 계산
## =========================

subtypes   <- sort(unique(df$subtype))
timepoints <- c("T0", "T1", "T2")

res_all <- purrr::map_dfr(sig_cols, function(feat) {
  purrr::map_dfr(timepoints, function(tp) {
    purrr::map_dfr(subtypes, ~ calc_log2fc_one(df, feat, tp, .x))
  })
})

res_all <- as_tibble(res_all)

if (is.list(res_all$p)) {
  res_all$p <- sapply(res_all$p, function(x) {
    # x가 numeric(1) 이거나 길이 1짜리 리스트라고 가정
    if (length(x) == 0 || all(is.na(x))) {
      return(NA_real_)
    } else {
      return(as.numeric(x)[1])
    }
  })
}

res_all <- res_all %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    signif = case_when(
      is.na(p_adj) ~ "",
      p_adj < 1e-4 ~ "****",
      p_adj < 1e-3 ~ "***",
      p_adj < 1e-2 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    ),
    col_key = paste0(timepoint, "_", subtype)   # 예: "T0_CD8_EFF"
  )
col_order <- c(
  grep("^T0_", unique(res_all$col_key), value = TRUE),
  grep("^T1_", unique(res_all$col_key), value = TRUE),
  grep("^T2_", unique(res_all$col_key), value = TRUE)
)

log2fc_mat <- res_all %>%
  select(feature, col_key, log2FC) %>%
  distinct(feature, col_key, .keep_all = TRUE) %>%
  pivot_wider(names_from = col_key, values_from = log2FC) %>%
  tibble::column_to_rownames("feature") %>%
  as.matrix()

star_mat <- res_all %>%
  select(feature, col_key, signif) %>%
  distinct(feature, col_key, .keep_all = TRUE) %>%
  pivot_wider(names_from = col_key, values_from = signif) %>%
  tibble::column_to_rownames("feature") %>%
  as.matrix()

col_order <- intersect(col_order, colnames(log2fc_mat))
log2fc_mat <- log2fc_mat[, col_order, drop = FALSE]
star_mat   <- star_mat[,   col_order, drop = FALSE]

## =========================
## 5) Heatmap (Timepoint로 column split)
## =========================

tp_levels <- c("T0", "T1", "T2")
col_tp <- factor(sub("_.*", "", colnames(log2fc_mat)), levels = tp_levels)

col_fun <- circlize::colorRamp2(
  c(-1.5, 0, 1.5),
  c("#2166ac", "white", "#b2182b")
)

grid::grid.newpage()
ht <- ComplexHeatmap::Heatmap(
  log2fc_mat,
  name = "log2FC\n(Resp / NonResp)",
  col  = col_fun,
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  column_split    = col_tp,             
  gap             = unit(5, "mm"),
  column_title    = "Signature log2FC by timepoint & CD8/Tgd subtypes",
  column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
  column_names_rot = 45,
  row_names_gp     = grid::gpar(fontsize = 11),
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontface = "bold"),
    at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
  ),
  cell_fun = function(j, i, x, y, w, h, fill) {
    lab <- star_mat[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      grid::grid.text(lab, x, y, gp = grid::gpar(fontsize = 10, fontface = "bold"))
    }
  },
  border = TRUE
)

ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "bottom")






















