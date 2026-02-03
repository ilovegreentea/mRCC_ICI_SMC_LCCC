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




#CD8 functional heatmap

library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(rlang)
T_cell_sig <- read_excel("T cell.xlsx")
#T cell atlas CD8 functional signature
#Chu, Y., Dai, E., Li, Y., Han, G., Pei, G., Ingram, D. R., ... & Wang, L. (2023). Pan-cancer T cell atlas links a cellular stress response state to immunotherapy resistance. Nature medicine, 29(6), 1550-1562.
PBMC_T <- AddModuleScore(PBMC_T,features = as.list(T_cell_sig),name = names(T_cell_sig))

df <- meta0 %>%
  mutate(
    TimePoint     = sub("_.*", "", T_R),
    ResponseGroup = if_else(grepl("Non_Responder", T_R), "Non_Responder", "Responder"),
    T_NK_cell        = as.character(T_NK_cell)
  ) %>%
  filter(
    TimePoint %in% c("T0","T1","T2"),
    grepl("CD8", T_NK_cell)
  )

sig_cols  <- colnames(meta0)[c(39:57)]
subtypes  <- sort(unique(df$T_NK_cell))
timepoints <- c("T0","T1","T2")

calc_meandiff_one <- function(dat, feature, tp, subtype) {
  
  sub_df <- dat %>%
    filter(TimePoint == tp, T_NK_cell == subtype) %>%
    mutate(val = .data[[feature]]) %>%
    select(ResponseGroup, val) %>%
    filter(!is.na(val), !is.na(ResponseGroup))
  
  # 데이터 부족할 때
  if (nrow(sub_df) == 0 || n_distinct(sub_df$ResponseGroup) < 2) {
    return(tibble(
      feature   = feature,
      timepoint = tp,
      subtype   = subtype,
      mean_diff = NA_real_,
      p         = NA_real_
    ))
  }
  
  mean_resp <- mean(sub_df$val[sub_df$ResponseGroup == "Responder"])
  mean_non  <- mean(sub_df$val[sub_df$ResponseGroup == "Non_Responder"])
  
  mean_diff <- mean_resp - mean_non   
  
  pval <- tryCatch(
    wilcox.test(val ~ ResponseGroup, data = sub_df)$p.value,
    error = function(e) NA_real_
  )
  
  tibble(
    feature   = feature,
    timepoint = tp,
    subtype   = subtype,
    mean_diff = mean_diff,
    p         = pval
  )
}

res_all <- map_dfr(sig_cols, function(feat) {
  map_dfr(timepoints, function(tp) {
    map_dfr(subtypes, ~ calc_meandiff_one(df, feat, tp, .x))
  })
})

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
    col_key = paste0(timepoint, "_", subtype)
  )

col_order <- c(
  grep("^T0_", unique(res_all$col_key), value = TRUE),
  grep("^T1_", unique(res_all$col_key), value = TRUE),
  grep("^T2_", unique(res_all$col_key), value = TRUE)
)

mean_diff_mat <- res_all %>%
  select(feature, col_key, mean_diff) %>%
  distinct(feature, col_key, .keep_all = TRUE) %>%
  pivot_wider(names_from = col_key, values_from = mean_diff) %>%
  column_to_rownames("feature") %>%
  as.matrix()

star_mat <- res_all %>%
  select(feature, col_key, signif) %>%
  distinct(feature, col_key, .keep_all = TRUE) %>%
  pivot_wider(names_from = col_key, values_from = signif) %>%
  column_to_rownames("feature") %>%
  as.matrix()

col_order <- intersect(col_order, colnames(mean_diff_mat))
mean_diff_mat <- mean_diff_mat[, col_order, drop = FALSE]
star_mat      <- star_mat[, col_order, drop = FALSE]

col_fun <- colorRamp2(c(-0.5, 0, 0.5),
                      c("#2166ac", "white", "#b2182b"))

tp_levels <- c("T0","T1","T2")
col_tp <- factor(sub("_.*","",colnames(mean_diff_mat)), levels = tp_levels)

grid.newpage()
ht <- Heatmap(
  mean_diff_mat,
  name  = "Δmean (R-NR)",
  col   = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = col_tp,
  gap = unit(5, "mm"),
  column_title = "Signature mean difference (Responder - Non_Responder)",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  column_names_rot = 45,
  heatmap_legend_param = list(
    title_gp = gpar(fontface = "bold"),
    at = c(-0.5, 0, 0.5)
  ),
  cell_fun = function(j, i, x, y, w, h, fill) {
    lab <- star_mat[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      grid.text(lab, x, y, gp = gpar(fontsize = 10, fontface = "bold"))
    }
  },
  border = TRUE
)
draw(ht)

#Lineplot Dynamic change across timepoints
library(dplyr)
library(tidyr)
library(ggplot2)
PBMC_CD8 <- subset(PBMC_T, Merged %in% c("CD8_EFF","CD8_pEx","CD8_Temra","CD8_MAIT","CD8_Naive"))

PBMC_CD8 <-
  AddModuleScore(
    object = PBMC_CD8,
    features = hallmark_list,
    name = names(hallmark_list),
    assay = 'RNA') 

sig_cols <- c(
  "IFN Response14",
  "MAPK Signaling12",
  "NFKB Signaling10"
)

cell_types <- c("CD8_EFF", "CD8_MAIT", "CD8_Temra", "CD8_Naive", "CD8_pEx")

sig_df <- PBMC_CD8@meta.data %>%
  filter(
    T_NK_cell %in% cell_types,
    Timepoint != "T3"
  ) %>%
  mutate(
    TimePoint = sub("_.*", "", T_R)
  ) %>%
  filter(TimePoint %in% c("T0", "T1", "T2")) %>%
  select(T_NK_cell, Response, TimePoint, all_of(sig_cols)) %>%
  pivot_longer(
    cols      = all_of(sig_cols),
    names_to  = "Signature",
    values_to = "value"
  ) %>%
  drop_na(value)

sig_df$TimePoint <- factor(sig_df$TimePoint, levels = c("T0", "T1", "T2"))

## 4) mean / sd / n → mean difference + 95% CI
summary_ci <- sig_df %>%
  group_by(Signature, T_NK_cell, TimePoint, Response) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    n    = dplyr::n(),
    .groups = "drop"
  ) %>%
  # Responder / Non_Responder를 wide로
  pivot_wider(
    names_from  = Response,
    values_from = c(mean, sd, n)
  ) %>%
  # mean diff & SE & 95% CI
  mutate(
    mean_diff_R_NR = mean_Responder - mean_Non_Responder,
    se_diff = sqrt(
      (sd_Responder^2 / n_Responder) +
        (sd_Non_Responder^2 / n_Non_Responder)
    ),
    ci_low  = mean_diff_R_NR - 1.96 * se_diff,
    ci_high = mean_diff_R_NR + 1.96 * se_diff
  )
ggplot(summary_ci,
       aes(x = TimePoint, y = mean_diff_R_NR,
           group = T_NK_cell, color = T_NK_cell)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.15,
    linewidth = 0.4
  ) +
  facet_wrap(~ Signature, scales = "free_y") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Signature mean difference (Responder − Non_Responder) with 95% CI",
    x     = "TimePoint",
    y     = "Δmean (Responder − Non_Responder)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )



p_stats <- sig_df %>%
  group_by(Signature, T_NK_cell, TimePoint) %>%
  summarise(
    p = tryCatch(
      wilcox.test(value ~ Response)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    label = case_when(
      is.na(p_adj) ~ "",
      p_adj < 1e-4 ~ "****",
      p_adj < 1e-3 ~ "***",
      p_adj < 1e-2 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

annot_df <- summary_ci %>%
  left_join(p_stats, by = c("Signature", "T_NK_cell", "TimePoint")) %>%
  filter(label != "")
ggplot(summary_ci,
       aes(x = TimePoint, y = mean_diff_R_NR,
           group = T_NK_cell, color = T_NK_cell)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.15,
    linewidth = 0.4
  ) +
  geom_text(
    data = annot_df,
    aes(label = label),
    vjust = -1.0,
    size  = 3,
    show.legend = FALSE
  ) +
  facet_wrap(~ Signature, scales = "free_y") +
  scale_color_brewer(palette = "Set3") +
  labs(
    title = "Signature mean difference (R − NR) with 95% CI and adj.p",
    x     = "TimePoint",
    y     = "Δmean (Responder − Non_Responder)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

dodge <- position_dodge(width = 0.4)
set3_colors <- c(
  CD8_EFF   = "#FDB462",
  CD8_MAIT  = "#B3DE69",
  CD8_Temra = "#BC80BD",
  CD8_Naive = "#FCCDE5",
  CD8_pEx   = "#D9D9D9"
)

p1 <- ggplot(summary_ci,
       aes(x = TimePoint, y = mean_diff_R_NR,
           group = T_NK_cell, color = T_NK_cell)) +
  geom_line(position = dodge, size = 1) +
  geom_point(position = dodge, size = 2) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    position = dodge,
    width = 0.15,
    linewidth = 0.4
  ) +
  facet_wrap(~ Signature, scales = "free_y") +
  scale_color_manual(values = set3_colors) +
  labs(
    title = "Signature mean difference (R − NR) with 95% CI",
    x     = "TimePoint",
    y     = "Δmean (Responder − Non_Responder)"
  ) +
  theme_bw(base_size = 11) 
p1/p2/p3

library(dplyr)
library(tidyr)
library(ggplot2)

# 1) 시그니처들
sig_cols <- c(
  "Pro-apoptosis18",
  "Senescence8"
)
PBMC_CD8@meta.data %>% colnames()
# 2) cell type 지정
cell_types <- c("CD8_EFF", "CD8_MAIT", "CD8_Temra", "CD8_Naive", "CD8_pEx")

# 3) 셀 레벨 데이터 (Response 제거, T0~T3 모두 사용)
sig_df <- PBMC_CD8@meta.data %>%
  filter(T_NK_cell %in% cell_types) %>%
  mutate(TimePoint = sub("_.*", "", T_R)) %>%      # 기존 로직 유지
  filter(TimePoint %in% c("T0", "T1", "T2", "T3")) %>%
  select(T_NK_cell, TimePoint, all_of(sig_cols)) %>%
  pivot_longer(
    cols = all_of(sig_cols),
    names_to = "Signature",
    values_to = "value"
  ) %>%
  drop_na(value)

sig_df$TimePoint <- factor(sig_df$TimePoint, levels = c("T0", "T1", "T2", "T3"))

summary_tp <- sig_df %>%
  group_by(Signature, T_NK_cell, TimePoint) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = dplyr::n(),
    se = sd / sqrt(n),
    ci_low = mean - 1.96 * se,
    ci_high = mean + 1.96 * se,
    .groups = "drop"
  )

set3_colors <- c(
  CD8_EFF = "#FDB462",
  CD8_MAIT = "#B3DE69",
  CD8_Temra = "#BC80BD",
  CD8_Naive = "#FCCDE5",
  CD8_pEx = "#D9D9D9"
)

p_mean <- ggplot(summary_tp, aes(x = TimePoint, y = mean, group = T_NK_cell, color = T_NK_cell)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 0.4) +
  facet_wrap(~Signature, scales = "free_y") +
  scale_color_manual(values = set3_colors) +
  labs(
    title = "Signature mean by TimePoint (95% CI)",
    x = "TimePoint",
    y = "Mean module score"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

p_mean

baseline_df <- summary_tp %>%
  select(Signature, T_NK_cell, TimePoint, mean, sd, n) %>%
  pivot_wider(names_from = TimePoint, values_from = c(mean, sd, n)) %>%
  mutate(
    d_T1_T0 = mean_T1 - mean_T0,
    se_T1_T0 = sqrt((sd_T1^2 / n_T1) + (sd_T0^2 / n_T0)),
    lo_T1_T0 = d_T1_T0 - 1.96 * se_T1_T0,
    hi_T1_T0 = d_T1_T0 + 1.96 * se_T1_T0,
    
    d_T2_T0 = mean_T2 - mean_T0,
    se_T2_T0 = sqrt((sd_T2^2 / n_T2) + (sd_T0^2 / n_T0)),
    lo_T2_T0 = d_T2_T0 - 1.96 * se_T2_T0,
    hi_T2_T0 = d_T2_T0 + 1.96 * se_T2_T0,
    
    d_T3_T0 = mean_T3 - mean_T0,
    se_T3_T0 = sqrt((sd_T3^2 / n_T3) + (sd_T0^2 / n_T0)),
    lo_T3_T0 = d_T3_T0 - 1.96 * se_T3_T0,
    hi_T3_T0 = d_T3_T0 + 1.96 * se_T3_T0
  ) %>%
  select(Signature, T_NK_cell, starts_with("d_"), starts_with("lo_"), starts_with("hi_")) %>%
  pivot_longer(
    cols = starts_with("d_"),
    names_to = "Contrast",
    values_to = "diff"
  ) %>%
  mutate(
    TimePoint = recode(
      Contrast,
      d_T1_T0 = "T1-T0",
      d_T2_T0 = "T2-T0",
      d_T3_T0 = "T3-T0"
    ),
    ci_low = case_when(
      Contrast == "d_T1_T0" ~ lo_T1_T0,
      Contrast == "d_T2_T0" ~ lo_T2_T0,
      Contrast == "d_T3_T0" ~ lo_T3_T0
    ),
    ci_high = case_when(
      Contrast == "d_T1_T0" ~ hi_T1_T0,
      Contrast == "d_T2_T0" ~ hi_T2_T0,
      Contrast == "d_T3_T0" ~ hi_T3_T0
    )
  ) %>%
  select(Signature, T_NK_cell, TimePoint, diff, ci_low, ci_high)

baseline_df$TimePoint <- factor(baseline_df$TimePoint, levels = c("T1-T0", "T2-T0", "T3-T0"))

p_base <- ggplot(baseline_df, aes(x = TimePoint, y = diff, group = T_NK_cell, color = T_NK_cell)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 0.4) +
  facet_wrap(~Signature, scales = "free_y") +
  scale_color_manual(values = set3_colors) +
  labs(
    title = "Baseline contrasts (T1-T0, T2-T0, T3-T0) with 95% CI",
    x = "Contrast",
    y = "Δmean (vs T0)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

p_base

consec_df <- summary_tp %>%
  arrange(Signature, T_NK_cell, TimePoint) %>%
  group_by(Signature, T_NK_cell) %>%
  mutate(
    mean_prev = lag(mean),
    sd_prev = lag(sd),
    n_prev = lag(n),
    diff = mean - mean_prev,
    se = sqrt((sd^2 / n) + (sd_prev^2 / n_prev)),
    ci_low = diff - 1.96 * se,
    ci_high = diff + 1.96 * se,
    Contrast = paste0(TimePoint, "-", lag(TimePoint))
  ) %>%
  ungroup() %>%
  filter(Contrast %in% c("T1-T0", "T2-T1", "T3-T2"))

consec_df$Contrast <- factor(consec_df$Contrast, levels = c("T1-T0", "T2-T1", "T3-T2"))

p_consec <- ggplot(consec_df, aes(x = Contrast, y = diff, group = T_NK_cell, color = T_NK_cell)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 0.4) +
  facet_wrap(~Signature, scales = "free_y") +
  scale_color_manual(values = set3_colors) +
  labs(
    title = "Consecutive contrasts (T1-T0, T2-T1, T3-T2) with 95% CI",
    x = "Contrast",
    y = "Δmean (consecutive)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

p_consec






library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

meta <- PBMC_CD8@meta.data %>%
  filter(sample_name !=' RCC19_B_T3') %>%  #Responder
  mutate(Timepoint = factor(Timepoint, levels = c("T0","T1","T2","T3")))

subset_var  <- "T_NK_cell"
subset_keep <- c("CD8_EFF", "CD8_Temra")

modules_keep <- c("Pro-apoptosis18", "Senescence8") 

con_tbl <- tibble(
  contrast = c("T1-T0", "T2-T1", "T3-T2"),
  t0 = c("T0","T1","T2"),
  t1 = c("T1","T2","T3")
)

## =========================================================
## 2) cell-level long format
## =========================================================
df_long_cell <- meta %>%
  filter(!is.na(Timepoint)) %>%
  filter(.data[[subset_var]] %in% subset_keep) %>%
  select(Timepoint, subset = all_of(subset_var), all_of(modules_keep)) %>%
  pivot_longer(cols = all_of(modules_keep),
               names_to = "module", values_to = "score") %>%
  filter(is.finite(score)) %>%
  mutate(
    subset = factor(subset, levels = subset_keep),
    module = factor(module, levels = modules_keep)
  )

## =========================================================
## 3) Welch 95% CI for mean difference
##    - Δmean = mean(x1) - mean(x0)
## =========================================================
welch_ci <- function(x0, x1) {
  x0 <- x0[is.finite(x0)]
  x1 <- x1[is.finite(x1)]
  n0 <- length(x0); n1 <- length(x1)
  if (n0 < 2 || n1 < 2) return(c(NA_real_, NA_real_))
  
  m0 <- mean(x0); m1 <- mean(x1)
  v0 <- var(x0);  v1 <- var(x1)
  
  se <- sqrt(v0/n0 + v1/n1)
  
  # Welch-Satterthwaite df
  df <- (v0/n0 + v1/n1)^2 / ((v0^2)/((n0^2)*(n0-1)) + (v1^2)/((n1^2)*(n1-1)))
  tcrit <- qt(0.975, df = df)
  
  delta <- m1 - m0
  c(delta - tcrit*se, delta + tcrit*se)
}

stats_delta_cell <- df_long_cell %>%
  group_by(subset, module) %>%
  group_modify(~{
    dat <- .x
    
    bind_rows(lapply(seq_len(nrow(con_tbl)), function(i){
      t0 <- con_tbl$t0[i]
      t1 <- con_tbl$t1[i]
      ct <- con_tbl$contrast[i]
      
      v0 <- dat %>% filter(Timepoint == t0) %>% pull(score)
      v1 <- dat %>% filter(Timepoint == t1) %>% pull(score)
      
      n0 <- sum(is.finite(v0)); n1 <- sum(is.finite(v1))
      if (n0 < 5 || n1 < 5) {
        return(tibble(
          contrast = ct, n0 = n0, n1 = n1,
          delta_mean = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_
        ))
      }
      
      delta_mean <- mean(v1, na.rm = TRUE) - mean(v0, na.rm = TRUE)
      ci <- welch_ci(v0, v1)
      
      pval <- suppressWarnings(
        wilcox.test(v0, v1, paired = FALSE, exact = FALSE)$p.value
      )
      
      tibble(
        contrast = ct, n0 = n0, n1 = n1,
        delta_mean = delta_mean, ci_low = ci[1], ci_high = ci[2], p = pval
      )
    }))
  }) %>%
  ungroup() %>%
  mutate(contrast = factor(contrast, levels = con_tbl$contrast)) %>%
  group_by(subset, module) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%  
  ungroup()

panel_range <- stats_delta_cell %>%
  group_by(subset, module) %>%
  summarise(
    ymin = min(ci_low, ci_high, delta_mean, na.rm = TRUE),
    ymax = max(ci_low, ci_high, delta_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(pad = 0.10 * (ymax - ymin + 1e-6))

stats_delta_cell <- stats_delta_cell %>%
  left_join(panel_range, by = c("subset","module")) %>%
  mutate(
    y_label = ymax + pad,
    label = ifelse(
      is.na(p.adj),
      "NA",
      paste0("adj p =", formatC(p.adj, format = "e", digits = 2),
             "\n(n=", n0, ",", n1, ")")
    )
  )

## =========================================================
## Barplot + 95% CI + p-value
## =========================================================

time_cols <- c(
  "T1-T0" = "#66C2A5",   # Set2 1
  "T2-T1" = "#FC8D62",   # Set2 2
  "T3-T2" = "#8DA0CB"    # Set2 3
)

p_delta_bar_cell <- ggplot(stats_delta_cell, aes(x = contrast, y = delta_mean, fill = contrast)) +
  geom_hline(yintercept = 0, linetype = 1, linewidth = 0.4,colour = 'red') +
  geom_col(width = 0.65, color = "black", alpha = 0.9, na.rm = TRUE) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.5, na.rm = TRUE) +
  geom_text(aes(y = y_label, label = label),
            size = 3.2, vjust = 0, na.rm = TRUE) +
  facet_grid(subset ~ module, scales = "free_y") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    x = "Consecutive comparison",
    y = "Δ mean (cell-level; mean(T1) − mean(T0)) with 95% Welch CI"
  )

p_delta_bar_cell+
  scale_fill_manual(values = time_cols)















