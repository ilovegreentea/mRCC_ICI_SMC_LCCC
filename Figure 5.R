library(Seurat)
library(msigdbr)
H <- msigdbr(species = "Homo sapiens", category = "H")

hallmark_list <-
  H |>
  distinct(gs_name, gene_symbol) |>
  group_by(gs_name) |>
  summarise(genes = list(gene_symbol), .groups = "drop") |>
  tibble::deframe()

DimPlot(tissue_Myeloid,group.by = 'Myeloid')


tissue_Myeloid@meta.data$subtypes <- tissue_Myeloid@meta.data$Myeloid

C1 <- subset(tissue_Myeloid, subtypes == 'C1QC_TAM')
DefaultAssay(C1) <- 'RNA'
C1 <- NormalizeData(C1)
C1_T0_DE <- FindMarkers(C1,ident.1 = 'T0_Responder',ident.2 = 'T0_Non_Responder',group.by = 'T_R')

library(EnhancedVolcano)
C1_T0_DE <- filter(C1_T0_DE, pct.1 > 0.1 & pct.2 > 0.1)
C1_T0_DE %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 0.01 )
EnhancedVolcano(C1_T0_DE,lab = rownames(C1_T0_DE),x = 'avg_log2FC',y = 'p_val_adj',FCcutoff = 0.5,pCutoff = 0.01,drawConnectors = T,selectLab = c("ALOX5AP",'NFKBIZ','FCGR2B','NEAT1','HSPA1B','IER5','IL1B','HSPH1',
                                                                                                                                                  'DNAJB1','CD300A','CXCL8','CX3CR1','JUN','HSPA1A','VSIG4','CCL3'))

C3 <- subset(tissue_Myeloid, subtypes == 'GPNMB_TAM')
DefaultAssay(C3) <- 'RNA'
C3 <- NormalizeData(C3)
C3_T0_DE <- FindMarkers(C3,ident.1 = 'T0_Responder',ident.2 = 'T0_Non_Responder',group.by = 'T_R')

library(EnhancedVolcano)
C3_T0_DE <- filter(C3_T0_DE, pct.1 > 0.1 & pct.2 > 0.1)
C3_T0_DE %>% filter(abs(avg_log2FC) > 0.5 & p_val_adj < 0.01 )
EnhancedVolcano(C3_T0_DE,lab = rownames(C3_T0_DE),x = 'avg_log2FC',y = 'p_val_adj',FCcutoff = 0.5,pCutoff = 0.01,drawConnectors = T)
CXCL <- subset(newseurat,TAM_type  == 'C2_GPNMB_TAM')
CXCL <- AddModuleScore(CXCL,features = hallmark_list,name = names(hallmark_list))



library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
library(fgsea)
library(ggplot2)
library(stringr)
library(tidytext)
library(scales)

cell_types <- c("C1QC_TAM", "GPNMB_TAM",'CCL3_TAM')
timepoints <- c("T0", "T2")

# 1) celltype + timepoint 하나에 대해 fgsea Top/Bottom 10 뽑는 함수
one_fgsea <- function(seu, celltype, tp, pathways, n = 10) {
  obj <- subset(seu, subtypes == celltype)
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  
  de <- FindMarkers(
    obj,
    ident.1 = paste0(tp, "_Responder"),
    ident.2 = paste0(tp, "_Non_Responder"),
    group.by = "T_R"
  )
  
  # 너가 하던 pct 컷
  de <- de %>% filter(pct.1 > 0.1, pct.2 > 0.1)
  
  # ranks
  ranks <- de$avg_log2FC
  names(ranks) <- rownames(de)
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  
  fg <- fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500) %>%
    as_tibble() %>%
    arrange(padj) %>%
    mutate(
      cell_type = celltype,
      tp = tp,
      logFDR = -log10(padj + 1e-300),
      pathway_clean = str_wrap(str_replace_all(pathway, "_", " "), 34)
    )
  
  top <- fg %>% filter(NES > 0) %>% slice_max(NES, n = n, with_ties = FALSE) %>% mutate(extreme = "Top (NES>0)")
  bot <- fg %>% filter(NES < 0) %>% slice_min(NES, n = n, with_ties = FALSE) %>% mutate(extreme = "Bottom (NES<0)")
  
  bind_rows(top, bot) %>%
    mutate(pathway_plot = tidytext::reorder_within(pathway_clean, NES, interaction(cell_type, tp, extreme)))
}

# 2) 전부 돌려서 합치기
df_all <- crossing(cell_type = cell_types, tp = timepoints) %>%
  mutate(res = map2(cell_type, tp, ~one_fgsea(tissue_Myeloid, .x, .y, hallmark_list, n = 10))) %>%
  unnest(res)

# 3) plot
max_abs <- max(abs(df_all$NES), na.rm = TRUE)

ggplot(df_all, aes(NES, pathway_plot)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_col(aes(fill = NES), width = 0.85) +
  geom_point(aes(size = logFDR, color = tp), shape = 21, fill = "white", stroke = 0.7) +
  facet_grid(cell_type ~ extreme, scales = "free_y", space = "free_y") +
  tidytext::scale_y_reordered() +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-max_abs, max_abs), oob = scales::squish,
    name = "NES"
  ) +
  scale_color_manual(values = c(T0 = "#66C2A5", T2 = "#FC8D62"), name = "Timepoint") +
  scale_size_continuous(name = expression(-log[10]("FDR")), range = c(1.6, 5.5)) +
  labs(x = "NES", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )




library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(scales)

# df_all에서 top/bottom 줄이기 (예: top5/bottom5)
df_plot <- df_all %>%
  group_by(cell_type, tp, extreme) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    y_lab = pathway %>%
      str_replace_all("^HALLMARK_", "") %>%
      str_replace_all("_", " ") %>%
      str_wrap(26)
  )

# cell_type별로 pathway 정렬(모든 tp를 합쳐 |NES| 큰 순으로)
ord_tbl <- df_plot %>%
  group_by(cell_type, y_lab) %>%
  summarise(ord = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
  arrange(cell_type, desc(ord))

df_plot <- df_plot %>%
  left_join(ord_tbl, by = c("cell_type","y_lab")) %>%
  group_by(cell_type) %>%
  mutate(y_lab = fct_reorder(y_lab, ord, .desc = TRUE)) %>%
  ungroup()

max_abs <- max(abs(df_plot$NES), na.rm = TRUE)

ggplot(df_plot, aes(x = NES, y = y_lab)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
  geom_col(aes(fill = NES), width = 0.75) +
  geom_point(aes(size = logFDR), shape = 21, fill = "white", color = "black", stroke = 0.6) +
  facet_grid(cell_type ~ tp, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-max_abs, max_abs), oob = scales::squish, name = "NES"
  ) +
  scale_size_continuous(name = expression(-log[10]("FDR")), range = c(1.5, 5.2)) +
  labs(x = "NES (Responder vs Non-responder)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )
















library(Seurat)
library(dplyr)
library(fgsea)
library(tibble)
library(ggplot2)
library(stringr)
library(tidytext)
library(scales)

## --------------------------
## 0) 공통: fgsea -> top/bottom 10 뽑는 작은 블록(함수X)
##    (아래에서 copy/paste로 4번 쓸거라서 코드 중복 허용)
## --------------------------

## ========== 1) C1QC_TAM / T0 ==========
C1 <- subset(tissue_Myeloid, subtypes == "C1QC_TAM")
DefaultAssay(C1) <- "RNA"
C1 <- NormalizeData(C1, verbose = FALSE)

C1_T0_DE <- FindMarkers(
  C1,
  ident.1 = "T0_Responder",
  ident.2 = "T0_Non_Responder",
  group.by = "T_R"
)
C1_T0_DE <- C1_T0_DE %>% filter(pct.1 > 0.1, pct.2 > 0.1)

ranks <- C1_T0_DE$avg_log2FC
names(ranks) <- rownames(C1_T0_DE)
ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)

set.seed(1)
fg <- fgseaMultilevel(pathways = hallmark_list, stats = ranks, minSize = 15, maxSize = 500) %>%
  as_tibble() %>%
  arrange(padj) %>%
  mutate(
    cell_type = "C1QC_TAM",
    tp = "T0",
    logFDR = -log10(padj + 1e-300),
    pathway_clean = str_wrap(str_replace_all(pathway, "_", " "), 34)
  )

C1_T0_top <- fg %>% filter(NES > 0) %>% slice_max(NES, n = 10, with_ties = FALSE) %>% mutate(extreme = "Top (NES>0)")
C1_T0_bot <- fg %>% filter(NES < 0) %>% slice_min(NES, n = 10, with_ties = FALSE) %>% mutate(extreme = "Bottom (NES<0)")
C1_T0_df  <- bind_rows(C1_T0_top, C1_T0_bot)


## ========== 2) C1QC_TAM / T2 ==========
C1_T2_DE <- FindMarkers(
  C1,
  ident.1 = "T2_Responder",
  ident.2 = "T2_Non_Responder",
  group.by = "T_R"
)
C1_T2_DE <- C1_T2_DE %>% filter(pct.1 > 0.1, pct.2 > 0.1)

ranks <- C1_T2_DE$avg_log2FC
names(ranks) <- rownames(C1_T2_DE)
ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)

set.seed(1)
fg <- fgseaMultilevel(pathways = hallmark_list, stats = ranks, minSize = 15, maxSize = 500) %>%
  as_tibble() %>%
  arrange(padj) %>%
  mutate(
    cell_type = "C1QC_TAM",
    tp = "T2",
    logFDR = -log10(padj + 1e-300),
    pathway_clean = str_wrap(str_replace_all(pathway, "_", " "), 34)
  )

C1_T2_top <- fg %>% filter(NES > 0) %>% slice_max(NES, n = 10, with_ties = FALSE) %>% mutate(extreme = "Top (NES>0)")
C1_T2_bot <- fg %>% filter(NES < 0) %>% slice_min(NES, n = 10, with_ties = FALSE) %>% mutate(extreme = "Bottom (NES<0)")
C1_T2_df  <- bind_rows(C1_T2_top, C1_T2_bot)


## ========== 3) GPNMB_TAM / T0 ==========
C3 <- subset(tissue_Myeloid, subtypes == "GPNMB_TAM")
DefaultAssay(C3) <- "RNA"
C3 <- NormalizeData(C3, verbose = FALSE)

C3_T0_DE <- FindMarkers(
  C3,
  ident.1 = "T0_Responder",
  ident.2 = "T0_Non_Responder",
  group.by = "T_R"
)
C3_T0_DE <- C3_T0_DE %>% filter(pct.1 > 0.1, pct.2 > 0.1)

ranks <- C3_T0_DE$avg_log2FC
names(ranks) <- rownames(C3_T0_DE)
ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)

set.seed(1)
fg <- fgseaMultilevel(pathways = hallmark_list, stats = ranks, minSize = 15, maxSize = 500) %>%
  as_tibble() %>%
  arrange(padj) %>%
  mutate(
    cell_type = "GPNMB_TAM",
    tp = "T0",
    logFDR = -log10(padj + 1e-300),
    pathway_clean = str_wrap(str_replace_all(pathway, "_", " "), 34)
  )

C3_T0_top <- fg %>% filter(NES > 0) %>% slice_max(NES, n = 10, with_ties = FALSE) %>% mutate(extreme = "Top (NES>0)")
C3_T0_bot <- fg %>% filter(NES < 0) %>% slice_min(NES, n = 10, with_ties = FALSE) %>% mutate(extreme = "Bottom (NES<0)")
C3_T0_df  <- bind_rows(C3_T0_top, C3_T0_bot)


## ========== 4) GPNMB_TAM / T2 ==========
C3_T2_DE <- FindMarkers(
  C3,
  ident.1 = "T2_Responder",
  ident.2 = "T2_Non_Responder",
  group.by = "T_R"
)
C3_T2_DE <- C3_T2_DE %>% filter(pct.1 > 0.1, pct.2 > 0.1)

ranks <- C3_T2_DE$avg_log2FC
names(ranks) <- rownames(C3_T2_DE)
ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)

set.seed(1)
fg <- fgseaMultilevel(pathways = hallmark_list, stats = ranks, minSize = 15, maxSize = 500) %>%
  as_tibble() %>%
  arrange(padj) %>%
  mutate(
    cell_type = "GPNMB_TAM",
    tp = "T2",
    logFDR = -log10(padj + 1e-300),
    pathway_clean = str_wrap(str_replace_all(pathway, "_", " "), 34)
  )

C3_T2_top <- fg %>% filter(NES > 0) %>% slice_max(NES, n = 10, with_ties = FALSE) %>% mutate(extreme = "Top (NES>0)")
C3_T2_bot <- fg %>% filter(NES < 0) %>% slice_min(NES, n = 10, with_ties = FALSE) %>% mutate(extreme = "Bottom (NES<0)")
C3_T2_df  <- bind_rows(C3_T2_top, C3_T2_bot)

## ========== 5) CCL3_TAM / T0 ==========
C4 <- subset(tissue_Myeloid, subtypes == "CCL3_TAM")
DefaultAssay(C4) <- "RNA"
C4 <- NormalizeData(C4, verbose = FALSE)

C4_T0_DE <- FindMarkers(
  C4,
  ident.1 = "T0_Responder",
  ident.2 = "T0_Non_Responder",
  group.by = "T_R"
)
C4_T0_DE <- C4_T0_DE %>% filter(pct.1 > 0.1, pct.2 > 0.1)

ranks <- C4_T0_DE$avg_log2FC
names(ranks) <- rownames(C4_T0_DE)
ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)

set.seed(1)
fg <- fgseaMultilevel(
  pathways = hallmark_list,
  stats = ranks,
  minSize = 15,
  maxSize = 500
) %>%
  as_tibble() %>%
  arrange(padj) %>%
  mutate(
    cell_type = "CCL3_TAM",
    tp = "T0",
    logFDR = -log10(padj + 1e-300),
    pathway_clean = str_wrap(str_replace_all(pathway, "_", " "), 34)
  )

C4_T0_top <- fg %>% filter(NES > 0) %>% slice_max(NES, n = 10, with_ties = FALSE) %>% mutate(extreme = "Top (NES>0)")
C4_T0_bot <- fg %>% filter(NES < 0) %>% slice_min(NES, n = 10, with_ties = FALSE) %>% mutate(extreme = "Bottom (NES<0)")
C4_T0_df  <- bind_rows(C4_T0_top, C4_T0_bot)


## ========== 6) CCL3_TAM / T2 ==========
C4_T2_DE <- FindMarkers(
  C4,
  ident.1 = "T2_Responder",
  ident.2 = "T2_Non_Responder",
  group.by = "T_R"
)
C4_T2_DE <- C4_T2_DE %>% filter(pct.1 > 0.1, pct.2 > 0.1)

ranks <- C4_T2_DE$avg_log2FC
names(ranks) <- rownames(C4_T2_DE)
ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)

set.seed(1)
fg <- fgseaMultilevel(
  pathways = hallmark_list,
  stats = ranks,
  minSize = 15,
  maxSize = 500
) %>%
  as_tibble() %>%
  arrange(padj) %>%
  mutate(
    cell_type = "CCL3_TAM",
    tp = "T2",
    logFDR = -log10(padj + 1e-300),
    pathway_clean = str_wrap(str_replace_all(pathway, "_", " "), 34)
  )

C4_T2_top <- fg %>% filter(NES > 0) %>% slice_max(NES, n = 10, with_ties = FALSE) %>% mutate(extreme = "Top (NES>0)")
C4_T2_bot <- fg %>% filter(NES < 0) %>% slice_min(NES, n = 10, with_ties = FALSE) %>% mutate(extreme = "Bottom (NES<0)")
C4_T2_df  <- bind_rows(C4_T2_top, C4_T2_bot)

df_all <- bind_rows(
  C1_T0_df, C1_T2_df,
  C3_T0_df, C3_T2_df,
  C4_T0_df, C4_T2_df
) %>%
  mutate(
    pathway_plot = tidytext::reorder_within(pathway_clean, NES, interaction(cell_type, tp, extreme)),
    tp = factor(tp, levels = c("T0", "T2")),
    cell_type = factor(cell_type, levels = c("C1QC_TAM", "GPNMB_TAM", "CCL3_TAM"))
  )

max_abs <- max(abs(df_all$NES), na.rm = TRUE)

library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(scales)

# 1) Top5/Bottom5만 남기기 (이미 df_plot 만들었으면 이 블록은 스킵 가능)
df_plot <- df_all %>%
  group_by(cell_type, tp, extreme) %>%
  slice_head(n = 5) %>%
  ungroup()

# 2) 라벨 만들기 (대소문자 유지)
df_plot2 <- df_plot %>%
  mutate(
    y_lab = pathway %>%
      str_replace_all("^HALLMARK_", "") %>%
      str_replace_all("_", " ") %>%
      str_wrap(24)
  )

# 3) y축 순서: cell_type 안에서 |NES| 큰 애들이 위로 오게 (한 번만 정렬)
ord <- df_plot2 %>%
  group_by(cell_type, y_lab) %>%
  summarise(ord = max(abs(NES), na.rm = TRUE), .groups = "drop")

df_plot2 <- df_plot2 %>%
  left_join(ord, by = c("cell_type", "y_lab")) %>%
  group_by(cell_type) %>%
  mutate(y_lab = fct_reorder(y_lab, ord)) %>%
  ungroup()
library(dplyr)
library(forcats)

ord_tbl <- df_plot2 %>%
  group_by(cell_type, y_lab) %>%
  summarise(ord = max(abs(NES), na.rm = TRUE), .groups = "drop")


pos <- position_dodge2(width = 0.7, preserve = "single")
library(dplyr)
library(forcats)


lvl_tbl <- df_plot2 %>%
  group_by(cell_type, y_lab) %>%
  summarise(ord = max(NES, na.rm = TRUE), .groups = "drop") %>%
  arrange(cell_type, desc(ord))

library(stringr)

df_plot2_fix <- df_plot2 %>%
  left_join(lvl_tbl, by = c("cell_type", "y_lab")) %>%
  mutate(
    y_key = paste0(cell_type, "___", y_lab)  # 패널별로 유니크 y축
  )

# y_key 레벨 순서(=정렬) 만들기
y_levels <- lvl_tbl %>%
  mutate(y_key = paste0(cell_type, "___", y_lab)) %>%
  pull(y_key)

df_plot2_fix <- df_plot2_fix %>%
  mutate(
    y_key = factor(y_key, levels = y_levels)
  )

pos <- position_dodge2(width = 0.7, preserve = "single")
max_abs <- max(abs(df_plot2_fix$NES), na.rm = TRUE)

ggplot(df_plot2_fix, aes(x = NES, y = y_key, group = tp)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_col(aes(fill = NES), width = 0.65, position = pos) +
  geom_point(aes(size = logFDR, color = tp,fill = tp),
             shape = 21, stroke = 0.7, position = pos) +
  facet_grid(cell_type ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*?___", "", x)) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-max_abs, max_abs), oob = scales::squish,
    name = "NES"
  ) +
  scale_color_manual(values = c(T0 = "#66C2A5", T2 = "#FC8D62"), name = "Timepoint") +
  scale_size_continuous(name = expression(-log[10]("FDR")), range = c(1.6, 5.5)) +
  labs(x = "NES", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )
library(ggnewscale)

ggplot(df_plot2_fix, aes(x = NES, y = y_key, group = tp)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  
  # 1) 막대: fill = NES (연속형)
  geom_col(aes(fill = NES), width = 0.65, position = pos) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-max_abs, max_abs), oob = scales::squish,
    name = "NES"
  ) +
  
  ggnewscale::new_scale_fill() +   # ✅ fill 스케일 리셋
  
  # 2) 점: fill = tp (범주형)
  geom_point(
    aes(size = logFDR, fill = tp),
    shape = 21, color = 'gray', stroke = 0.7, position = pos
  ) +
  scale_fill_manual(values = c(T0 = "#66C2A5", T2 = "#FC8D62"), name = "Timepoint") +
  
  facet_grid(cell_type ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*?___", "", x)) +
  scale_size_continuous(name = expression(-log[10]("FDR")), range = c(1.6, 5.5)) +
  labs(x = "NES", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )













tissue_Myeloid <- JoinLayers(tissue_Myeloid)
saveRDS(tissue_Myeloid,'tissue_Myeloid.rds')






library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)


tissue_Myeloid@meta.data$subtypes %>% table()
tam_inhibitory_sig <- list(
  TAM_inhibitory = c("CD274", "PDCD1LG2", "CD47", "SIRPA", "VSIR", "HAVCR2", "LGALS9")
)

# object에 존재하는 gene만 남기기
tam_inhibitory_sig <- lapply(
  tam_inhibitory_sig,
  \(x) intersect(x, rownames(tissue_Myeloid))
)

# (체크) 너무 적게 남으면 경고
print(tam_inhibitory_sig)

tissue_Myeloid <- AddModuleScore(
  tissue_Myeloid,
  features = tam_inhibitory_sig,
  name = "TAM_inhibitory"
)

# 생성된 score 컬럼 확인 (보통 TAM_inhibitory1)
grep("^TAM_inhibitory", colnames(tissue_Myeloid@meta.data), value = TRUE)



features_fun <- c(
  "attract", "suppress",
  "M1_Curated", "M2_Curated",
  "IFNg.stimulated.Monocyte.derived.Macrophage",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "GO_RESPONSE_TO_TYPE_I_INTERFERON",
  "GO_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "HALLMARK_COMPLEMENT",
  "KEGG_PROTEASOME",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "GO_FC_RECEPTOR_SIGNALING_PATHWAY",
  "HALLMARK_TGF_BETA_SIGNALING",
  "Angiogenesis",
  "Phagocytosis",
  "TAM_inhibitory1"
)

# (선택) 없는 컬럼 있으면 제외
features_fun <- intersect(features_fun, colnames(tissue_Myeloid@meta.data))
features_fun


sub_keep <- c("CCL3_TAM", "C1QC_TAM", "GPNMB_TAM", "ABCA1_TAM", "Monocyte", "Prolif_TAM")

meta <- tissue_Myeloid@meta.data

mat <- meta %>%
  filter(subtypes %in% sub_keep) %>%
  select(subtypes, all_of(features_fun)) %>%
  group_by(subtypes) %>%
  summarise(across(all_of(features_fun), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("subtypes") %>%
  as.matrix() %>%
  t()   # ✅ rows=features, cols=subtypes


mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0




rownames(mat_z) <- rownames(mat_z) %>%
  str_replace_all("^HALLMARK_", "") %>%
  str_replace_all("^GO_", "") %>%
  str_replace_all("^KEGG_", "") %>%
  str_replace_all("_", " ") %>%
  str_wrap(28)

colnames(mat_z) <- colnames(mat_z) %>% str_wrap(12)



fun_group <- case_when(
  grepl("INTERFERON|IFN", rownames(mat_z), ignore.case = TRUE) ~ "IFN",
  grepl("IL6|STAT3|TGF", rownames(mat_z), ignore.case = TRUE) ~ "Cytokine",
  grepl("PROTEASOME|OXIDATIVE|FATTY|XENOBIOTIC", rownames(mat_z), ignore.case = TRUE) ~ "Metabolism",
  grepl("ANTIGEN|PRESENTATION", rownames(mat_z), ignore.case = TRUE) ~ "Antigen",
  grepl("PHAGOCYTOSIS|FC RECEPTOR", rownames(mat_z), ignore.case = TRUE) ~ "Phagocytosis",
  grepl("inhibitory", rownames(mat_z), ignore.case = TRUE) ~ "Checkpoint",
  TRUE ~ "Other"
)

# 1) 그룹 순서 지정 (원하는대로)
group_order <- c("Checkpoint", "Phagocytosis", "Antigen", "IFN", "Cytokine", "Metabolism", "Other")

# fun_group <- factor(fun_group, levels = group_order)
# 
# # 2) 행을 그룹 기준으로 정렬 (그룹 안에서는 현재 순서 유지)
# ord <- order(fun_group)
# mat_z2 <- mat_z[ord, , drop = FALSE]
# fun_group2 <- fun_group[ord]
# 
# # 3) row annotation도 같은 순서로
# row_ha2 <- rowAnnotation(
#   Group = fun_group2,
#   col = list(Group = group_cols),
#   annotation_name_gp = gpar(fontface = "bold"),
#   simple_anno_size = unit(3, "mm")
# )

# 4) heatmap (row clustering OFF)
ht <- Heatmap(
  mat_z2,
  name = "Row Z",
  col  = col_fun,
  # left_annotation = row_ha2,
  cluster_rows = FALSE,          # ✅ 핵심
  cluster_columns = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 9),
  rect_gp = gpar(col = "grey90", lwd = 0.5)
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")


library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)

# ---- (옵션) TAM inhibitory score를 DC에도 그대로 계산하고 싶으면 유지 ----
tam_inhibitory_sig <- list(
  TAM_inhibitory = c("CD274", "PDCD1LG2", "CD47", "SIRPA", "VSIR", "HAVCR2", "LGALS9")
)

tam_inhibitory_sig <- lapply(
  tam_inhibitory_sig,
  \(x) intersect(x, rownames(tissue_Myeloid))
)

print(tam_inhibitory_sig)

tissue_Myeloid <- AddModuleScore(
  tissue_Myeloid,
  features = tam_inhibitory_sig,
  name = "TAM_inhibitory"
)

grep("^TAM_inhibitory", colnames(tissue_Myeloid@meta.data), value = TRUE)

# ---- feature set ----
features_fun <- c(
  "attract", "suppress",
  "M1_Curated", "M2_Curated",
  "IFNg.stimulated.Monocyte.derived.Macrophage",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "GO_RESPONSE_TO_TYPE_I_INTERFERON",
  "GO_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "HALLMARK_COMPLEMENT",
  "KEGG_PROTEASOME",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "GO_FC_RECEPTOR_SIGNALING_PATHWAY",
  "HALLMARK_TGF_BETA_SIGNALING",
  "Angiogenesis",
  "Phagocytosis",
  "TAM_inhibitory1"
)

features_fun <- intersect(features_fun, colnames(tissue_Myeloid@meta.data))
features_fun

# ---- ✅ DC subtypes로 교체 ----
sub_keep <- c(
  "cDC1",
  "cDC2",
  "pDC",
  "mregDC_LAMP3(high)",
  "mregDC_LAMP3(low)"
)

meta <- tissue_Myeloid@meta.data

mat <- meta |>
  filter(subtypes %in% sub_keep) |>
  select(subtypes, all_of(features_fun)) |>
  group_by(subtypes) |>
  summarise(across(all_of(features_fun), \(x) mean(x, na.rm = TRUE)), .groups = "drop") |>
  column_to_rownames("subtypes") |>
  as.matrix() |>
  t()

mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0

rownames(mat_z) <- rownames(mat_z) |>
  str_replace_all("^HALLMARK_", "") |>
  str_replace_all("^GO_", "") |>
  str_replace_all("^KEGG_", "") |>
  str_replace_all("_", " ") |>
  str_wrap(28)

colnames(mat_z) <- colnames(mat_z) |> str_wrap(12)

# ---- row grouping (그대로 사용) ----
fun_group <- case_when(
  grepl("INTERFERON|IFN", rownames(mat_z), ignore.case = TRUE) ~ "IFN",
  grepl("IL6|STAT3|TGF", rownames(mat_z), ignore.case = TRUE) ~ "Cytokine",
  grepl("PROTEASOME|OXIDATIVE|FATTY|XENOBIOTIC", rownames(mat_z), ignore.case = TRUE) ~ "Metabolism",
  grepl("ANTIGEN|PRESENTATION", rownames(mat_z), ignore.case = TRUE) ~ "Antigen",
  grepl("PHAGOCYTOSIS|FC RECEPTOR", rownames(mat_z), ignore.case = TRUE) ~ "Phagocytosis",
  grepl("inhibitory", rownames(mat_z), ignore.case = TRUE) ~ "Checkpoint",
  TRUE ~ "Other"
)

group_order <- c("Checkpoint", "Phagocytosis", "Antigen", "IFN", "Cytokine", "Metabolism", "Other")
fun_group <- factor(fun_group, levels = group_order)

ord <- order(fun_group)
mat_z2 <- mat_z[ord, , drop = FALSE]
fun_group2 <- fun_group[ord]

# ---- (주의) 아래 두 객체는 너 코드에서 이미 정의돼 있어야 함 ----
# group_cols : named vector (levels(group_order) -> colors)
# col_fun    : colorRamp2(...) 같은 heatmap color function

row_ha2 <- rowAnnotation(
  Group = fun_group2,
  col = list(Group = group_cols),
  annotation_name_gp = gpar(fontface = "bold"),
  simple_anno_size = unit(3, "mm")
)
library(circlize)

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

ht <- Heatmap(
  mat_z2,
  name = "Row Z",
  col = col_fun,
  left_annotation = row_ha2,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 9),
  rect_gp = gpar(col = "grey90", lwd = 0.5)
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")





library(dplyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(patchwork)
tissue_Myeloid@meta.data$subtypes %>% table()

#-------------------------
# 0) proportion df 만들기 (앞에서 만든 subtype_prop 재사용 가능)
#-------------------------
meta <- tissue_Myeloid@meta.data %>%
  filter(Timepoint %in% c("T0","T2") & subtypes != 'pDC' & subtypes != 'mregDC_LAMP3(low)') %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("T0","T2")),
    Response  = factor(Response)
  )

sample_counts <- meta %>%
  group_by(sample_name) %>%
  summarise(total_cells = n(), .groups = "drop")

subtype_counts <- meta %>%
  group_by(sample_name, subtypes) %>%
  summarise(cell_count = n(), .groups = "drop")

df <- subtype_counts %>%
  left_join(sample_counts, by = "sample_name") %>%
  mutate(proportion = cell_count / total_cells) %>%
  left_join(
    meta %>% select(sample_name, Response, Timepoint, patients) %>% distinct(),
    by = "sample_name"
  ) %>%
  filter(!is.na(proportion), !is.na(Response), !is.na(Timepoint)) %>%
  droplevels()

time_cols <- c("T0" = "#66C2A5", "T2" = "#FC8D62")
x_ref <- "T0"   # geom_blank용

#-------------------------
library(dplyr)
library(rstatix)

# df는 sample × subtypes proportion 테이블 (이미 만든 df)
df <- df %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("T0","T2")),
    Response  = factor(Response)
  )

cap <- df %>%
  group_by(subtypes, Response) %>%
  summarise(y_max = max(proportion, na.rm = TRUE), .groups = "drop")

stat_all <- df %>%
  group_by(subtypes, Response) %>%
  filter(n_distinct(Timepoint) == 2) %>%
  wilcox_test(proportion ~ Timepoint, paired = FALSE) %>%
  ungroup() %>%
  group_by(Response,subtypes) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  left_join(cap, by = c("subtypes","Response")) %>%
  mutate(
    xmin = 1, xmax = 2,                         # ✅ 숫자로 고정 (T0=1, T2=2)
    y.position = y_max * 1.15,
    label = paste0("p=", signif(p_adj, 2))
  )


library(ggplot2)
library(ggpubr)
library(patchwork)

time_cols <- c("T0" = "#66C2A5", "T2" = "#FC8D62")

make_resp_plot <- function(resp_level, show_x = TRUE) {
  
  df_sub <- df %>% filter(Response == resp_level)
  
  # 유의한 것만 표시 (없을 수도 있음!)
  stat_sub <- stat_all %>%
    filter(Response == resp_level, !is.na(p_adj), p < 0.5)
  
  p <- ggplot(df_sub, aes(x = Timepoint, y = proportion, fill = Timepoint)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, color = "black") +
    geom_jitter(
      aes(color = Timepoint),
      position = position_jitter(width = 0.12, height = 0),
      size = 1.1, alpha = 0.8, show.legend = FALSE
    ) +
    facet_wrap(~ subtypes, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = time_cols) +
    scale_color_manual(values = time_cols) +
    theme_bw(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(face = "bold")
    ) +
    labs(
      title = resp_level,
      x = if (show_x) "Timepoint" else NULL,
      y = "Proportion"
    )
  
  # ✅ 유의한 결과가 있을 때만 브래킷 올리기 (없으면 안 올림 → 에러 방지)
  if (nrow(stat_sub) > 0) {
    p <- p +
      ggpubr::stat_pvalue_manual(
        data = stat_sub,
        label = "label",
        xmin = "xmin", xmax = "xmax",
        y.position = "y.position",
        tip.length = 0.01,
        bracket.size = 0.4,
        size = 3
      )
  }
  
  if (!show_x) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  p
}

p_non  <- make_resp_plot("Non_Responder", show_x = FALSE)
p_resp <- make_resp_plot("Responder",     show_x = TRUE)

p_non / p_resp

# 패널별 y_max (브래킷 위치용)
cap <- df %>%
  group_by(subtypes, Timepoint) %>%
  summarise(y_max = max(proportion, na.rm = TRUE), .groups = "drop")

# Timepoint별로 subtypes 안에서 NR vs R
stat_resp_tp <- df %>%
  group_by(Timepoint, subtypes) %>%
  filter(n_distinct(Response) == 2) %>%
  wilcox_test(proportion ~ Response, paired = FALSE) %>%
  ungroup() %>%
  group_by(Timepoint,subtypes) %>%                      # ✅ timepoint 내 BH
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  left_join(cap, by = c("subtypes", "Timepoint")) %>%
  mutate(
    xmin = "Non_Responder",
    xmax = "Responder",
    y.position = y_max * 1.15,
    label = paste0("p=", signif(p_adj, 2))
  )
time_cols_resp <- c("Non_Responder" = "#FC8D62", "Responder" = "#66C2A5")

df_T0 <- df %>% filter(Timepoint == "T0")
st_T0 <- stat_resp_tp %>% filter(Timepoint == "T0", !is.na(p_adj), p < 1)

# 브래킷이 잘리지 않게 headroom (유의한 결과 있을 때만)
blank_T0 <- st_T0 %>%
  dplyr::group_by(subtypes) %>%
  dplyr::summarise(y_blank = max(y.position, na.rm = TRUE) * 1.03, .groups = "drop") %>%
  mutate(Response = "Non_Responder")

p_T0 <- ggplot(df_T0, aes(x = Response, y = proportion, fill = Response)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, color = "black") +
  geom_jitter(width = 0.1, size = 1.1, alpha = 0.75) +
  facet_wrap(~ subtypes, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = time_cols_resp) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  ) +
  labs(title = "T0", x = NULL, y = "Proportion")

# 유의한 것 있을 때만 geom_blank + p-value
if (nrow(st_T0) > 0) {
  p_T0 <- p_T0 +
    geom_blank(data = blank_T0, aes(x = Response, y = y_blank), inherit.aes = FALSE) +
    ggpubr::stat_pvalue_manual(
      st_T0,
      label = "label",
      xmin = "xmin", xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.4,
      size = 3
    )
}

p_T0
df_T2 <- df %>% filter(Timepoint == "T2")
st_T2 <- stat_resp_tp %>% filter(Timepoint == "T2", !is.na(p_adj), p < 1)

blank_T2 <- st_T2 %>%
  group_by(subtypes) %>%
  summarise(y_blank = max(y.position, na.rm = TRUE) * 1.03, .groups = "drop") %>%
  mutate(Response = "Non_Responder")

p_T2 <- ggplot(df_T2, aes(x = Response, y = proportion, fill = Response)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, color = "black") +
  geom_jitter(width = 0.1, size = 1.1, alpha = 0.75) +
  facet_wrap(~ subtypes, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = time_cols_resp) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  ) +
  labs(title = "T2", x = "Response", y = "Proportion")

if (nrow(st_T2) > 0) {
  p_T2 <- p_T2 +
    geom_blank(data = blank_T2, aes(x = Response, y = y_blank), inherit.aes = FALSE) +
    ggpubr::stat_pvalue_manual(
      st_T2,
      label = "label",
      xmin = "xmin", xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.4,
      size = 3
    )
}

p_T2
p_T0/p_T2





