library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggpubr)
library(gt)

setwd("~/Library/CloudStorage/OneDrive-ArcusBiosciences,Inc/MZ/BFX intern")
fig.dir <- "~/Library/CloudStorage/OneDrive-ArcusBiosciences,Inc/MZ/BFX intern/results/integration/noGene_filter/slide_figure"
obj <- readRDS("~/Library/CloudStorage/OneDrive-ArcusBiosciences,Inc/MZ/BFX intern/processed_obj/GSE235711_CRS_annotated_seurat_obj.rds")
View(obj@meta.data)

# visualize major cell types
p <- DimPlot(
  obj,
  reduction = "umap.harmony",
  group.by = "celltype3",
  label = T,
  raster = TRUE
) +
  ggtitle("Major cell types")
p

ggsave(
  filename = file.path(fig.dir, "UMAP_major_celltypes.png"),
  plot = p,
  width = 8,
  height = 6
)

markers <- c(
  "CD3E","CD3D","TRAC","KLRD1","GNLY","KLRB1",
  "MZB1","JCHAIN","IGHG1","IGKC",
  "CLEC4C","IL3RA","LILRB4",
  "CXCR2","FCGR3B","IL1R2",
  "ITGAX","CPA3","TPSB2","TPSAB1",
  "CMA1","MS4A2","KIT",
  "COL1A1","COL1A2",
  "KRT7","EPCAM","KRT18",
  "VWF","PLVAP","PECAM1",
  "MS4A1","CD79A","CD19"
)

cell_order <- c(
  "T/NK cells",
  "Plasma cells",
  "pDCs",
  "Neutrophil",
  "Myeloid",
  "Mast cells",
  "Fibroblasts",
  "Epithelial cells",
  "Endothelial cells",
  "B cells"
)

Idents(obj) <- "celltype"

## dotplots: markers by celltype
p <- DotPlot(
  obj,
  features = markers
  #cols = c("thistle2", "blue")
) +
  scale_y_discrete(limits = rev(cell_order)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    #panel.background = element_rect(fill = "grey92"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.title = element_blank()
  )
p

##----cell type composition across condition----
df_summary <- obj@meta.data %>%
  group_by(condition2, orig.ident, celltype) %>%
  summarise(n = n(), .groups = "drop")

df_summary <- df_summary %>%
  group_by(condition2, orig.ident) %>%
  mutate(pct = n / sum(n) * 100)

p <- ggplot(df_summary, aes(x = orig.ident, y = pct, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.9) +
  facet_wrap(~condition2, nrow = 1, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "Patient", y = "Percent", fill = NULL) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey80", color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p

ggsave(
  filename = file.path(fig.dir, "celltype_composition_bycondition.png"),
  plot = p,
  width = 8,
  height = 4
)

##----mast cell proportion across condition----
df_summary <- obj@meta.data %>%
  group_by(orig.ident, condition2) %>%
  summarise(
    total_cells = n(),
    mast_cells = sum(celltype == "Mast cells"),
    pct_mast = mast_cells / total_cells * 100,
    .groups = "drop"
  )

p <- ggplot(df_summary, aes(x = condition2, y = pct_mast, fill = condition2)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, color = "black") +
  scale_fill_manual(values = c("Control"="#ffadad",
                               "CRSsNP_Eth"="#b6e2d3",
                               "CRSwNP_Eth"="#809bce", 
                               "CRSwNP_NP" = "#fbaf87"
  ))+
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     comparisons = list(
                       c("Control", "CRSsNP_Eth"),
                       c("CRSsNP_Eth", "CRSwNP_Eth"),
                       c("CRSwNP_Eth", "CRSwNP_NP")
                     ))+
  labs(
    title = "Mast Cell Proportion Across Conditions",
    x = "Condition",
    y = "Mast Cell Percentage"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
p

ggsave(
  filename = file.path(fig.dir, "mast_composition_bycondition.png"),
  plot = p,
  width = 6.5,
  height = 4
)

##----mast cell subtype proportion across condition----
df_mast <- obj@meta.data %>%
  filter(celltype == "Mast cells") 

# Count per subtype per condition
df_mast_summary <- df_mast %>%
  group_by(celltype3, condition2) %>%
  summarise(n = n(), .groups = "drop")

# Calculate percentage within each condition
df_mast_summary <- df_mast_summary %>%
  group_by(condition2) %>%
  mutate(pct = n / sum(n) * 100)

df_mast_summary <- df_mast_summary %>%
  mutate(label = paste0(n, " (", sprintf("%.1f", pct), "%)"))

###----summary table----
table_wide <- df_mast_summary %>%
  dplyr::select(celltype3, condition2, label) %>%
  tidyr::pivot_wider(names_from = condition2, values_from = label)

total_counts <- df_mast %>%
  group_by(celltype3) %>%
  summarise(total_n = n(), .groups = "drop")

table_final <- table_wide %>%
  left_join(total_counts, by = "celltype3") %>%
  dplyr::rename(`Mast subtype` = celltype3,
         `Total mast cells` = total_n)

gt_tbl <- table_final %>%
  gt() %>%
  #tab_header(
  #title = "Mast Cell Subtype Distribution Across Conditions")
  cols_align(
    align = "center",
    -`Mast subtype`
  ) %>%
  cols_label(
    `Total mast cells` = "Total"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  )
gt_tbl

mat <- matrix(
  c(
    205, 135, 639, 15,  # MCT
    29,  13,  90,  1    # MCTC
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(mat) <- c("MCT", "MCTC")
colnames(mat) <- c("CRSsNP_Eth", "CRSwNP_Eth", "CRSwNP_NP", "Control")

chisq.test(mat)

###----pie chart----
p <- ggplot(df_mast_summary, aes(x = "", y = pct, fill = celltype3)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~condition2, nrow = 1) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_stack(vjust = 0.5),
            size = 3.5) +
  theme_void() +
  labs(title = "",
       fill = "Subtype")
p
ggsave(
  filename = file.path(fig.dir, "Mast_subtype_pct_piechart.png"),
  plot = p,
  width = 8,
  height = 4
)

##----X2 signature in mast cells----
X2_signature <- c('ARL5B', 'CREM', 'LDLR', 'MAPK6', 'NR4A1', 'RRM2', 
                  'SGK1', 'SOCS3', 'TBX20', 'PNOC', 'IL13', 'SEC14L2', 'ARC')

p <- VlnPlot(obj, features = X2_signature, group.by = 'celltype', ncol = 5)
p

ggsave(
  filename = file.path(fig.dir, "X2_signature_allcell.png"),
  plot = p,
  width = 20,
  height = 12
)

##----IL13 expression----
df <- obj@meta.data
df$IL13 <- FetchData(obj, vars = "IL13")[,1]

df_summary <- df %>%
  group_by(celltype, condition2) %>%
  summarise(IL13_pct_expr = mean(IL13 > 0) * 100)

p <- ggplot(df_summary, aes(x = condition2, y = IL13_pct_expr, fill = condition2)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_wrap(~ celltype, scales = "free_x", nrow = 2) +
  scale_fill_manual(values = c("Control"="#ffadad",
                               "CRSsNP_Eth"="#b6e2d3",
                               "CRSwNP_Eth"="#809bce", 
                               "CRSwNP_NP" = "#fbaf87"
  ))+
  theme_classic() +
  labs(x = "Condition", y = "Percent expressing IL13",
       fill = "Condition")+
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),  
    legend.text  = element_text(size = 12)
  )
p

df_summary <- df %>% filter(celltype=="Mast cells") %>%
  group_by(sample_id) %>%
  summarise(IL13_pct_expr = mean(IL13 > 0) * 100)
p <- ggplot(df_summary, aes(x = sample_id, y = IL13_pct_expr, fill = sample_id)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_classic() +
  labs(x = "Condition", y = "Percent expressing IL13",
       fill = "Condition")+
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),  
    legend.text  = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p

ggsave(
  filename = file.path(fig.dir, "IL13_pct_mast_bycondition.png"),
  plot = p,
  width = 15,
  height = 5
)

mast_obj <- subset(obj, subset = celltype=="Mast cells")
df <- FetchData(
  mast_obj,
  vars = c("IL13", "condition2", "sample_id")
)

df <- df %>% mutate(IL13_positive = IL13>0)

mast_il13_sample <- df %>%
  group_by(sample_id, condition2) %>%
  summarise(n_mast_cells = n(),
            n_IL13_pos = sum(IL13_positive, na.rm=T),
            percent_IL13_pos = 100*n_IL13_pos / n_mast_cells,
            mean_IL13_expr = mean(IL13, na.rm=T),
            .group='drop')
mast_il13_sample

pairwise.wilcox.test(
  x = mast_il13_sample$percent_IL13_pos,
  g = mast_il13_sample$condition2,
  method = "BH"
)

##----mean X2 score---- 
###----mast cells vs others----
df <- obj@meta.data[,c("sample_id","celltype", "celltype3","condition", "condition2", "X2_signature")]
df$is_mast <- ifelse(df$celltype=='Mast cells', 'Mast cells', 'Others')

df_patient <- df %>% group_by(sample_id, is_mast) %>%
  summarise(n_mast_cells = n(),
            mean_X2_signature = mean(X2_signature, na.rm=T),
            median_X2_signature = median(X2_signature, na.rm=T),
            .groups = "drop")
df_patient

p <- ggviolin(
  df_patient,
  x = "is_mast",
  y = "mean_X2_signature",
  fill = "is_mast",
  add = "boxplot",
  add.params = list(width = 0.15)   # control jitter spread
) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +  # more control
  stat_compare_means(label.x.npc = "right")
p

p <- ggplot(df, aes(x=condition2, y=X2_signature, fill=condition2))+
  geom_violin()+
  facet_wrap(~is_mast)+
  stat_compare_means(comparisons = list(
    c("Control", "CRSsNP_Eth"),
    c("Control", "CRSwNP_Eth"),
    c("Control", "CRSwNP_NP"),
    c("CRSsNP_Eth", "CRSwNP_Eth"),
    c("CRSsNP_Eth", "CRSwNP_NP"),
    c("CRSwNP_Eth", "CRSwNP_NP")
  ))+
  theme_classic()
p

df_crswnp_patient <- df %>% filter(condition=="CRSwNP") %>%
  group_by(sample_id, is_mast) %>%
  summarise(n_mast_cells = n(),
            mean_X2_signature = mean(X2_signature, na.rm=T),
            median_X2_signature = median(X2_signature, na.rm=T),
            .groups = "drop")
df_crswnp_patient

p <- ggviolin(
  df_crswnp_patient,
  x = "is_mast",
  y = "mean_X2_signature",
  fill = "is_mast",
  #color = "condition2",
  add = "boxplot",
  #position = position_dodge(0.8)
  add.params = list(width = 0.15)   # control jitter spread
) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +  # more control
  stat_compare_means(label.x.npc = "right")
p

df_crswnp_eth_patient <- df %>% filter(condition2=="CRSwNP_Eth") %>%
  group_by(sample_id, is_mast) %>%
  summarise(n_mast_cells = n(),
            mean_X2_signature = mean(X2_signature, na.rm=T),
            median_X2_signature = median(X2_signature, na.rm=T),
            .groups = "drop")
df_crswnp_eth_patient

p <- ggviolin(
  df_crswnp_eth_patient,
  x = "is_mast",
  y = "mean_X2_signature",
  title = "CRSwNP_Eth",
  fill = "is_mast",
  #color = "condition2",
  add = "boxplot",
  #position = position_dodge(0.8)
  add.params = list(width = 0.15)   # control jitter spread
) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +  # more control
  stat_compare_means(label.x.npc = "right")
p

df_crswnp_np_patient <- df %>% filter(condition2=="CRSwNP_NP") %>%
  group_by(sample_id, is_mast) %>%
  summarise(n_mast_cells = n(),
            mean_X2_signature = mean(X2_signature, na.rm=T),
            median_X2_signature = median(X2_signature, na.rm=T),
            .groups = "drop")
df_crswnp_np_patient

p <- ggviolin(
  df_crswnp_np_patient,
  x = "is_mast",
  y = "mean_X2_signature",
  title = "CRSwNP_NP",
  fill = "is_mast",
  #color = "condition2",
  add = "boxplot",
  #position = position_dodge(0.8)
  add.params = list(width = 0.15)   # control jitter spread
) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +  # more control
  stat_compare_means(label.x.npc = "right")
p

df_patient <- df %>%
  filter(condition == "CRSwNP") %>%
  group_by(sample_id, condition2, is_mast) %>%
  summarise(
    n_cells = n(),
    mean_X2_signature = mean(X2_signature, na.rm = TRUE),
    median_X2_signature = median(X2_signature, na.rm = TRUE),
    .groups = "drop"
  )

p <- ggplot(df_patient, aes(x = is_mast, y = mean_X2_signature, fill = is_mast)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  facet_wrap(~ condition2) +   # 
  theme_classic()+
  stat_compare_means(
  method = "wilcox.test",
  label = "p.format",
  label.x.npc = "right"
)
p

ggsave(
  filename = file.path(fig.dir, "mean_X2_signature_CRSwNP.png"),
  plot = p,
  width = 5.5,
  height = 4
)


###----mast cells: compared by condition----
mast_df <- df[df$celltype=="Mast cells",]

df_patient_mast <- mast_df %>% group_by(sample_id, condition2, celltype3) %>%
  summarise(n_mast_cells = n(),
            mean_X2_signature = mean(X2_signature, na.rm=T),
            median_X2_signature = median(X2_signature, na.rm=T),
            .groups = "drop")
df_patient_mast

p <- ggplot(mast_df, aes(x=condition2, y=X2_signature, fill=condition2))+
  geom_violin()+
  geom_boxplot(width=0.1, outlier.shape = NA)+
  scale_fill_manual(values = c("Control"="#ffadad",
                               "CRSsNP_Eth"="#b6e2d3",
                               "CRSwNP_Eth"="#809bce", 
                               "CRSwNP_NP" = "#fbaf87"
  ))+
  #stat_compare_means(method = "wilcox.test", comparisons = list(
  # c("Control", "CRSsNP"),
  # c("CRSwNP", "CRSsNP"),
  #c("Control", "CRSwNP")
  #))+
  theme_classic(base_size = 12)+
  labs(title = "MRGPRX2 signature in Mast cells by condition",
       x = "Condition",
       fill="Condition")+
  theme(plot.title = element_text(hjust = 0.5))
p

ggsave(
  filename = file.path(fig.dir, "X2_signature_mast_bycondition.png"),
  plot = p,
  width = 8,
  height = 5
)

p <- ggplot(df_patient_mast, aes(x=condition2, y=mean_X2_signature, fill=condition2))+
  geom_boxplot(width=0.5, outlier.shape = NA, alpha = 0.6)+
  geom_jitter(width = 0.15, size = 3, alpha = 0.9)+
  scale_fill_manual(values = c("Control"="#ffadad",
                               "CRSsNP_Eth"="#b6e2d3",
                               "CRSwNP_Eth"="#809bce", 
                               "CRSwNP_NP" = "#fbaf87"
  ))+
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     comparisons = list(
                       c("Control", "CRSsNP_Eth"),
                       c("CRSsNP_Eth", "CRSwNP_Eth"),
                       c("CRSwNP_Eth", "CRSwNP_NP"),
                       c("Control", "CRSwNP_Eth"),
                       c("CRSsNP_Eth", "CRSwNP_NP"),
                       c("Control", "CRSwNP_NP")
                     ))+
  # stat_compare_means(method = "kruskal.test",
  # label.y = max(df.patient.mast$median_X2_signature, na.rm=TRUE)*1.25)+
  theme_classic(base_size = 12)+
  labs(title = "Mean MRGPRX2 signature in Mast cells by condition",
       x = "Condition",
       fill="Condition")+
  theme(plot.title = element_text(hjust = 0.5))
p

ggsave(
  filename = file.path(fig.dir, "X2_signature_mast_bycondition_sample_level.png"),
  plot = p,
  width = 7,
  height = 5
)


ggplot(
  mast_df,
  aes(
    x = celltype3,
    y = X2_signature,
    fill = celltype3
  )
) +
  geom_violin(
    trim = FALSE,
    color = "black"
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white"
  ) +
  facet_wrap(
    ~ condition2,
    ncol = 2
  ) +
  scale_fill_manual(
    values = c(
      "MCT" = "#F8766D",
      "MCTC"     = "#619CFF"
    )
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_bw() +
  labs(
    x = "Mast cell subtype",
    y = "X2 signature",
    fill = "Celltype"
  )

####----explore which genes drive X2 score increase----
mast_obj <- subset(obj, subset = celltype == "Mast cells")
genes <- c('ARL5B', 'CREM', 'LDLR', 'MAPK6', 'NR4A1', 'RRM2', 'SGK1', 'SOCS3', 'TBX20', 'PNOC', 'IL13', 'SEC14L2', 'ARC')
print(genes)
length(genes)

Idents(mast_obj) <- 'condition'
de <- FindMarkers(mast_obj, ident.1 = "CRSwNP", ident.2 = "Control",
                  features = genes, logfc.threshold = 0, min.pct = 0)
de <- de[order(-de$avg_log2FC),]
de

mast_crsnp <- subset(mast_obj, subset = condition == "CRSwNP")
mast_ctr <- subset(mast_obj, subset = condition == "Control")
expr <- GetAssayData(mast_obj, layer="data")[genes, ]
expr_crsnp <- GetAssayData(mast_crsnp, layer="data")[genes, ]
expr_ctr <- GetAssayData(mast_ctr, layer="data")[genes, ]
X2_score <- mast_obj$X2_signature
X2_score_crsnp <- mast_crsnp$X2_signature
X2_score_ctr <- mast_ctr$X2_signature

cor_vals <- apply(expr, 1, function(x) cor(x, X2_score))
sort(cor_vals, decreasing=T)
cor_vals_crsnp <- apply(expr_crsnp, 1, function(x) cor(x, X2_score_crsnp))
sort(cor_vals_crsnp, decreasing=T)
cor_vals_ctr <- apply(expr_ctr, 1, function(x) cor(x, X2_score_ctr))
sort(cor_vals_ctr, decreasing=T)

DotPlot(mast_obj, features = genes, group.by = 'condition2') + RotatedAxis()

all_genes <- rownames(mast_obj)
mast_obj <- ScaleData(mast_obj, features = all_genes)
DoHeatmap(mast_obj, features = genes, group.by = 'condition2')

###----myeloid cells: compared by condition----
myeloid_df <- df[df$celltype=="Myeloid",]
df_patient_myeloid <- myeloid_df %>% group_by(sample_id, condition2) %>%
  summarise(n_myeloid_cells = n(),
            mean_X2_signature = mean(X2_signature, na.rm=T),
            median_X2_signature = median(X2_signature, na.rm=T),
            .groups = "drop")
df_patient_myeloid

p <- ggplot(myeloid_df, aes(x=condition2, y=X2_signature, fill=condition2))+
  geom_violin()+
  geom_boxplot(width=0.1, outlier.shape = NA)+
  scale_fill_manual(values = c("Control"="#ffadad",
                               "CRSsNP_Eth"="#b6e2d3",
                               "CRSwNP_Eth"="#809bce", 
                               "CRSwNP_NP" = "#fbaf87"
  ))+
  #stat_compare_means(method = "wilcox.test", comparisons = list(
  # c("Control", "CRSsNP"),
  # c("CRSwNP", "CRSsNP"),
  #c("Control", "CRSwNP")
  #))+
  theme_classic(base_size = 12)+
  labs(title = "MRGPRX2 signature in Myeloid cells by condition",
       x = "Condition",
       fill="Condition")+
  theme(plot.title = element_text(hjust = 0.5))
p

ggsave(
  filename = file.path(fig.dir, "X2_signature_Myeloid_bycondition.png"),
  plot = p,
  width = 8,
  height = 5
)

p <- ggplot(df_patient_myeloid, aes(x=condition2, y=mean_X2_signature, fill=condition2))+
  geom_boxplot(width=0.5, outlier.shape = NA, alpha = 0.6)+
  geom_jitter(width = 0.15, size = 3, alpha = 0.9)+
  scale_fill_manual(values = c("Control"="#ffadad",
                               "CRSsNP_Eth"="#b6e2d3",
                               "CRSwNP_Eth"="#809bce", 
                               "CRSwNP_NP" = "#fbaf87"
  ))+
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     comparisons = list(
                       c("Control", "CRSsNP_Eth"),
                       c("CRSsNP_Eth", "CRSwNP_Eth"),
                       c("CRSwNP_Eth", "CRSwNP_NP"),
                       c("Control", "CRSwNP_Eth"),
                       c("CRSsNP_Eth", "CRSwNP_NP"),
                       c("Control", "CRSwNP_NP")
                     ))+
  # stat_compare_means(method = "kruskal.test",
  # label.y = max(df.patient.mast$median_X2_signature, na.rm=TRUE)*1.25)+
  theme_classic(base_size = 12)+
  labs(title = "Mean Myeloid signature in Mast cells by condition",
       x = "Condition",
       fill="Condition")+
  theme(plot.title = element_text(hjust = 0.5))
p

ggsave(
  filename = file.path(fig.dir, "X2_signature_myeloid_bycondition_sample_level.png"),
  plot = p,
  width = 7,
  height = 5
)

###----epithelial cells: compared by condition----
epi_df <- df[df$celltype=="Epithelial cells",]
df_patient_epi <- epi_df %>% group_by(sample_id, condition2) %>%
  summarise(n_epi_cells = n(),
            mean_X2_signature = mean(X2_signature, na.rm=T),
            median_X2_signature = median(X2_signature, na.rm=T),
            .groups = "drop")
df_patient_epi

p <- ggplot(epi_df, aes(x=condition2, y=X2_signature, fill=condition2))+
  geom_violin()+
  geom_boxplot(width=0.1, outlier.shape = NA)+
  scale_fill_manual(values = c("Control"="#ffadad",
                               "CRSsNP_Eth"="#b6e2d3",
                               "CRSwNP_Eth"="#809bce", 
                               "CRSwNP_NP" = "#fbaf87"
  ))+
  #stat_compare_means(method = "wilcox.test", comparisons = list(
  # c("Control", "CRSsNP"),
  # c("CRSwNP", "CRSsNP"),
  #c("Control", "CRSwNP")
  #))+
  theme_classic(base_size = 12)+
  labs(title = "MRGPRX2 signature in Epithelial cells by condition",
       x = "Condition",
       fill="Condition")+
  theme(plot.title = element_text(hjust = 0.5))
p

ggsave(
  filename = file.path(fig.dir, "X2_signature_epi_bycondition.png"),
  plot = p,
  width = 8,
  height = 5
)

p <- ggplot(df_patient_epi, aes(x=condition2, y=mean_X2_signature, fill=condition2))+
  geom_boxplot(width=0.5, outlier.shape = NA, alpha = 0.6)+
  geom_jitter(width = 0.15, size = 3, alpha = 0.9)+
  scale_fill_manual(values = c("Control"="#ffadad",
                               "CRSsNP_Eth"="#b6e2d3",
                               "CRSwNP_Eth"="#809bce", 
                               "CRSwNP_NP" = "#fbaf87"
  ))+
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     comparisons = list(
                       c("Control", "CRSsNP_Eth"),
                       c("CRSsNP_Eth", "CRSwNP_Eth"),
                       c("CRSwNP_Eth", "CRSwNP_NP"),
                       c("Control", "CRSwNP_Eth"),
                       c("CRSsNP_Eth", "CRSwNP_NP"),
                       c("Control", "CRSwNP_NP")
                     ))+
  # stat_compare_means(method = "kruskal.test",
  # label.y = max(df.patient.mast$median_X2_signature, na.rm=TRUE)*1.25)+
  theme_classic(base_size = 12)+
  labs(title = "Mean Epithelial cell signature in Mast cells by condition",
       x = "Condition",
       fill="Condition")+
  theme(plot.title = element_text(hjust = 0.5))
p

ggsave(
  filename = file.path(fig.dir, "X2_signature_epi_bycondition_sample_level.png"),
  plot = p,
  width = 7,
  height = 5
)


## test low_mast_related genes expression level
low_mast_related <- c('IL1RL1', 'IL18R1', 'HPGDS', 'ALOX5', 'LTC4S',
                      'SIGLEC6', 'CD200R1', 'GATA2', 'PTGDR2', 'ADORA3')
low_mast_related <- c('CD3D', 'MMS4A1', 'CD14', 'MPO', 'FOXP3', 'GATA3')
mast_pos <- c('TPSAB1', 'TPSB2', 'CPA3', 'MS4A2')
gene_panel <- c(low_mast_related, "MRGPRX2", mast_pos)

mast <- subset(obj, subset = celltype == "Mast cells")
counts <- GetAssayData(mast, assay = "RNA", layer = "counts")
dim(counts)

gene_panel <- intersect(gene_panel, rownames(counts))
gene_summary <- data.frame(
  gene=gene_panel,
  pct_detected=Matrix::rowMeans(counts[gene_panel,,drop=F]>0)*100,
  mean_counts=Matrix::rowMeans(counts[gene_panel,,drop=F])
)
gene_summary <- gene_summary %>% arrange(pct_detected)
gene_summary

gene_summary <- gene_summary %>%
  mutate(
    category = case_when(
      gene %in% c('CD3D', 'MMS4A1', 'CD14', 'MPO', 'FOXP3', 'GATA3') ~ "Background",
      gene == "MRGPRX2" ~ "Target gene",
      TRUE ~ "Positive"
    )
  )

gene_summary <- gene_summary %>%
  mutate(
    gene = factor(gene, levels = gene[order(pct_detected)])
  )


p <- ggplot(gene_summary, aes(x=pct_detected, y=gene, fill=category))+
  geom_col(width=0.7)+
  geom_text(aes(label=sprintf("%.1f%%", pct_detected)),
             hjust=-0.05, size=4)+
  scale_fill_manual(
    values=c("Background" = "#F4A582",
    "Target gene" = "#4393C3",
    "Positive" = "#D6604D"
  ))+
  scale_x_continuous()+
  labs(title = "Gene detection rate in mast cells",
       x="% cells with detected expression",
       y=NULL,
       fill=NULL)+
  theme_minimal(base_size = 14)+
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )
p

ggsave(
  filename = file.path(fig.dir, "gene_dection_rate.png"),
  plot = p,
  width = 10,
  height = 5
)

##----Fibroblastsubclustering----
fib <- subset(obj, celltype == "Fibroblasts")

umap_fib <- Embeddings(fib, "umap.harmony")

fib$island <- ifelse(
  umap_fib[,1] < 13, "Fib1", "Fib2"
)

DimPlot(
  fib,
  reduction = "umap.harmony",
  group.by = "island",
  label = T,
  #raster = TRUE
) 

table(fib$island, fib$sample_id)
table(fib$island, fib$condition2)

df_summary <- fib@meta.data %>%
  group_by(condition2, sample_id, island) %>%
  summarise(n = n(), .groups = "drop")

df_summary <- df_summary %>%
  group_by(condition2, sample_id) %>%
  mutate(pct = n / sum(n) * 100)

p <- ggplot(df_summary, aes(x = sample_id, y = pct, fill = island)) +
  geom_bar(stat = "identity", width = 0.9) +
  facet_wrap(~condition2, nrow = 1, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "Patient", y = "Percent", fill = NULL) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey80", color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p

FeaturePlot(fib, features = c("COL1A1", "COL1A2", "POSTN", "CXCL14"))
FeaturePlot(fib, features = c("SAA1", "PRG4", "MFAP5"))
FeaturePlot(fib, features = c("RPLP0", "RPS3", "HSPA1A", "JUN", "FOS", "ATF3", "DDIT3"))
FeaturePlot(fib, features = c("RGS5", "MCAM", "CSPG4", "PDGFRB","THY1", "DES","ACTA2"))
FeaturePlot(fib, features = c("COL3A1", "COL6A1", "COL6A2", "DCN", "LUM", "FBLN1"))

VlnPlot(fib,
        features = c(
          "nCount_RNA",
          "nFeature_RNA",
          "percent.mt",
          "percent.ribo"
        ),
        group.by = 'island',
        pt.size=0, ncol = 2)


VlnPlot(obj,
        features = c(
          "nCount_RNA",
          "nFeature_RNA",
          "percent.mt",
          "percent.ribo"
        ),
        group.by = 'sample_id',
        pt.size=0)

VlnPlot(fib,
        features = c(
          "nCount_RNA",
          "nFeature_RNA",
          "percent.mt"
        ),
        group.by = 'sample_id',
        pt.size=0)

VlnPlot(obj,
        features = c(
          "nCount_RNA",
          "nFeature_RNA",
          "percent.mt"
        ),
        group.by = 'sample_id',
        pt.size=0)

Idents(fib) <- 'island'
deg_fib <- FindMarkers(
  fib,
  ident.1 = 'Fib2',
  ident.2 = 'Fib1',
  min.pct = 0.1
)
deg_fib$gene <- rownames(deg_fib)

deg_fib <- deg_fib %>%
  mutate(
    comparison = "All cells: Fib2 vs Fib1",
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC >= 0.5 ~"Higher in Fib2",
      p_val_adj < 0.05 & avg_log2FC <= -0.5 ~"Higher in Fib1",
      TRUE ~ "Not significant"
    )
  )

deg_fib %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC)>=0.5) %>%
           arrange(desc(abs(avg_log2FC))) %>%
           select(gene, avg_log2FC, pct.1, pct.2) %>%
           head(30)
          
deg_fib <- deg_fib %>%
  mutate(
    plot_p = pmax(p_val_adj, 1e-300),
    minus_log10_padj = -log10(plot_p)
  )

label_genes <- bind_rows(
  deg_fib %>%
    filter(significance == "Higher in Fib2") %>%
    arrange(p_val_adj, desc(avg_log2FC)) %>%
    slice_head(n=50),
  
  deg_fib %>%
    filter(significance == "Higher in Fib1") %>%
    arrange(p_val_adj, desc(avg_log2FC)) %>%
    slice_head(n=30)
) %>%
  distinct(gene, .keep_all = TRUE)

library(ggrepel)
p_vol <- ggplot(
  deg_fib,
  aes(x=avg_log2FC, y=minus_log10_padj)
)+
  geom_point(
    aes(color = significance),
    alpha = 0.65,
    size = 1.5
  )+
  geom_vline(
    xintercept = c(-0.5, 0.5),
    linetype = "dashed"
  )+
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  )+
  geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    min.segment.length = 0
  )+
  scale_color_manual(
    values = c(
      "Higher in Fib2" = "#D55E00",
      "Higher in Fib1" = "#0072B2",
      "Not significant" = "grey75"
    )
  )+theme_classic()+
  labs(x = "Average log2 fold change",
       y = expression(-log[10]("adjusted p-value")),
       color=NULL)
p_vol

genes <- c("COL1A1", "COL3A1", "COL5A2", "MFAP5", "TNXB", "MMP23B", "PRG4", "SAA1", "CREB3L1", "DCN", "LUM", "COL6A1")

DotPlot(fib, features = genes, group.by = 'island')+RotatedAxis()

FeaturePlot(fib,
            features = c("COL1A1", "COL3A1"),
            blend = TRUE)

fib_eth3 <- subset(fib, subset = sample_id == "CRSwNP_Eth_3")

VlnPlot(fib_eth3,
        features = c(
          "nCount_RNA",
          "nFeature_RNA",
          "percent.mt",
          "percent.ribo"
        ),
        group.by = 'island',
        pt.size=0, ncol = 2)

FeaturePlot(fib_eth3, features = c("COL1A1", "COL1A2", "POSTN", "CXCL14"))
FeaturePlot(fib_eth3, features = c("SAA1", "PRG4", "MFAP5", "TNXB", "ZEB2", "MMP23B"))
FeaturePlot(fib_eth3, features = c("RPLP0", "RPS3", "HSPA1A", "JUN", "FOS", "ATF3", "DDIT3"))
FeaturePlot(fib_eth3, features = c("RGS5", "MCAM", "CSPG4", "PDGFRB","THY1", "DES","ACTA2"))
FeaturePlot(fib_eth3, features = c("COL3A1", "COL6A1", "COL6A2", "DCN", "LUM", "FBLN1"))

deg_fib_eth3 <- FindMarkers(
  fib_eth3,
  ident.1 = 'Fib2',
  ident.2 = 'Fib1',
  min.pct = 0.1
)
deg_fib_eth3$gene <- rownames(deg_fib_eth3)

deg_fib_eth3 <- deg_fib_eth3 %>%
  mutate(
    comparison = "All cells: Fib2 vs Fib1",
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC >= 0.5 ~"Higher in Fib2",
      p_val_adj < 0.05 & avg_log2FC <= -0.5 ~"Higher in Fib1",
      TRUE ~ "Not significant"
    )
  )

deg_fib_eth3 %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC)>=0.5) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  select(gene, avg_log2FC, pct.1, pct.2) %>%
  head(30)

deg_fib_eth3 <- deg_fib %>%
  mutate(
    plot_p = pmax(p_val_adj, 1e-300),
    minus_log10_padj = -log10(plot_p)
  )

label_genes <- bind_rows(
  deg_fib_eth3 %>%
    filter(significance == "Higher in Fib2") %>%
    arrange(p_val_adj, desc(avg_log2FC)) %>%
    slice_head(n=30),
  
  deg_fib_eth3 %>%
    filter(significance == "Higher in Fib1") %>%
    arrange(p_val_adj, desc(avg_log2FC)) %>%
    slice_head(n=30)
) %>%
  distinct(gene, .keep_all = TRUE)

ecm.genes <- c("COL1A1", "COL1A2", "COL3A1", "COL5A2", "LUM", "CXCL14", "PRG4", "SAA1", "MFAP5", "CREB3L1", "KLF4")
fib <- AddModuleScore(fib, features = list(ecm.genes), name='ECM_score')
VlnPlot(fib, features = "ECM_score1")

DoHeatmap(fib, features = ecm.genes)

library(clusterProfiler)
library(org.Hs.eg.db)

genes_up <- deg_fib %>%
  filter(
    avg_log2FC > 1,
    p_val_adj < 0.05
  )%>%
  filter(!grepl("^RP[SL]", gene),
         !grepl("^MT-", gene))%>% pull(gene)

entrez <- bitr(
  genes_up,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ego <- enrichGO(
  gene = entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP"
)

dotplot(ego)

write.csv(deg_fib, file="CRS_Fibroblast_DEG.csv")
deg_fib2 <- deg_fib
deg_fib <- deg_fib %>%
  filter(
    avg_log2FC > 1,
    p_val_adj < 0.05
  )%>%
  filter(!grepl("^RP[SL]", gene),
         !grepl("^MT-", gene))

gene.df <- bitr(
  deg_fib$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)


deg2 <- left_join(
  deg_fib,
  gene.df,
  by = c("gene" = "SYMBOL")
)
geneList <- deg2$avg_log2FC
names(geneList) <- deg2$ENTREZID

geneList <- sort(geneList, decreasing = TRUE)

geneList <- geneList[!is.na(names(geneList))]
geneList <- geneList[!duplicated(names(geneList))]

gsea_go <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

library(enrichplot)
terms_keep <- c('extracellular matrix organization',
                'extracellular structure organization',
                'external encapsulating structure organization',
                'collagen fibril organization',
                'cell-substrate adhesion',
                'small GTPase-mediated signal transduction')

gsea_go_subset <- gsea_go
gsea_go_subset@result <- gsea_go_subset@result %>%
  filter(Description %in% terms_keep)
dotplot(gsea_go_subset, showCategory = 8)

gseaplot2(
  gsea_go,
  geneSetID = 1
)

library(msigdbr)

hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_list <- hallmark %>%
  dplyr::select(gs_name, entrez_gene)

gsea_h <- GSEA(
  geneList,
  TERM2GENE = hallmark_list,
  pvalueCutoff = 0.05
)

dotplot(gsea_h, showCategory = 20)
