library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(SpatialDecon)
library(Seurat)
library(patchwork)

set.seed(42)

## ============ read GeoMx data ============
# normalized & batch-corrected count matrix (gene x ROI), for downstream analysis
count_df <- read.csv("~/Downloads/countFile_normalized_and_batch_effect_corrected(ClusterBasedAnnotated).csv",
                     row.names = 1,
                     check.names = FALSE
)
row.names(count_df) <- count_df$TargetName
count_df$TargetName <- NULL
norm_author <- as.matrix(count_df)

# raw counts (gene x ROI), row=gene/probe(check negative probe),col=ROI/AOI segment ID
count_raw <- read.csv("~/Downloads/countFile_raw(ClusterBasedAnnotated).csv",
                      row.names = 1,
                      check.names = FALSE
)
raw_counts <- as.matrix(count_raw)

# ROI/AOI annotation
segment_meta <- read.csv("~/Downloads/sampleAnnoFile(ClusterBasedAnnotated).csv")
rownames(segment_meta) <- segment_meta$SegmentDisplayName

# sample alignment
common_samples <- Reduce(intersect, list(colnames(raw_counts), colnames(norm_author), segment_meta$SegmentDisplayName))
print(length(common_samples))

raw_counts   <- raw_counts[, common_samples, drop = FALSE]
norm_author  <- norm_author[, common_samples, drop = FALSE]
segment_meta <- segment_meta[common_samples, , drop = FALSE]

# identify negative probe
neg_idx <- grepl("negprobe", rownames(raw_counts), ignore.case = TRUE)
has_negprobe <- any(neg_idx)
cat(sprintf("%d negative probe is detected\n", sum(neg_idx)))

# perform Q3 normalization
q3_normalize <- function(counts_mat) {
  q3 <- apply(counts_mat, 2, function(x) quantile(x[x > 0], 0.75, na.rm = TRUE))
  q3_geomean <- exp(mean(log(q3[q3 > 0])))
  norm_factors <- q3 / q3_geomean
  list(norm = sweep(counts_mat, 2, norm_factors, "/"), factors = norm_factors)
}

raw_genes  <- raw_counts[!neg_idx, , drop = FALSE]
q3_out     <- q3_normalize(raw_genes)
norm_q3    <- q3_out$norm

if (has_negprobe) {
  neg_raw  <- raw_counts[neg_idx, , drop = FALSE]
  neg_norm <- sweep(neg_raw, 2, q3_out$factors, "/")
  
  norm_for_bg <- rbind(norm_q3, neg_norm)
  probepool <- rep(1, nrow(norm_for_bg))
  
  bg <- derive_GeoMx_background(norm = norm_for_bg, probepool = probepool, negnames = rownames(neg_norm))
  bg <- bg[rownames(norm_q3), , drop = FALSE]
} else {
  message("negative probes were not detected; a simplified background approximation (low-quantile estimation) is used.")
  bg <- matrix(rep(apply(norm_q3, 2, quantile, probs = 0.1), each = nrow(norm_q3)),
               nrow = nrow(norm_q3), dimnames = dimnames(norm_q3))
}

## ============ scRNA-seq data reference ============
scref <- readRDS("~/Library/CloudStorage/OneDrive-ArcusBiosciences,Inc/MZ/BFX intern/processed_obj/GSE235711_CRS_annotated_seurat_obj.rds")

sc_counts <- as.matrix(GetAssayData(scref, layer = "counts"))

cell_annots <- data.frame(
  cellID   = colnames(sc_counts),
  celltype = as.character(scref$celltype),
  stringsAsFactors = FALSE
)
print(table(cell_annots$celltype))

# generate cell profile matrix (gene x celltype)
custom_profile_mtx <- create_profile_matrix(
  mtx              = sc_counts,
  cellAnnots       = cell_annots,
  cellTypeCol      = "celltype",
  cellNameCol      = "cellID",
  matrixName       = "custom_scRNA_profile",
  outDir           = NULL,     
  normalize        = FALSE,   
  minCellNum       = 20,       
  minGenes         = 10,
  scalingFactor    = 5,       
  discardCellTypes = TRUE     
)

## ============ run SpatialDecon ============
shared_genes <- Reduce(intersect, list(rownames(norm_q3), rownames(bg), rownames(custom_profile_mtx)))

decon_res <- spatialdecon(
  norm        = norm_q3,
  bg          = bg,
  X           = custom_profile_mtx,
  raw         = raw_genes,    
  align_genes = TRUE
)

### prop_of_all: the proportion of the total deconvolved signal attributed to each cell type within an AOI.
### not equal to cell proportion
View(decon_res$prop_of_all)
View(decon_res$prop_of_nontumor)

## ============ save results ============
cell_fraction <- as.data.frame(t(decon_res$prop_of_all))
cell_fraction$SampleID <- rownames(cell_fraction)

cell_fraction_annot <- cell_fraction %>%
  left_join(segment_meta %>% mutate(SampleID = rownames(segment_meta)), by = "SampleID")

write.csv(cell_fraction_annot, "GeoMx_deconvolution_cell_proportions.csv", row.names = FALSE)

## ============ Visualization ============
prop.df <- as.data.frame(t(decon_res$prop_of_all))
prop.df$ROI <- rownames(prop.df)

prop.df <- left_join(
  prop.df,
  segment_meta,
  by = c("ROI" = "SegmentDisplayName")
)

plot_marker <- function(gene){
  
  df <- data.frame(
    ROI = colnames(norm_q3),
    Expression = as.numeric(norm_q3[gene,]),
    SegmentLabel = segment_meta[colnames(norm_q3),
                                "SegmentLabel"]
  )
  
  comparisons <- list(
    c("EPI", "IMM"),
    c("EPI", "MAC"),
    c("IMM", "MAC")
  )
  
  ggplot(
    df,
    aes(
      SegmentLabel,
      Expression,
      fill = SegmentLabel
    )
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
      width = 0.20,
      alpha = 0.5,
      size = 1.5
    ) +
    stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      label = "p.signif"
    ) +
    theme_classic() +
    ggtitle(gene)
}

## pie chart: overall cell composition
rm <- rowMeans(
  decon_res$prop_of_all,
  na.rm = TRUE
)

pie(
  rm,
  col = rainbow(length(rm)),
  main = "Average SpatialDecon Cell-Type Composition"
)

legend(
  "topright",
  legend = names(rm),
  fill = rainbow(length(rm)),
  cex = 0.8
)

pct <- round(rm / sum(rm) * 100, 1)

pie(
  rm,
  col = rainbow(length(rm)),
  labels = paste0(
    names(rm),
    "\n",
    pct,
    "%"
  ),
  main = "Average SpatialDecon Cell-Type Composition"
)

prop.df <- as.data.frame(t(decon_res$prop_of_all))
prop.df$ROI <- rownames(prop.df)

prop.df <- dplyr::left_join(
  prop.df,
  segment_meta,
  by = c("ROI" = "SegmentDisplayName")
)

seg.comp <- prop.df %>%
  dplyr::group_by(SegmentLabel) %>%
  dplyr::summarise(
    across(
      c(
        Myeloid,
        Fibroblasts,
        `Epithelial cells`,
        `T/NK cells`,
        Neutrophil,
        `Endothelial cells`,
        `B cells`,
        `Mast cells`,
        `Plasma cells`,
        pDCs
      ),
      mean,
      na.rm = TRUE
    )
  )

seg.long <- seg.comp %>%
  pivot_longer(
    -SegmentLabel,
    names_to = "CellType",
    values_to = "Proportion"
  )

seg.long2 <- seg.long %>%
  group_by(SegmentLabel) %>%
  arrange(desc(Proportion), .by_group = TRUE) %>%
  mutate(
    pct = paste0(round(Proportion * 100, 1), "%"),
    ypos = cumsum(Proportion) - 0.5 * Proportion
  )

ggplot(
  seg.long2,
  aes(
    x = "",
    y = Proportion,
    fill = CellType
  )
) +
  geom_col(width = 1) +
  geom_text(
    aes(
      y = ypos,
      label = paste0(CellType, "\n", pct)
    ),
    size = 3
  ) +
  coord_polar(theta = "y") +
  facet_wrap(~ SegmentLabel) +
  theme_void() +
  labs(
    title = "Average Cell-Type Composition by GeoMx Segment"
  )

library(dplyr)
library(ggplot2)

seg.long2 <- seg.long %>%
  group_by(SegmentLabel) %>%
  arrange(desc(Proportion), .by_group = TRUE) %>%
  mutate(
    pct = round(Proportion * 100, 1),
    label = paste0(CellType, "\n", pct, "%"),
    ypos = cumsum(Proportion) - 0.5 * Proportion
  )


## ============ GeoMx mask validation ============
epi_markers <- c(
  "EPCAM",
  "KRT8",
  "KRT18"
  #"KRT19",
  #"KRT17"
)

immune_markers <- c(
  "PTPRC",
  "CD3D",
  "CD3E"
)

macro_markers <- c(
  "CD68",
  "C1QA",
  "C1QB"
)

mast_markers <- c(
  "TPSAB1",
  "CPA3",
  "CMA1"
)


## Plot marker panels
wrap_plots(
  lapply(epi_markers, plot_marker),
  ncol = 3
)

wrap_plots(
  lapply(immune_markers, plot_marker),
  ncol = 3
)

wrap_plots(
  lapply(macro_markers, plot_marker),
  ncol = 3
)

wrap_plots(
  lapply(mast_markers, plot_marker),
  ncol = 3
)

## ============ SpatialDecon validation ============
## check whether the deconvolution results is reasonable

## Myeloid validation
myeloid_genes <- c(
  "CD68",
  "CD163",
  "AIF1",
  "TYROBP",
  "FCER1G",
  "C1QA",
  "C1QB"
)

## examine whether AOIs with higher inferred myeloid proportitons also have higher CD68 expression
myeloid_cor <- lapply(
  myeloid_genes,
  function(g){
    
    test <- cor.test(
      decon_res$prop_of_all["Myeloid", ],
      norm_q3[g, ],
      method = "spearman"
    )
    
    data.frame(
      Gene = g,
      rho = unname(test$estimate),
      pvalue = test$p.value
    )
  }
) %>% bind_rows()

myeloid_cor %>%
  arrange(desc(rho))

## Mast-cell validation
mast_validation_genes <- c(
  "TPSAB1",
  "CPA3",
  "CMA1",
  "KIT",
  "HPGDS"
)

mast_cor <- lapply(
  mast_validation_genes,
  function(g){
    
    test <- cor.test(
      decon_res$prop_of_all["Mast cells", ],
      norm_q3[g, ],
      method = "spearman"
    )
    
    data.frame(
      Gene = g,
      rho = unname(test$estimate),
      pvalue = test$p.value
    )
  }
) %>% bind_rows()

mast_cor %>%
  arrange(desc(rho))

## Average cell-type composition
avg_comp <- data.frame(
  CellType = rownames(decon_res$prop_of_all),
  MeanProp = rowMeans(
    decon_res$prop_of_all,
    na.rm = TRUE
  )
)

avg_comp %>%
  arrange(desc(MeanProp))

## ============  Segment enrichment ============ 
## calculate the median myeloid prop within each segment
aggregate(
  Myeloid ~ SegmentLabel,
  prop.df,
  median
)

kruskal.test(
  Myeloid ~ SegmentLabel,
  data = prop.df
)

aggregate(
  `Mast cells` ~ SegmentLabel,
  prop.df,
  median
)

kruskal.test(
  `Mast cells` ~ SegmentLabel,
  data = prop.df
)

aggregate(
  `Epithelial cells` ~ SegmentLabel,
  prop.df,
  median
)

kruskal.test(
  `Epithelial cells` ~ SegmentLabel,
  data = prop.df
)

## Boxplots
ggplot(
  prop.df,
  aes(
    SegmentLabel,
    Myeloid,
    fill = SegmentLabel
  )
) +
  geom_boxplot() +
  theme_classic()+
  ggtitle("Myeloid proportion across GeoMx segments")+
  stat_compare_means(method = "kruskal.test")

ggplot(
  prop.df,
  aes(
    SegmentLabel,
    `Epithelial cells`,
    fill = SegmentLabel
  )
) +
  geom_boxplot() +
  theme_classic()+
  ggtitle("Epithelial cell proportion across GeoMx segments")+
  stat_compare_means(method = "kruskal.test")

ggplot(
  prop.df,
  aes(
    SegmentLabel,
    `Mast cells`,
    fill = SegmentLabel
  )
) +
  geom_boxplot() +
  theme_classic()+
  ggtitle("Mast proportion across GeoMx segments")+
  stat_compare_means(method = "kruskal.test")


## ============ MRGPRX2 evaluation ============
mrg_genes <- c(
  "TPSAB1",
  "CPA3",
  "CMA1",
  "KIT",
  "HPGDS"
)

mrg_cor <- lapply(
  mrg_genes,
  function(g){
    
    test <- cor.test(
      norm_q3["MRGPRX2", ],
      norm_q3[g, ],
      method = "spearman"
    )
    
    data.frame(
      Gene = g,
      rho = unname(test$estimate),
      pvalue = test$p.value
    )
  }
) %>% bind_rows()

mrg_cor %>%
  arrange(desc(rho))

## MRGPRX2 vs mast abundance
cor.test(
  decon_res$prop_of_all["Mast cells", ],
  norm_q3["MRGPRX2", ],
  method = "spearman"
)

mrg.df <- data.frame(
  ROI = colnames(norm_q3),
  MRGPRX2 = as.numeric(norm_q3["MRGPRX2", ])
)

mrg.df <- left_join(
  mrg.df,
  segment_meta,
  by = c("ROI" = "SegmentDisplayName")
)

ggplot(
  mrg.df,
  aes(
    x = Working_group,
    y = MRGPRX2,
    fill = Working_group
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  theme_classic() +
  labs(
    title = "MRGPRX2 expression across disease groups",
    x = "",
    y = "Q3-normalized MRGPRX2 expression"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

mrg.df <- data.frame(
  ROI = colnames(norm_q3),
  MRGPRX2 = as.numeric(norm_q3["MRGPRX2", ])
)

mrg.df <- left_join(
  mrg.df,
  segment_meta,
  by = c("ROI" = "SegmentDisplayName")
)

ggplot(
  mrg.df,
  aes(
    x = SegmentLabel,
    y = MRGPRX2,
    fill = SegmentLabel
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  facet_wrap(~ Working_group, scales = "free_y") +
  theme_classic() +
  labs(
    title = "MRGPRX2 expression across GeoMx segments within each disease group",
    x = "Segment",
    y = "Q3-normalized MRGPRX2"
  )

mrg.df <- data.frame(
  ROI = colnames(norm_q3),
  MRGPRX2 = as.numeric(norm_q3["MRGPRX2", ])
)

mrg.df <- left_join(
  mrg.df,
  segment_meta,
  by = c("ROI" = "SegmentDisplayName")
)

# Remove "ImmuneTissue" samples
mrg.df2 <- subset(
  mrg.df,
  Working_group != "ImmuneTissue"
)

# Replace CRSwNP_UNC with CRSwNP
mrg.df2$Working_group[
  mrg.df2$Working_group == "CRSwNP_UNC"
] <- "CRSwNP"

library(ggpubr)

ggplot(
  mrg.df2,
  aes(
    x = Working_group,
    y = MRGPRX2,
    fill = Working_group
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  facet_wrap(~ SegmentLabel) +
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

## ============ MRGPRX2 signature genes ============
x2_genes <- c('ARL5B', 'CREM', 'LDLR', 'MAPK6', 'NR4A1', 'RRM2', 
              'SGK1', 'SOCS3', 'TBX20', 'PNOC', 'IL13', 'SEC14L2', 'ARC')

x2_genes <- intersect(
  x2_genes,
  rownames(norm_q3)
)

expr.df <- as.data.frame(t(norm_q3[x2_genes, ]))
expr.df$ROI <- rownames(expr.df)

expr.df <- left_join(
  expr.df,
  segment_meta,
  by = c("ROI" = "SegmentDisplayName")
)

expr.long <- expr.df %>%
  pivot_longer(
    cols = all_of(x2_genes),
    names_to = "Gene",
    values_to = "Expression"
  )

ggplot(
  expr.long,
  aes(
    x = Working_group,
    y = Expression,
    fill = Working_group
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  facet_grid(
    Gene ~ SegmentLabel,
    scales = "free_y"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

heat.df <- expr.long %>%
  group_by(
    Working_group,
    SegmentLabel,
    Gene
  ) %>%
  summarise(
    MeanExpr = mean(Expression),
    .groups = "drop"
  )

heat.mat <- heat.df %>%
  unite(
    Group,
    Working_group,
    SegmentLabel,
    sep = "_"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = MeanExpr
  )

heat.mat <- as.data.frame(heat.mat)
rownames(heat.mat) <- heat.mat$Gene

heat.mat$Gene <- NULL
heat.mat <- as.matrix(heat.mat)

pheatmap::pheatmap(
  heat.mat,
  scale = "row"
)

genes_to_plot <- c(
  "MRGPRX2",
  "IL13",
  "CREM",
  "PNOC",
  "CPA3",
  "TPSAB1"
)

genes_to_plot <- intersect(
  genes_to_plot,
  rownames(norm_q3)
)

library(dplyr)
library(tidyr)
library(ggplot2)

expr.df <- as.data.frame(
  t(norm_q3[genes_to_plot, ])
)

expr.df$ROI <- rownames(expr.df)

expr.df <- left_join(
  expr.df,
  segment_meta,
  by = c("ROI" = "SegmentDisplayName")
)

expr.long <- expr.df %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "Gene",
    values_to = "Expression"
  )

ggplot(
  expr.long,
  aes(
    x = Working_group,
    y = Expression,
    fill = Working_group
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.15,
    alpha = 0.4,
    size = 0.8
  ) +
  facet_grid(
    Gene ~ SegmentLabel,
    scales = "free_y"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = "X2-related gene expression across working groups and segments",
    x = NULL,
    y = "Q3-normalized expression"
  )


x2_score <- colMeans(
  norm_q3[x2_genes, ],
  na.rm = TRUE
)

x2.df <- data.frame(
  ROI = colnames(norm_q3),
  X2_score = x2_score
)

x2.df <- left_join(
  x2.df,
  segment_meta,
  by = c("ROI" = "SegmentDisplayName")
)

ggplot(
  x2.df,
  aes(
    Working_group,
    X2_score,
    fill = Working_group
  )
) +
  geom_boxplot() +
  facet_wrap(~ SegmentLabel) +
  theme_classic() +
  labs(
    title = "MRGPRX2-associated mast cell signature"
  )

library(ggplot2)

df <- data.frame(
  MRGPRX2 = as.numeric(norm_q3["MRGPRX2", ]),
  CMA1    = as.numeric(norm_q3["CMA1", ])
)

cor.test(
  df$MRGPRX2,
  df$CMA1,
  method = "spearman"
)

ggplot(df, aes(x = MRGPRX2, y = CMA1)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  theme_classic() +
  labs(
    title = "MRGPRX2 vs CMA1",
    x = "MRGPRX2 expression",
    y = "CMA1 expression"
  )