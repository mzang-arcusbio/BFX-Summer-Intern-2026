rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc()

setwd("~/project/results_nofilter_gene/")
fig_dir <- "~/project/results_nofilter_gene/integration_harmony/myeloid_subset/"

.libPaths("~/R/intern")
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(patchwork)

#----Myeloid cell subclustering----
#obj <- readRDS('./integration_harmony/CRS_annotated_obj.rds')
obj <- readRDS('./integration_harmony/GSE235711_CRS_annotated_seurat_obj.rds')

set.seed(42)
myeloid.obj <- subset(obj, subset = celltype=="Myeloid")
myeloid.obj <- NormalizeData(myeloid.obj)
myeloid.obj <- FindVariableFeatures(myeloid.obj, selection.method = "vst", nfeatures = 3000)
myeloid.obj <- ScaleData(myeloid.obj, features = rownames(myeloid.obj))
myeloid.obj <- RunPCA(myeloid.obj, features = VariableFeatures(object = myeloid.obj), npcs = 50)

# Harmony batch correction by sample_id
set.seed(42)
seu.harmony <- RunHarmony(
  myeloid.obj,
  "sample_id"
)

# Check reductions
print(Reductions(seu.harmony))

#UMAP and clustering after Harmony
set.seed(42)
seu.harmony <- RunUMAP(
  seu.harmony,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap.harmony"
)

seu.harmony <- FindNeighbors(
  seu.harmony,
  dims = 1:30
)

set.seed(42)
seu.harmony <- FindClusters(
  seu.harmony,
  resolution = 0.05
)

p <- DimPlot(
  seu.harmony,
  reduction = "umap.harmony"
) 
p

marker_list <- list(
  "Monocytes" = c("VCAN", "FCN1", "CD14", "FCGR3A"),
  "Macrophages" = c("C1QA", "C1QB", "C1QC"),
  "cDC1" = c("CLEC9A"),
  "cDC2" = c("CD1C"),
  "pDC" = c("CLEC4C", "IL3RA"),
  "Cycling" = c("MKI67"),
  "mregDC" = c("LAMP3", "CCR7"),
  "T-cells" = c("CD3E", "CD3D", "CD8A", "CD4", "TRAC"),
  "NK-cells" = c("KLRD1", "GNLY", "KLRB1")
)

marker_list <- list(
  "Monocytes" = c("VCAN", "FCN1", "CD14", "FCGR3A"),
  "Macrophages" = c("C1QA", "C1QB", "C1QC")
  #"Others"=c("CD16", "", "MRC1", "VEGFA", "CCL4L2")
)

p <- DotPlot(
  object = seu.harmony,
  features = marker_list,
  #group.by = "harmony_snn_res.0.1",
  scale = F)
p

ggsave(
  filename = file.path(fig_dir, "Myeloid_marker_dotplot.png"),
  plot = p,
  width = 7.5,
  height = 4.5
)

p <- FeaturePlot(seu.harmony, 
            reduction = "umap.harmony",
            #features = c("C1QA", "C1QB", "C1QC", "VCAN", "FCN1", "CD14")
            features = c(
              "FCN1", "S100A8", "S100A9", "CCR2",
              "VCAN", "FCN1", "CD14", "FCGR3A", "C1QA", "C1QB", "C1QC", "APOE", "CD68")
                         #, "MRC1", "VEGFA", "CCL4L2")
            )
p
ggsave(
  filename = file.path(fig_dir, "Myeloid_marker_featurplot.png"),
  plot = p,
  width = 12,
  height = 10
)

de_0_vs_2 <- FindMarkers(
  seu.harmony,
  ident.1 = "0",
  ident.2 = "2"
)

head(de_0_vs_2)

table(Idents(seu.harmony))

set.seed(42)
obj0 <- subset(seu.harmony, idents = "0")
obj0 <- FindVariableFeatures(obj0)
obj0 <- ScaleData(obj0)
obj0 <- RunPCA(obj0)
obj0 <- FindNeighbors(obj0, reduction = "pca", dims = 1:15)
obj0 <- FindClusters(obj0, resolution = 0.05)
obj0 <- RunUMAP(obj0, reduction = "pca", dims = 1:15)
DimPlot(
  obj0,
  reduction = "umap"
) 

DotPlot(
  obj0,
  features = c(
    "FCN1","S100A8","S100A9","LYZ",
    "FCER1G","CTSB","IFI30",
    "CD68","C1QA","C1QB","C1QC",
    "APOE","APOC1","MSR1"
  )
) + RotatedAxis()

table(Idents(obj0))
FindMarkers(obj0, ident.1 = "0", ident.2 = "1")
FindMarkers(obj0, ident.1 = "1", ident.2 = "0")

FeaturePlot(
  obj0,
  features = c(
    "CD3D",
    "TRAC",
    "FCN1",
    "VCAN",
    "CD14",
    "C1QA"
  ),
  reduction = 'umap'
)

table(seu.harmony$seurat_clusters)
seu.harmony$celltype3 <- as.character(seu.harmony$seurat_clusters)
seu.harmony$celltype3[seu.harmony$seurat_clusters == "0"] <- 'Mono_Macro'
seu.harmony$celltype3[seu.harmony$seurat_clusters == "1"] <- 'Monocytes'
seu.harmony$celltype3[seu.harmony$seurat_clusters == "2"] <- 'Macrophages'

obj$celltype3 <- obj$celltype2
obj$celltype3[colnames(seu.harmony)] <- seu.harmony$celltype3
table(obj$celltype3)

myeloid.obj <- subset(obj, celltype %in% c('Myeloid', 'Neutrophil', 'pDCs'))

p <- DimPlot(
  myeloid.obj,
  reduction = "umap.harmony",
  group.by = "celltype3",
  label = T,
  raster = TRUE
) +
  ggtitle("Myeloid cells")

ggsave(
  filename = file.path(fig_dir, "Myeloid_umap.png"),
  plot = p,
  width = 7,
  height = 5.5
)

#----check X2 activation signature genes----
X2.signature <- c('ARL5B', 'CREM', 'LDLR', 'MAPK6', 'NR4A1', 'RRM2', 
                  'SGK1', 'SOCS3', 'TBX20', 'PNOC', 'IL13', 'SEC14L2', 'ARC')

p <- VlnPlot(myeloid.obj, features = X2.signature, group.by = 'celltype3', ncol = 5)
ggsave(
  filename = file.path(fig_dir, "Myeloid_X2_vlnplot.png"),
  plot = p,
  width = 15,
  height = 10
)

plots <- VlnPlot(myeloid.obj, features = "X2_signature", split.by = "condition2", group.by = "celltype3",
                 cols = c("Control"="#ffadad",
                        "CRSsNP_Eth"="#b6e2d3",
                        "CRSwNP_Eth"="#809bce", 
                        "CRSwNP_NP" = "#fbaf87"), pt.size = 0, combine = FALSE)
wrap_plots(plots = plots)

myeloid.obj@meta.data %>% group_by(celltype3, condition2) %>% summarise(n=n())
myeloid.obj@meta.data %>% group_by(celltype, condition2) %>% summarise(n=n())

df_summary <- myeloid.obj@meta.data %>% 
  group_by(celltype3, orig.ident, condition2) %>% 
  summarise(n = n(), .groups = "drop")

df_summary <- df_summary %>%
  group_by(condition2, orig.ident) %>%
  mutate(pct = n / sum(n) * 100)

p <- ggplot(df_summary, aes(x = orig.ident, y = pct, fill = celltype3)) +
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
  filename = file.path(fig.dir, "myeloid_composition_bycondition.png"),
  plot = p,
  width = 8,
  height = 4)


X2.signature <- c('ARL5B', 'CREM', 'LDLR', 'MAPK6', 'NR4A1', 'RRM2', 
                  'SGK1', 'SOCS3', 'TBX20', 'PNOC', 'IL13', 'SEC14L2', 'ARC')

expr <- GetAssayData(myeloid.obj, assay = 'RNA', slot = 'data')[X2.signature, ] 
pct_expr <- sapply(unique(myeloid.obj$celltype3), function(ct){
  cells <- colnames(myeloid.obj)[myeloid.obj$celltype3 == ct]
  Matrix::rowMeans(expr[, cells] > 0) * 100
})
pct_expr <- as.data.frame(pct_expr)
pct_expr$gene <- rownames(pct_expr)
pct_expr

