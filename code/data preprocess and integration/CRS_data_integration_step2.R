rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc()
file.remove(".Rhistory")

setwd("~/project/results_nofilter_gene/")

.libPaths("~/R/intern")
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)

# Output directory
out_dir <- "~/project/results_nofilter_gene/integration_harmony"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- "~/project/results_nofilter_gene/integration_harmony/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

file.dir <- '~/project/results_nofilter_gene/CRS_batch_outputs/afterQC'
rdata_files <- list.files(
  path = file.dir,
  pattern = "\\.RData$",
  full.names = TRUE
)

# Check files
length(rdata_files)
print(rdata_files)

# ----Load all Seurat objects into a list----
set.seed(42)
seurat_list <- list()

for (f in rdata_files) {
  
  message("Loading: ", f)
  
  # Load into a temporary environment to avoid overwriting objects
  e <- new.env()
  load(f, envir = e)
  obj <- e$seurat_obj_qc
  
  # Get sample_id from metadata
  sample_name <- unique(as.character(obj$sample_id))
  message("Sample ID: ", sample_name)
  message("Number of cells: ", ncol(obj))
  
  # Rename cell barcodes to avoid duplicated cell names after merge
  obj <- RenameCells(
    object = obj,
    add.cell.id = sample_name
  )
  
  # Store in list
  seurat_list[[sample_name]] <- obj
}

# Check number of cells per sample
cell_counts <- sapply(seurat_list, ncol)
print(cell_counts)

#----Merge all samples----
combined <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1]
)
combined
# total number of cells: 65,913

# Check duplicated cell names
dup_n <- anyDuplicated(colnames(combined))
print(dup_n)

# Check sample distribution
print(table(combined$sample_id))

# Optional: check condition if you already have it
if ("condition" %in% colnames(combined@meta.data)) {
  print(table(combined$condition))
}

#----tandard workflow before batch correction----
DefaultAssay(combined) <- "RNA"

combined <- NormalizeData(
  combined,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

combined <- FindVariableFeatures(
  combined,
  selection.method = "vst",
  nfeatures = 3000
)

combined <- ScaleData(
  combined,
  features = VariableFeatures(combined)
)

combined <- RunPCA(
  combined,
  features = VariableFeatures(combined),
  npcs = 50
)

# Save elbow plot
p_elbow <- ElbowPlot(combined, ndims = 50)
p_elbow
ggsave(
  filename = file.path(fig_dir, "ElbowPlot_raw_PCA.pdf"),
  plot = p_elbow,
  width = 7,
  height = 5
)

#----UMAP before Harmony----
set.seed(42)
combined <- RunUMAP(
  combined,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap.raw"
)

p <- DimPlot(
  combined,
  reduction = "umap.raw",
  group.by = "sample_id",
  raster = TRUE
) +
  ggtitle("Before Harmony: sample_id")
p

ggsave(
  filename = file.path(fig_dir, "UMAP_before_Harmony_by_sample_id.pdf"),
  plot = p,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(fig_dir, "UMAP_before_Harmony_by_sample_id.png"),
  plot = p,
  width = 8,
  height = 6
)

p <- DimPlot(
  combined,
  reduction = "umap.raw",
  group.by = "condition",
  raster = TRUE
) +
  ggtitle("Before Harmony: condition")
p

ggsave(
  filename = file.path(fig_dir, "UMAP_before_Harmony_by_condition.pdf"),
  plot = p,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(fig_dir, "UMAP_before_Harmony_by_condition.png"),
  plot = p,
  width = 8,
  height = 6
)

saveRDS(combined, file = file.path(out_dir, "CRS_integrated_obj.rds")) 

#----Harmony batch correction by sample_id----
set.seed(42)
seu.harmony <- RunHarmony(
  combined,
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

p <- DimPlot(
  seu.harmony,
  reduction = "umap.harmony",
  group.by = "sample_id",
  raster = TRUE
) +
  ggtitle("After Harmony: sample_id")
p

ggsave(
  filename = file.path(fig_dir, "UMAP_after_Harmony_by_sample_id.pdf"),
  plot = p,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(fig_dir, "UMAP_after_Harmony_by_sample_id.png"),
  plot = p,
  width = 8,
  height = 6
)

p <- DimPlot(
  seu.harmony,
  reduction = "umap.harmony",
  group.by = "condition",
  raster = TRUE
) +
  ggtitle("After Harmony: condition")
p

ggsave(
  filename = file.path(fig_dir, "UMAP_after_Harmony_by_condition.pdf"),
  plot = p,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(fig_dir, "UMAP_before_Harmony_by_condition.png"),
  plot = p,
  width = 8,
  height = 6
)

seu.harmony <- FindNeighbors(
  seu.harmony,
  reduction = "harmony",
  dims = 1:30,
  graph.name = "harmony_snn"
)

for (res in c(0.05, 0.1, 0.2, 0.4, 0.5)) {
  seu.harmony <- FindClusters(
    seu.harmony,
    graph.name = "harmony_snn",
    resolution = res
  )
}

p <- DimPlot(seu.harmony, reduction = "umap.harmony", group.by = "harmony_snn_res.0.05")
ggsave(filename = file.path(fig_dir, "umap_res0_05.png"), plot = p,width = 8, height = 6)
ggsave(filename = file.path(fig_dir, "umap_res0_05.pdf"), plot = p,width = 8, height = 6)
p <- DimPlot(seu.harmony, reduction = "umap.harmony", group.by = "harmony_snn_res.0.1")
ggsave(filename = file.path(fig_dir, "umap_res0_1.png"), plot = p,width = 8, height = 6)
ggsave(filename = file.path(fig_dir, "umap_res0_1.pdf"), plot = p,width = 8, height = 6)
p <- DimPlot(seu.harmony, reduction = "umap.harmony", group.by = "harmony_snn_res.0.2")
ggsave(filename = file.path(fig_dir, "umap_res0_2.png"), plot = p,width = 8, height = 6)
ggsave(filename = file.path(fig_dir, "umap_res0_2.pdf"), plot = p,width = 8, height = 6)
p <- DimPlot(seu.harmony, reduction = "umap.harmony", group.by = "harmony_snn_res.0.4")
ggsave(filename = file.path(fig_dir, "umap_res0_4.png"), plot = p,width = 8, height = 6)
ggsave(filename = file.path(fig_dir, "umap_res0_4.pdf"), plot = p,width = 8, height = 6)
p <- DimPlot(seu.harmony, reduction = "umap.harmony", group.by = "harmony_snn_res.0.5")
ggsave(filename = file.path(fig_dir, "umap_res0_5.png"), plot = p,width = 8, height = 6)
ggsave(filename = file.path(fig_dir, "umap_res0_5.pdf"), plot = p,width = 8, height = 6)
p <- DimPlot(seu.harmony, reduction = "umap.harmony", group.by = "harmony_snn_res.0.6")
ggsave(filename = file.path(fig_dir, "umap_res0_6.png"), plot = p,width = 8, height = 6)
ggsave(filename = file.path(fig_dir, "umap_res0_6.pdf"), plot = p,width = 8, height = 6)
p <- DimPlot(seu.harmony, reduction = "umap.harmony", group.by = "harmony_snn_res.0.8")
ggsave(filename = file.path(fig_dir, "umap_res0_8.png"), plot = p,width = 8, height = 6)
ggsave(filename = file.path(fig_dir, "umap_res0_8.pdf"), plot = p,width = 8, height = 6)
p <- DimPlot(seu.harmony, reduction = "umap.harmony", group.by = "harmony_snn_res.1")
ggsave(filename = file.path(fig_dir, "umap_res1.png"), plot = p,width = 8, height = 6)
ggsave(filename = file.path(fig_dir, "umap_res1.pdf"), plot = p,width = 8, height = 6)

saveRDS(seu.harmony, file = file.path(out_dir, "CRS_harmony_obj.rds")) 

marker_list <- list(
  "Epithelial" = c("KRT7", "EPCAM", "KRT18"),
  "T-cells" = c("CD3E", "CD3D", "CD8A", "CD4", "TRAC"),
  "NK-cells" = c("KLRD1", "GNLY", "KLRB1"),
  "Myeloid" = c("LYZ", "ITGAM", "ITGAX"),
  "B-cells"= c("MS4A1", "CD79A", "CD19"),
  "pDCs" = c("CLEC4C", "IL3RA", "LILRA4"),
  "Fibroblasts" = c("COL1A1", "COL1A2", "DCN1"),
  "Plasma-Cells" = c("MZB1", "JCHAIN", "IGHG1", "IGKC", "IGHA1"),
  "Cycling-Cells" = c("MKI67", "STMN1", "TOP2A"),
  "Mast-Cells" = c("KIT", "GATA2", "TPSB2", "TPSAB1"),
  "Endothelial-Cells" = c("VWF", "PLVAP", "PECAM1"),
  "Erythroid" = c("HBA1", "HBA2", "ALAS2"),
  "Neutrophil" = c("CXCR2", "FCGR3B", "IL1R2")
)

#marker_vec <- unlist(marker_list, use.names = F)
p <- DotPlot(
  object = seu.harmony,
  features = marker_list,
  group.by = "harmony_snn_res.0.1",
  scale = T)
p

seu.harmony$celltype <- as.character(seu.harmony$harmony_snn_res.0.1)

seu.harmony$celltype[seu.harmony$harmony_snn_res.0.1 %in% c(1, 2, 8, 13)] <- 'Epithelial cells'
seu.harmony$celltype[seu.harmony$harmony_snn_res.0.1 == "0"] <- 'T/NK cells'
seu.harmony$celltype[seu.harmony$harmony_snn_res.0.1 == "3"] <- 'Myeloid'
seu.harmony$celltype[seu.harmony$harmony_snn_res.0.1 == "6"] <- 'B cells'
seu.harmony$celltype[seu.harmony$harmony_snn_res.0.1 == "12"] <- 'pDCs'
seu.harmony$celltype[seu.harmony$harmony_snn_res.0.1 %in% c(4, 10)] <- 'Fibroblasts'
seu.harmony$celltype[seu.harmony$harmony_snn_res.0.1 == "9"] <- 'Plasma cells'
seu.harmony$celltype[seu.harmony$harmony_snn_res.0.1 == "11"] <- 'Mast cells'
seu.harmony$celltype[seu.harmony$harmony_snn_res.0.1 == "5"] <- 'Endothelia cells'
seu.harmony$celltype[seu.harmony$harmony_snn_res.0.1 == "7"] <- 'Neutrophil'

p <- DimPlot(
  seu.harmony,
  reduction = "umap.harmony",
  group.by = "celltype",
  label = T,
  raster = TRUE
) +
  ggtitle("Major cell types")
p

ggsave(
  filename = file.path(fig_dir, "UMAP_major_celltypes.png"),
  plot = p,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(fig_dir, "UMAP_major_celltypes.pdf"),
  plot = p,
  width = 8,
  height = 6
)

FeaturePlot(seu.harmony, features = c("TPSAB1", "TPSB2", "CMA1","MRGPRX2"), reduction = "umap.harmony")
DefaultAssay(seu.harmony)
length(rownames(seu.harmony))
"MRGPRX2" %in% rownames(seu.harmony)

table(seu.harmony$celltype)
