.libPaths("~/R/intern")
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)

# Output directory
out_dir <- "~/project/results/integration_harmony"
fig_dir <- "~/project/results/integration_harmony/figures"
marker_dir <- "~/project/results/integration_harmony/markers"

set.seed(42)
seu.harmony <- readRDS(file.path(out_dir, "CRS_harmony_obj.rds"))
DefaultAssay(seu.harmony) <- "RNA"
resolutions <- c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0)

for (res in resolutions) {
  
  Idents(seu.harmony) <- paste0("harmony_snn_res.", res)
  
  # Find markers
  markers <- FindAllMarkers(
    seu.harmony,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  
  # Save to CSV
  file_name <- paste0("markers_res_", res, ".csv")
  
  write.csv(
    markers,
    file = file.path(marker_dir, file_name),
    row.names = FALSE
  )
}