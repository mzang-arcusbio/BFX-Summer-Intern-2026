.libPaths("~/R/intern")
options(future.globals.maxSize = 32 * 1024^5)

library(CellChat)
library(Seurat)

obj <- readRDS('~/project/results_nofilter_gene/integration_harmony/GSE235711_CRS_annotated_seurat_obj.rds')
condition_list <- SplitObject(obj, split.by = "condition2")

cellchat.list <- lapply(condition_list, function(seu) {
  DefaultAssay(seu) <- "RNA"
  data.input <- GetAssayData(seu, slot = "data")
  meta <- seu@meta.data[, c("celltype3"), drop = FALSE]
  createCellChat(data.input, meta = meta, group.by = "celltype3")
})

cellchat.list <- lapply(cellchat.list, function(cellchat) {
  
  # set database
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB
  
  # subset signaling genes
  cellchat <- subsetData(cellchat)
  future::plan("multisession", workers = 4)
  
  # identify overexpressed genes/interactions
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # compute communication probability
  cellchat <- computeCommunProb(cellchat)
  
  # filter low-cell interactions
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # compute pathway level communication
  cellchat <- computeCommunProbPathway(cellchat)
  
  # aggregate network
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
})

saveRDS(cellchat.list, '~/project/results_nofilter_gene/integration_harmony/cellchat/cellchat_obj_by_condition_celltype3.rds')
