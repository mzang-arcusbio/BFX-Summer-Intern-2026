library(UCell)
library(Seurat)

#----UCell----
signature <- list(
  X2_signature=c('ARL5B', 'CREM', 'LDLR', 'MAPK6', 'NR4A1', 'RRM2', 
                 'SGK1', 'SOCS3', 'TBX20', 'PNOC', 'IL13', 'SEC14L2', 'ARC')
)

obj <- AddModuleScore_UCell(obj, 
                            features=signature, name=NULL)

head(obj[[]])
