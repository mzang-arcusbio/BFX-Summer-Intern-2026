rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc()

setwd("~/project/results_nofilter_gene/")
fig_dir <- "~/project/results_nofilter_gene/cellchat/"

.libPaths("~/R/intern")
options(future.globals.maxSize = 32 * 1024^5)

library(CellChat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)
library(purrr)
library(ComplexHeatmap)

# ----CellChat Comparative Analysis Across CRS Conditions----

cellchat.list <- readRDS('~/project/results_nofilter_gene/integration_harmony/cellchat/cellchat_obj_by_condition_celltype3.rds')

cellchat.CL   <- cellchat.list$Control
cellchat.CSE  <- cellchat.list$CRSsNP_Eth
cellchat.CWE  <- cellchat.list$CRSwNP_Eth
cellchat.CWNP <- cellchat.list$CRSwNP_NP

object.list <- list(
  Control    = cellchat.CL,
  CRSsNP_Eth = cellchat.CSE,
  CRSwNP_Eth = cellchat.CWE,
  CRSwNP_NP  = cellchat.CWNP
)

## ----Merge CellChat objects for condition comparison----
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)


## ----Compare interaction number and interaction strength----
png(
  filename = paste0(fig_dir, "compareInteractions.png"),
  width = 1500,
  height = 1000,
  res = 300
)

gg1 <- compareInteractions(
  cellchat,
  group = c(1,2,3,4),
  angle.x = 45,
  show.legend = FALSE
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gg2 <- compareInteractions(
  cellchat,
  group = c(1,2,3,4),
  measure = "weight",
  show.legend = FALSE
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gg1 + gg2

dev.off()

## ----Differential network comparison----
# CRSwNP_NP vs Control
par(mfrow = c(1,2), xpd = TRUE)

netVisual_diffInteraction(
  cellchat,
  comparison = c(1,4),
  weight.scale = TRUE
)

netVisual_diffInteraction(
  cellchat,
  comparison = c(1,4),
  weight.scale = TRUE,
  measure = "weight"
)

# CRSwNP_Eth vs Control
par(mfrow = c(1,2), xpd = TRUE)

netVisual_diffInteraction(
  cellchat,
  comparison = c(1,3),
  weight.scale = TRUE
)

netVisual_diffInteraction(
  cellchat,
  comparison = c(1,3),
  weight.scale = TRUE,
  measure = "weight"
)

# CRSwNP_Eth vs CRSwNP_NP
par(mfrow = c(1,2), xpd = TRUE)

netVisual_diffInteraction(
  cellchat,
  comparison = c(3,4),
  weight.scale = TRUE
)

netVisual_diffInteraction(
  cellchat,
  comparison = c(3,4),
  weight.scale = TRUE,
  measure = "weight"
)

## ----Heatmap comparison----
gg1 <- netVisual_heatmap(
  cellchat,
  comparison = c(1,4)
)

gg2 <- netVisual_heatmap(
  cellchat,
  comparison = c(1,4),
  measure = "weight"
)

gg1 + gg2

## ----Compute centrality once----
object.list <- lapply(
  object.list,
  function(x){
    netAnalysis_computeCentrality(
      x,
      slot.name = "netP"
    )
  }
)


## ----Signaling role scatter plot for all cell types----
png(
  filename = paste0(fig_dir, "signalingrole_scatter.png"),
  width = 3200,
  height = 2400,
  res = 300
)

num.link <- sapply(
  object.list,
  function(x){
    rowSums(x@net$count) +
      colSums(x@net$count) -
      diag(x@net$count)
  }
)

weight.MinMax <- c(
  min(num.link),
  max(num.link)
)

gg <- list()

for(i in seq_along(object.list)){
  
  gg[[i]] <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    title = names(object.list)[i],
    weight.MinMax = weight.MinMax
  )
}

patchwork::wrap_plots(gg)

dev.off()

## ----Subset myeloid populations----
myeloid.celltypes <- c(
  "pDCs",
  "Monocytes",
  "Macrophages",
  "Mono_Macro",
  "Neutrophi"
)

object.list.sub <- lapply(
  object.list,
  function(x){
    
    x.sub <- subsetCellChat(
      x,
      idents.use = myeloid.celltypes
    )
    
    x.sub <- netAnalysis_computeCentrality(
      x.sub,
      slot.name = "netP"
    )
    
    return(x.sub)
  }
)

## ----Extract signaling role summary----

extractSignalingRoleDF <- function(object){
  
  centr <- object@netP$centr
  
  outgoing <- matrix(
    0,
    nrow = nlevels(object@idents),
    ncol = length(centr)
  )
  
  incoming <- matrix(
    0,
    nrow = nlevels(object@idents),
    ncol = length(centr)
  )
  
  dimnames(outgoing) <- list(
    levels(object@idents),
    names(centr)
  )
  
  dimnames(incoming) <- dimnames(outgoing)
  
  for(i in seq_along(centr)){
    
    outgoing[,i] <- centr[[i]]$outdeg
    incoming[,i] <- centr[[i]]$indeg
  }
  
  interaction.count <- rowSums(object@net$count) +
    colSums(object@net$count) -
    diag(object@net$count)
  
  data.frame(
    celltype = rownames(outgoing),
    outgoing_strength = rowSums(outgoing),
    incoming_strength = rowSums(incoming),
    interaction_count = interaction.count
  )
}

## ----Build combined dataframe----
df.all <- do.call(
  rbind,
  lapply(
    names(object.list),
    function(nm){
      
      df <- extractSignalingRoleDF(
        object.list[[nm]]
      )
      
      df$condition <- nm
      
      df
    }
  )
)

## ----Function to generate signaling role scatter----
plot_signaling_role <- function(df.all,
                                celltypes.use,
                                title){
  
  df.sel <- df.all[
    df.all$celltype %in% celltypes.use,
  ]
  
  weight.MinMax <- c(
    min(df.all$interaction_count),
    max(df.all$interaction_count)
  )
  
  ggplot(
    df.sel,
    aes(
      x = outgoing_strength,
      y = incoming_strength,
      colour = condition
    )
  ) +
    geom_point(
      aes(size = interaction_count),
      alpha = 0.75
    ) +
    geom_text_repel(
      aes(label = celltype),
      show.legend = FALSE,
      max.overlaps = Inf
    ) +
    scale_size_continuous(
      range = c(3,10),
      limits = weight.MinMax,
      name = "Interaction count"
    ) +
    labs(
      x = "Outgoing interaction strength",
      y = "Incoming interaction strength",
      title = title,
      colour = "Condition"
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

p.mast <- plot_signaling_role(
  df.all,
  c("MCT","MCTC"),
  "Mast Cells"
)

p.myeloid <- plot_signaling_role(
  df.all,
  c(
    "pDCs",
    "Monocytes",
    "Macrophages",
    "Mono_Macro",
    "Neutrophi"
  ),
  "Myeloid Cells"
)

p.bt <- plot_signaling_role(
  df.all,
  c(
    "T/NK cells",
    "B cells"
  ),
  "Lymphocytes"
)

p.structural <- plot_signaling_role(
  df.all,
  c(
    "Fibroblasts",
    "Epithelial cells",
    "Endothelial cells"
  ),
  "Structural Cells"
)

(p.mast | p.myeloid) /
  (p.bt | p.structural)


## -----Functional signaling similarity----
cellchat <- computeNetSimilarityPairwise(
  cellchat,
  type = "functional"
)

cellchat <- netEmbedding(
  cellchat,
  type = "functional"
)

cellchat <- netClustering(
  cellchat,
  type = "functional"
)

netVisual_embeddingPairwise(
  cellchat,
  type = "functional",
  label.size = 3.5
)

## ----Rank signaling pathways----
png(
  filename = paste0(fig_dir, "rankNet.png"),
  width = 3000,
  height = 3000,
  res = 300
)

gg1 <- rankNet(
  cellchat,
  mode = "comparison",
  comparison = c(1,2,3,4),
  measure = "weight",
  stacked = TRUE,
  do.stat = TRUE
)

gg2 <- rankNet(
  cellchat,
  mode = "comparison",
  comparison = c(1,2,3,4),
  measure = "weight",
  stacked = FALSE,
  do.stat = TRUE
)

gg1 + gg2

dev.off()

## ----Top pathways by total information flow----
df.rank <- rankNet(
  cellchat,
  mode = "comparison",
  comparison = c(1,2,3,4),
  measure = "weight",
  stacked = TRUE,
  return.data = TRUE
)

top.pathways <- df.rank %>%
  group_by(name) %>%
  summarise(
    total_flow = sum(contribution)
  ) %>%
  arrange(desc(total_flow)) %>%
  slice_head(n = 15) %>%
  pull(name)

rankNet(
  cellchat,
  mode = "comparison",
  comparison = c(1,2,3,4),
  measure = "weight",
  stacked = TRUE,
  signaling = top.pathways
)

## ----Signaling role heatmaps----
pathway.union <- Reduce(
  union,
  lapply(
    object.list,
    function(x) x@netP$pathways
  )
)

heatmaps <- lapply(
  seq_along(object.list),
  function(i){
    
    netAnalysis_signalingRole_heatmap(
      object.list[[i]],
      pattern = "outgoing",
      signaling = pathway.union,
      title = names(object.list)[i],
      width = 5,
      height = 6
    )
  }
)

draw(
  Reduce(`+`, heatmaps),
  ht_gap = unit(0.5, "cm")
)

## ----Mast cell to myeloid bubble plot----

png(
  filename = paste0(fig_dir, "Mast_Bubble.png"),
  width = 2400,
  height = 2400,
  res = 300
)

netVisual_bubble(
  cellchat,
  sources.use = c(5,6),
  targets.use = c(2,3,4),
  comparison = c(2,3,4),
  angle.x = 45,
  remove.isolate = TRUE
)

dev.off()

