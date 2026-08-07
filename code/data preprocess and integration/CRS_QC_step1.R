rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc()

setwd("~/project/results_nofilter_gene/")

.libPaths("~/R/intern")
library(Seurat)
library(dplyr)
library(DoubletFinder)
library(gt)

set.seed(42)
# ----global output----
out_root <- file.path(getwd(), "CRS_batch_outputs2")
dir.create(out_root, showWarnings = FALSE)

## ----QC summary dir----
summary_dir <- file.path(out_root, "QC_summary_tables2")
dir.create(summary_dir, showWarnings = FALSE)

## ----Seurat objects dir----
raw_seurat_dir <- file.path(out_root, "raw_seurat2")
dir.create(raw_seurat_dir, showWarnings = FALSE)

## ----noDoublet dir----
seurat_noDB_dir <- file.path(out_root, "noDoublet2")
dir.create(seurat_noDB_dir, showWarnings = FALSE)

## ----afterQC dir----
seurat_afterQC_dir <- file.path(out_root, "afterQC2")
dir.create(seurat_afterQC_dir, showWarnings = FALSE)

## ----calculate QC----
add_qc_metrics <- function(seurat_obj, species = "human") {
  DefaultAssay(seurat_obj) <- "RNA"
  
  # mitochondria
  mt_pattern <- ifelse(species == "mouse", "^mt-", "^MT-")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  
  # ribosome
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPL|^RPS")
  
  # hemoglobin (human-oriented)
  hb.pattern <- "^HB[ABDEGMQZ]|^HBA|^HBB|^HBD|^HBE|^HBG|^HBM"
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = hb.pattern)
  
  seurat_obj
}

make_qc_summary <- function(md, sample_id, N_doublets = NA_integer_) {
  md2 <- md %>%
    mutate(
      nUMI  = nCount_RNA,
      nGene = nFeature_RNA,
      Log10GenePerUMI = log10(pmax(nGene, 1)) / log10(pmax(nUMI, 1))
    )
  
  tibble(
    sample_id = sample_id,
    N_total_cells = nrow(md2),
    N_doublets = as.integer(N_doublets),
    
    pct_nUMI_le_500  = mean(md2$nUMI  <= 500, na.rm = TRUE)  * 100,
    pct_nUMI_le_800  = mean(md2$nUMI  <= 800, na.rm = TRUE)  * 100,
    pct_nGene_le_250 = mean(md2$nGene <= 250, na.rm = TRUE) * 100,
    
    pct_MT_gt_10 = mean(md2$percent.mt > 10, na.rm = TRUE) * 100,
    pct_MT_gt_15 = mean(md2$percent.mt > 15, na.rm = TRUE) * 100,
    pct_MT_gt_20 = mean(md2$percent.mt > 20, na.rm = TRUE) * 100,
    pct_MT_gt_25 = mean(md2$percent.mt > 25, na.rm = TRUE) * 100,
    
    pct_Log10GenePerUMI_le_0_8 = mean(md2$Log10GenePerUMI <= 0.8, na.rm = TRUE) * 100,
    
    pct_ribo_gt_10 = mean(md2$percent.ribo > 10, na.rm = TRUE) * 100,
    pct_ribo_gt_15 = mean(md2$percent.ribo > 15, na.rm = TRUE) * 100,
    pct_ribo_gt_20 = mean(md2$percent.ribo > 20, na.rm = TRUE) * 100,
    pct_ribo_gt_25 = mean(md2$percent.ribo > 25, na.rm = TRUE) * 100,
    
  )
}

format_percent_display <- function(df) {
  df %>%
    mutate(
      across(
        where(is.numeric) & !c(N_total_cells, N_doublets),
        ~ sprintf("%.2f%%", .x)
      )
    )
}

# ----find all samples----
sample_dirs <- list.dirs("~/project/data", recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("_filtered$", basename(sample_dirs))]

all_initial_qc <- list()
all_afterDB_qc <- list()

# for loop:QC + DoubletFinder for each sample
# Check the gene exists (important!)
for (sd in sample_dirs) {

  folder_name <- basename(sd)
  sample_id <- sub("^[^_]+_(.+)_filtered$", "\\1", folder_name)
  message("\n============================")
  message("Processing: ", sample_id)
  message("============================")
  
  data <- Read10X(data.dir = sd)
  seurat_obj <- CreateSeuratObject(
    counts = data[["Gene Expression"]],
    project = sample_id,
    min.cells = 0,
    min.features = 0
  )
# Extract expression vector
expr <- GetAssayData(seurat_obj, slot = "counts")["MRGPRX2", ]

# Count expressing cells
n_expr_cells <- sum(expr > 0)

# Total cells
n_cells <- ncol(seurat_obj)

# Fraction
pct_expr <- n_expr_cells / n_cells * 100

# Print result
message("MRGPRX2+ cells: ", n_expr_cells, "/", n_cells,
        " (", round(pct_expr, 2), "%)")

}

grep("MRGPR", rownames(seurat_obj), value=T)
grep("X2", rownames(seurat_obj), value=T)

for (sd in sample_dirs) {
  folder_name <- basename(sd)
  sample_id <- sub("^[^_]+_(.+)_filtered$", "\\1", folder_name)
  message("\n============================")
  message("Processing: ", sample_id)
  message("============================")
  
  data <- Read10X(data.dir = sd)
  seurat_obj <- CreateSeuratObject(
    counts = data[["Gene Expression"]],
    project = sample_id,
    min.cells = 0,
    min.features = 0
  )
  
  save(seurat_obj, file = file.path(raw_seurat_dir, paste0(sample_id, "_raw_SeuratObject.RData")))
  
  # initial QC summary
  seurat_obj <- add_qc_metrics(seurat_obj, species = "human")
  initial_qc <- make_qc_summary(seurat_obj@meta.data, sample_id = sample_id, N_doublets = NA_integer_)
  all_initial_qc[[sample_id]] <- initial_qc
  
  # ---------- DoubletFinder pre-processing(LogNormalize) ----------
  sObject <- seurat_obj
  sObject <- NormalizeData(sObject)
  sObject <- FindVariableFeatures(sObject, selection.method = "vst", nfeatures = 2000)
  sObject <- ScaleData(sObject)
  sObject <- RunPCA(sObject)
  
  stdv <- sObject[["pca"]]@stdev
  percent_stdv <- (stdv / sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
  co2 <- sort(which((percent_stdv[-length(percent_stdv)] - percent_stdv[-1]) > 0.1), decreasing = TRUE)[1] + 1
  min_pc <- min(co1, co2)
  if (is.na(min_pc) || min_pc < 10) min_pc <- 10  
  set.seed(42)
  sObject <- RunUMAP(sObject, dims = 1:min_pc)
  sObject <- FindNeighbors(sObject, dims = 1:min_pc)
  sObject <- FindClusters(sObject, resolution = 0.5)
  
  sweep_list <- paramSweep(sObject, PCs = 1:min_pc, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats)
  optimal.pk <- bcmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    pull(pK) %>%
    as.character() %>%
    as.numeric()
  
  annotations <- sObject@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  
  # multiplet rate table
  multiplet_rates_10x <- data.frame(
    Multiplet_rate = c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
    Loaded_cells   = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
    Recovered_cells= c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
  )
  
  multiplet_rate <- multiplet_rates_10x %>%
    filter(Recovered_cells < nrow(sObject@meta.data)) %>%
    slice(which.max(Recovered_cells)) %>%
    pull(Multiplet_rate) %>%
    as.numeric()
  
  if (length(multiplet_rate) == 0 || is.na(multiplet_rate)) {
    multiplet_rate <- 0.008
  }
  
  nExp.poi <- round(multiplet_rate * nrow(sObject@meta.data))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # DoubletFinder
  set.seed(42)
  sObject <- doubletFinder(
    seu  = sObject,
    PCs  = 1:min_pc,
    pK   = optimal.pk,
    nExp = nExp.poi.adj
  )
  
  df_col <- grep("^DF.classifications", colnames(sObject@meta.data), value = TRUE)
  colnames(sObject@meta.data)[colnames(sObject@meta.data) == df_col] <- "doublet_finder"
  
  N_doublets <- sum(sObject@meta.data$doublet_finder != "Singlet", na.rm = TRUE)
  NoDBcount <- sObject@assays$RNA@counts[, sObject@meta.data$doublet_finder=='Singlet']
  
  seurat_object_noDB <- CreateSeuratObject(
    counts = NoDBcount,
    project = sample_id
  )
  
  save(seurat_object_noDB, file = file.path(seurat_noDB_dir, 
                                            paste0(sample_id, "_SeuratObject_noDB.RData")))
  
  
  # after-DB QC summary
  qc_afterDB <- make_qc_summary(sObject@meta.data, sample_id = sample_id, N_doublets = N_doublets)
  all_afterDB_qc[[sample_id]] <- qc_afterDB
  all_afterDB_qc
  
  seurat_object_noDB <- add_qc_metrics(seurat_object_noDB, species = "human")
  seurat_obj_qc <- subset(
    seurat_object_noDB,
    subset = (nCount_RNA > 500) &
      (nCount_RNA < 50000) &
      (nFeature_RNA > 200) &
      (nFeature_RNA < 6000) &
      (percent.mt < 25)
  )
  
  seurat_obj_qc$sample_id <- sample_id
  seurat_obj_qc$condition <- sub("_.*$", "", sample_id)
  
  save(seurat_obj_qc, file = file.path(seurat_afterQC_dir, 
                                            paste0(sample_id, "_SeuratObject_afterQC.RData"))) 
  
}

# combine sample QC summary to 1 CSV file
initial_qc_all <- bind_rows(all_initial_qc)
afterDB_qc_all <- bind_rows(all_afterDB_qc)

# display
initial_qc_all_display <- format_percent_display(initial_qc_all)
afterDB_qc_all_display <- format_percent_display(afterDB_qc_all)

write.csv(initial_qc_all_display, file.path(summary_dir, "CRS_initial_QC_summary_ALL_display.csv"), row.names = FALSE)
write.csv(afterDB_qc_all_display, file.path(summary_dir, "CRS_afterDB_QC_summary_ALL_display.csv"), row.names = FALSE)

message("\nAll done! Summary tables saved to: ", summary_dir)

df <- afterDB_qc_all_display[, c(1,2,3,4,5,6,10,11,15)]
colnames(afterDB_qc_all_display)

df <- df %>%
  mutate(across(starts_with("pct"),
                ~ as.numeric(sub("%", "", .))))

tbl <- df %>%
  gt(rowname_col = "sample_id") %>%
  tab_header(
    title = "CRS Sample QC Summary"
  ) %>%
  
  # Rename columns (more readable)
  cols_label(
    N_total_cells = "Total Cells",
    N_doublets = "Doublets",
    pct_nUMI_le_500 = "<500 UMI",
    pct_nUMI_le_800 = "<800 UMI",
    pct_nGene_le_250 = "<250 Genes",
    pct_MT_gt_25 = "MT >25%",
    pct_Log10GenePerUMI_le_0_8 = "Complexity <0.8",
    pct_ribo_gt_25 = "Ribo >25%"
  ) %>%
  
  # Format numbers
  fmt_number(columns = c(N_total_cells, N_doublets), decimals = 0) %>%
  data_color(
    columns = starts_with("pct"),
    fn = function(x) {
      ifelse(
        x > 25, "red",
        ifelse(x > 15, "orange", "white")
      )
    }
  )%>%
  
  # Style
  cols_align(align = "center", -sample_id) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  )


tbl <- tbl %>%
  tab_footnote(
    footnote = "Complexity: log10(genes detected) / log10(total UMIs) per cell",
    locations = cells_column_labels(columns = pct_Log10GenePerUMI_le_0_8)
  )

tbl

# add mean_nUMI and mean_nGene
afterDB_qc_all_display <- read.csv("~/project/results_nofilter_gene/CRS_batch_outputs/QC_summary_tables/CRS_afterDB_QC_summary_ALL_display_update.csv")
colnames(afterDB_qc_all_display)
df <- afterDB_qc_all_display[, c(1,2,16,17,3,4,6,10,11,15)]

df <- df %>%
  mutate(across(starts_with("pct"),
                ~ as.numeric(sub("%", "", .))))

tbl <- df %>%
  gt(rowname_col = "sample_id") %>%
  tab_header(
    title = "CRS Sample QC Summary"
  ) %>%
  
  # Rename columns (more readable)
  cols_label(
    N_total_cells = "Total Cells",
    mean_nUMI = "mean_nUMI",
    mean_nGene = "mean_nGene",
    N_doublets = "Doublets",
    pct_nUMI_le_500 = "<500 UMI",
    pct_nGene_le_250 = "<250 Genes",
    pct_MT_gt_25 = "MT >25%",
    pct_Log10GenePerUMI_le_0_8 = "Complexity <0.8",
    pct_ribo_gt_25 = "Ribo >25%"
  ) %>%
  
  # Format numbers
  fmt_number(columns = c(N_total_cells, N_doublets), decimals = 0) %>%
  
  data_color(
    columns = starts_with("pct"),
    fn = function(x) {
      ifelse(
        x > 25, "red", "white"
      )
    }
  )%>%
  data_color(
    columns = mean_nGene,
    fn = function(x) {
      ifelse(
        x < 1000, "red", 
        ifelse(x < 1500, "orange", "white")
      )
    }
  )%>%
  data_color(
    columns = mean_nUMI,
    fn = function(x) {
      ifelse(
        x < 3000, "orange", "white"
      )
    }
  )%>%
  # Style
  cols_align(align = "center", -sample_id) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  )


tbl <- tbl %>%
  tab_footnote(
    footnote = "Complexity: log10(genes detected) / log10(total UMIs) per cell",
    locations = cells_column_labels(columns = pct_Log10GenePerUMI_le_0_8)
  )

tbl