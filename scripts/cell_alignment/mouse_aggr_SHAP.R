# This script is used to aggregate SHAP scrores across marker peaks for each cell of a mouse cell type
# The script will remove peaks that were used for model training and validation before aggregating SHAP

suppressMessages({
  library("Seurat")
  library("Signac")
  library("GenomeInfoDb")
  library("matrixStats")
  library("DESeq2")
  library("Matrix")
  library("Matrix.utils")
  library(GenomicFeatures)
  library(GenomicRanges)
  library(ggplot2)
  library(pheatmap)
  library(data.table)
  library(tidyr)
  library(dplyr)
  library("gprofiler2")
  library("readxl")
  library(parallel)
})

## Function to aggregate SHAP for the accessible marker peaks of each cell
parallel_meanSHAP <- function(cell_ids, shap_values, OUTFILE, num_cores) {
  # Run the function in parallel for each cell_id
  result_list <- mclapply(cell_ids, function(cell_id) {
    # Marker peaks in cell
    CELL_peaks <- peak_ids[cell_marker_indices[[cell_id]], "trimmed_id"]
    
    # Compute mean SHAP values
    CELL_shap_mean <- sapply(names(shap_values), function(model_class) {
      shap_matrix <- shap_values[[model_class]]
      mean_shap <- colMeans(shap_matrix[CELL_peaks, , drop = FALSE])
    }, simplify = FALSE)
    
    # Add CellID to the results
    CELL_shap_mean <- as.data.frame(CELL_shap_mean)
    CELL_shap_mean$CellID <- cell_id
    CELL_shap_mean$motif <- rownames(CELL_shap_mean)
    rownames(CELL_shap_mean) <- 1:nrow(CELL_shap_mean)
    
    return(CELL_shap_mean)
  }, mc.cores = num_cores)
  
  # Combine all results into a single data frame
  final_result <- do.call(rbind, result_list)
  write.table(final_result, file = OUTFILE, sep = ",", row.names = FALSE,
              col.names = !file.exists(OUTFILE), append = file.exists(OUTFILE), quote = FALSE)
}

# The script will be executed for the cells of the cell type provided as argument
args <- commandArgs(trailingOnly=TRUE)

CellType <- args[1]
my_number_cores <- as.integer(args[2])


# Motif counts for union peaks already filtered to keep only peaks overlapping marker peaks
counts <- read.table(file = "mouse_union_ov_markerAllct_ysedMarker.matched_counts.txt"
                     , sep = '\t', header = T, stringsAsFactors = F)

# Reading SHAP scores for each model class
shap_files <- paste0("union_allMarkers_by_ysed.markerP_matched_SHAP_", 1:7, ".txt")
shap <- lapply(shap_files, fread, header = F)
shap <- lapply(shap, as.data.frame)

# Add column and row names to SHAP
for(i in 1:length(shap)){
  colnames(shap[[i]]) <- colnames(counts)
  rownames(shap[[i]]) <- sub(":", "-", rownames(counts))
}


# Read Seurat object with ATAC counts and subset to young_sed condition
ysed_file <- "2024-03-01_Full_ATAC_MatchAnnotated.rds"
ysed.se <- readRDS(ysed_file)
ysed.se <- subset(ysed.se, subset = age_condition == "Young_Sed")
DefaultAssay(ysed.se) <- "ATAC" 

### reading peak coordinates (same as atac matrix rows) and trimmed peaks coordinates
peaks <- read.table(file = "YoungSed_min1rawcount.bed"
                    , sep = '\t', header = F, stringsAsFactors = F)
trimmed <- read.table(file = "YoungSed_min1rawcount_500bp.bed"
                      , sep = '\t', header = F, stringsAsFactors = F)
# peaks and trimmed peaks have the same order
peaks$trimmed_id <- trimmed$V4
peaks$V4 <- sub(":", "-", peaks$V4)
peaks$trimmed_id <- sub(":", "-", peaks$trimmed_id)
peaks$trimmed_id <- sub("chr", "", peaks$trimmed_id)

# Keep peak coordinates for those overlapping marker peaks
peaks <- peaks[peaks$trimmed_id %in% rownames(shap[[1]]),]

# -------------------------------------------------------- #
### Remove peaks used for model training

# Read motif counts including training, validation and test peaks (marker peaks of matched cell types)
markerP_motifs <- read.table(file = "ysed_markerP_matched_500bp_gimme_q0.5.txt", header = T
                             , stringsAsFactors = F, sep = '\t')
# Read peaks in test set
test_peaks <- read.table(file = "ysed_markerP_matched_500bp_gimme_q0.5_6915T_pred.txt"
                         , header = T, stringsAsFactors = F) 
# Keep training and validation peak ids
train_val.peaks <- setdiff(rownames(markerP_motifs), rownames(test_peaks))#
train_val.peaks <- as.data.frame(train_val.peaks)
colnames(train_val.peaks) <- "peakid"
train_val.peaks$chr <- sub(":.*", "", train_val.peaks$peakid)
train_val.peaks$start <- sub("-.*", "", sub(".*:", "", train_val.peaks$peakid))
train_val.peaks$end <- sub(".*-", "", sub(".*:", "", train_val.peaks$peakid))
train_val.peaks$start <- as.integer(train_val.peaks$start)
train_val.peaks$end <- as.integer(train_val.peaks$end)

# Remove peaks that overlap training and validation peaks
train_val.peaks_gr <- with(train_val.peaks, GRanges(chr, IRanges(start + 1, end)))
peaks_gr <- with(peaks, GRanges(V1, IRanges(V2, V3)))# 1-based
x <- as.data.frame(findOverlaps(peaks_gr, train_val.peaks_gr))
peaks <- peaks[-unique(x$queryHits),]

# -------------------------------------------------------- #

# reorder peaks in shap matrices keeping only the filtered peaks
for(i in 1:length(shap)){
  shap[[i]] <- shap[[i]][peaks$trimmed_id,]
}

peak_ids <- peaks


# -------------------------------------------------------- #
# For each cell, keep accessible marker peaks in the same cluster
# Read file that indicates the overlap of peaks with marker peaks and cluster
peaks_ov_markers <- read.table(file = "union_ov_ysed_markerPeaks", header = T
                               , stringsAsFactors = F, sep = '\t')

# -------------------------------------------------------- #
# Subset atac matrix to keep filtered peaks
ysed.se <- ysed.se[peak_ids$V4,] 

# -------------------------------------------------------- #
# Get peaks with non-zero atac count for each cell
row_indices <- ysed.se@assays$ATAC@counts@i + 1  # Row indices (zero-based, so add 1 for R indexing)
column_pointers <- ysed.se@assays$ATAC@counts@p  # Column pointers

# Initialize list to store non-zero indices per column
non_zero_indices <- vector("list", ncol(ysed.se@assays$ATAC@counts))

# Loop through columns using column pointers to get non-zero indices
for(j in 1:ncol(ysed.se@assays$ATAC@counts)) {
  start <- column_pointers[j] + 1
  end <- column_pointers[j + 1]
  non_zero_indices[[j]] <- row_indices[start:end]
}

names(non_zero_indices) <- colnames(ysed.se@assays$ATAC@counts)

# -------------------------------------------------------- #
# Subset to marker peaks in the same cluster

# add an index to peaks in peaks_ov_markers. The index reflects their position in the atac matrix
peaks_ov_markers$idx <- sapply(peaks_ov_markers$peakid
                               , function(x) which(rownames(ysed.se@assays$ATAC@counts) == x))

# Function to keep peaks that are marker peaks in the same cluster as the query cell
keep_marker_sameCellType <- function(cell, allCellsPeaks){
  cellpeaks <- allCellsPeaks[[cell]]
  celltype <- ysed.se@meta.data[ysed.se@meta.data$cellID == cell, "CellType"]
  result <- cellpeaks[cellpeaks %in% peaks_ov_markers[peaks_ov_markers$cluster == celltype, "idx"]]
  return(result)
}

cell_marker_indices <- lapply(names(non_zero_indices)
                              , function(x) keep_marker_sameCellType(cell = x
                                                                     , allCellsPeaks = non_zero_indices))
names(cell_marker_indices) <- names(non_zero_indices)
cell_marker_indices <- cell_marker_indices[lapply(cell_marker_indices, length) > 0]

# model classes as names for shap list
names(shap) <- c("Cardiomyocyte", "Endothelial", "Fibroblast"
                 , "Lymphocyte", "Macrophage", "Nervous", "Smooth_Muscle")

# extract cell ids for query cell type
cells <- ysed.se@meta.data[ysed.se@meta.data$CellType == CellType, "cellID"]
cells <- cells[cells %in% names(cell_marker_indices)]

# aggregate shap for cell set
parallel_meanSHAP(cell_ids = cells
                  , shap_values = shap
                  , OUTFILE = paste0("m.union_allMark_by_ysMarkerP_matched_", CellType, "cells_mSHAP_ho.txt")
                  , num_cores = my_number_cores)
