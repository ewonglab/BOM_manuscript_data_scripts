
suppressMessages({
  library(reshape2)
  library(umap)
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
  library(proxy)
  library(vegan)# This is vegan 2.5-7
  library("fda")
  library("CCA")
  library("ggpubr", lib = "/g/data/zk16/qing/tools/R/4.0")
  library(reshape2)
  library(ggbeeswarm)
})

# Function used to calculate cell type prediction statistics
calculate_metrics <- function(df) {
  # Get unique cell types
  cell_types <- unique(df$CellType.1)
 
  # Initialize a list to store results for each cell type
  metrics_list <- list()
 
  for (cell_type in cell_types) {
    # True Positives, False Positives, False Negatives, True Negatives for the current cell type
    tp <- sum(df$maxclass == cell_type & df$CellType.1 == cell_type)
    fp <- sum(df$maxclass == cell_type & df$CellType.1 != cell_type)
    fn <- sum(df$maxclass != cell_type & df$CellType.1 == cell_type)
    tn <- sum(df$maxclass != cell_type & df$CellType.1 != cell_type)
   
    # Accuracy
    accuracy <- (tp + tn) / (tp + tn + fp + fn)
   
    # Precision (avoid division by zero)
    precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
   
    # Recall (Sensitivity)
    recall <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
   
    # F1 Score
    f1 <- ifelse(!is.na(precision) & !is.na(recall) & (precision + recall) > 0,
                 2 * (precision * recall) / (precision + recall), NA)
   
    tp <- as.numeric(tp)
    tn <- as.numeric(tn)
    fp <- as.numeric(fp)
    fn <- as.numeric(fn)
   
    mcc_numerator <- (tp * tn) - (fp * fn)
    mcc_denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc <- ifelse(mcc_denominator > 0, mcc_numerator / mcc_denominator, NA)
   
    # Store metrics for the current cell type
    metrics_list[[cell_type]] <- data.frame(
      CellType = cell_type,
      Accuracy = accuracy,
      Precision = precision,
      Recall = recall,
      F1_Score = f1,
      MCC = mcc
    )
  }
 
  # Combine all metrics into one data frame
  metrics_df <- bind_rows(metrics_list)
 
  return(metrics_df)
}

# Function used to normalize SHAP per cell
normalize_shap <- function(cell) {
  result <- cell / sum(abs(cell))  # Divide by sum of absolute values
  return(result)
}

### Function used to normalize SHAP scores for all cells
norm_shap_all_cells <- function(x, nclasses){
  x <- as.data.frame(x)

  # define names of columns containing SHAP scores
  model_classes <- paste0("class_", 1:nclasses)

  # separate shap values into a list where every list contains the values of a model class
  shap_by_class <- lapply(model_classes, function(model_class)
    x[,c(model_class, "CellID", "motif")])

  # cells as columns
  shap_by_class <- lapply(shap_by_class, function(model_class)
    tidyr::spread(model_class, CellID, colnames(model_class)[1]))

  # Motifs as row names
  for(i in 1:nclasses){
    rownames(shap_by_class[[i]]) <- shap_by_class[[i]]$motif
    shap_by_class[[i]]$motif <- NULL
  }

  # normalize shap by cell
  shap_by_class <- lapply(shap_by_class, function(model_class)
    apply(model_class, 2, normalize_shap))
  shap_by_class <- lapply(shap_by_class, as.data.frame)

  ## motif name as a column
  for(i in 1:nclasses){
    shap_by_class[[i]]$motif <- rownames(shap_by_class[[i]])
  }

  # cells as rows
  shap_by_class <- lapply(shap_by_class, function(model_class)
    reshape2::melt(model_class, varnames = "motif"))

  # add model class as column
  for(i in 1:nclasses){
    shap_by_class[[i]]$model_class <- model_classes[i]
  }

  # rbind across model classes
  shap_by_class <- do.call("rbind", shap_by_class)

  # model classes as columns:
  shap_by_class <- tidyr::spread(shap_by_class, model_class, value)
  colnames(shap_by_class) <- sub("variable", "CellID", colnames(shap_by_class))
  setcolorder(shap_by_class, c(model_classes, "CellID", "motif"))
  return(shap_by_class)
}


# Human cell types
CellTypes <- c("aCM", "MAC", "EC", "NER", "FB", "SM", "LC", "vCM", "AD")

# Reading SHAP values for the cells of each cell type
shap_files <- paste0("HumanHeart_by_ysedMarkerP_matched_", CellTypes, "cells_meanSHAP.txt")
shap <- lapply(shap_files, fread)

# get normalised shap - normalised by cell
shap <- lapply(shap, function(celltype) norm_shap_all_cells(x = celltype, nclasses = 7))

# add a column with the cell type of every cells set
for(i in 1:length(shap)){
  shap[[i]]$CellType <- CellTypes[i]
}

shap <- do.call("rbind", shap)

model_classes <- c("Cardiomyocyte", "Endothelial", "Fibroblast"
                   , "Lymphocyte", "Macrophage", "Nervous", "Smooth_Muscle")
colnames(shap)[1:7] <-  model_classes

shap$motif <- NULL

# Normalize sum of positive SHAP by the sum of absolute SHAP
shap_sum.pos <- aggregate(. ~ CellID + CellType, shap
                          , function(x) sum(x[x > 0], na.rm = TRUE))#78729     9
shap_sum.abs <- aggregate(. ~ CellID + CellType, shap
                          , function(x) sum(abs(x), na.rm = TRUE))# 78729     9

shap_sum.pos <- melt(shap_sum.pos, varnames = c("CellID", "CellType"))
shap_sum.abs <- melt(shap_sum.abs, varnames = c("CellID", "CellType"))

colnames(shap_sum.pos)[3:4] <- c("ModelClass", "sum.pos.shap")
colnames(shap_sum.abs)[3:4] <- c("ModelClass", "sum.abs.shap")

norm_shap <- merge(shap_sum.pos, shap_sum.abs, by = c("CellID", "CellType", "ModelClass"))# 551103      5
norm_shap$norm_shap <- with(norm_shap, sum.pos.shap/sum.abs.shap)

norm_shap <- tidyr::spread(norm_shap[,c("CellID", "CellType", "ModelClass", "norm_shap")]
                           , ModelClass, norm_shap)# 78729     9

# zscore per model class of normalized SHAP
norm_shap.z <- cbind(norm_shap[,c("CellID", "CellType")]
                     , scale(x = norm_shap[,3:ncol(norm_shap)], center = T, scale = T))

modelClasses <- c("Cardiomyocyte", "Endothelial", "Fibroblast"
                  , "Lymphocyte", "Macrophage", "Nervous", "Smooth_Muscle")

# Identify model class with highers score for each cell (predicted labels)
norm_shap.z$maxclass <- apply(norm_shap.z[,modelClasses], 1, function(x) modelClasses[which.max(x)])

# saving normalized SHAP with predicted labels
saveRDS(object = norm_shap.z, file = "human_by_ysedMarkerP_matched_byCellNorm_z.norm.shap.rds")

# Reading the number of marker peaks per cell
n_markers <- readRDS(file = "human_heart_n_markers_byCell.rds")
# keep cells with more than 20 marker peaks
n_markers <- n_markers[n_markers$n_markers > 20,,drop = F]

norm_shap.z <- norm_shap.z[norm_shap.z$CellID %in% rownames(n_markers), ]

# Confusion matrix: frequency of actual and predicted cell types
conf.mat <- as.data.frame(table(norm_shap.z[,c("CellType", "maxclass")]))
colnames(conf.mat)[3] <- "Freq_predicted"

n.actual <- as.data.frame(table(norm_shap.z[,c("CellType")]))
colnames(n.actual) <- c("CellType", "n")

# Calculate ratio between the predicted cell types and sthe total number of cells by cell type
conf.mat <- merge(conf.mat, n.actual, by = "CellType", all.x = T)
conf.mat$pred_ratio <- with(conf.mat, Freq_predicted/n) 

# Save heatmap
pdf("human_min20MarkP_by_ysedMarker_matched_bycellnormSHAP_confMat.pdf")
ggplot(conf.mat, aes(x = maxclass, y = CellType, fill = pred_ratio)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels 90 degrees
dev.off()

## Adjust cell type names
norm_shap.z$CellType.1 <- norm_shap.z$CellType
norm_shap.z$CellType.1 <- sub("EC", "Endothelial", norm_shap.z$CellType.1)
norm_shap.z$CellType.1 <- sub("FB", "Fibroblast", norm_shap.z$CellType.1)
norm_shap.z$CellType.1 <- sub("SM", "Smooth_Muscle", norm_shap.z$CellType.1)
norm_shap.z$CellType.1 <- sub("aCM", "Cardiomyocyte", norm_shap.z$CellType.1)
norm_shap.z$CellType.1 <- sub("vCM", "Cardiomyocyte", norm_shap.z$CellType.1)
norm_shap.z$CellType.1 <- sub("LC", "Lymphocyte", norm_shap.z$CellType.1)
norm_shap.z$CellType.1 <- sub("MAC", "Macrophage", norm_shap.z$CellType.1)
norm_shap.z$CellType.1 <- sub("NER", "Nervous", norm_shap.z$CellType.1)

# unique(norm_shap.z$CellType.1)
# [1] "Cardiomyocyte" "Fibroblast"    "Macrophage"    "Endothelial"  
# [5] "Smooth_Muscle" "AD"            "Nervous"

# Calculate cell type prediction performance metrics
human_stats <- calculate_metrics(norm_shap.z)
#       CellType  Accuracy  Precision    Recall   F1_Score       MCC
# 1 Cardiomyocyte 0.9456214 0.99352185 0.8690021 0.92709958 0.8890874
# 2    Macrophage 0.9639087 0.86203373 0.6955267 0.76988021 0.7556327
# 3   Endothelial 0.9369788 0.50012300 0.5770650 0.53584607 0.5037055
# 4       Nervous 0.9689368 0.05013928 0.7438017 0.09394572 0.1880211
# 5    Fibroblast 0.9361736 0.95481474 0.8666471 0.90859705 0.8621220
# 6 Smooth_Muscle 0.9398418 0.66868262 0.5474704 0.60203598 0.5731828
# 7            AD 0.9990695         NA 0.0000000         NA        NA

# Save satatistics
# Save prediction metrics beeswarm plot
