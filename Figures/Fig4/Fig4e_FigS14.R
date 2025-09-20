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
  library(proxy)
  library(vegan)# This is vegan 2.5-7
  library("fda")
  library("CCA")
  library("ggpubr")
  library(ggbeeswarm)
})

alphabet <- function(n = 26) {
  if(n > 26){
    message("Only 26 colors are available with 'alphabet'")
    n <- 26
  }
  pal <- c("#F0A0FF","#0075DC","#993F00","#4C005C","#191919","#005C31",
           "#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00",
           "#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010",
           "#5EF1F2","#00998F","#E0FF66","#740AFF","#990000","#FFFF80",
           "#FFE100","#FF5005")
  names(pal) <- c("amethyst","blue","caramel","damson","ebony","forest",
                  "green","honeydew","iron","jade","khaki","lime","magenta",
                  "navy","orange","pink","quagmire","red","sky","turquoise",
                  "uranium","violet","wine","xanthin","yellow","zinnia")
 
  return(pal[1:n])
}

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

normalize_shap <- function(cell) {
  result <- cell / sum(abs(cell))  # Divide by sum of absolute values
  return(result)
}

setwd("human_by_ysed/ysed_padj0.05")

# aggregated shap files
CellTypes <- c("aCM", "MAC", "EC", "NER", "FB", "SM", "LC", "vCM", "AD")

shap_files <- paste0("HumanHeart_by_ysed.posMarkerP_padj0.05_matched_", CellTypes, "cells_meanSHAP.txt")

shap <- lapply(shap_files, fread)

sapply(shap, ncol)
# [1] 9 9 9 9 9 9 9 9 9
sapply(shap, nrow)
# [1] 12127960  9471105 11183800  1059370 33306930 10253690   932920 31325880
# [9]   952590
sapply(shap, nrow)/1405 # 1405 number of motifs in this model
# [1]  8632  6741  7960   754 23706  7298   664 22296   678

### normalising SHAP by cell to test whether cell prediction improves
norm_shap_all_cells <- function(x, nclasses){
  x <- as.data.frame(x)
  # define names of columns containing columns
  model_classes <- paste0("class_", 1:nclasses)
  # separate shap values into a list where every list contains the values of a model class
  shap_by_class <- lapply(model_classes, function(model_class)
    x[,c(model_class, "CellID", "motif")])
  # cells as columns
  shap_by_class <- lapply(shap_by_class, function(model_class)
    tidyr::spread(model_class, CellID, colnames(model_class)[1]))
  # shap as row names
  for(i in 1:nclasses){
    rownames(shap_by_class[[i]]) <- shap_by_class[[i]]$motif
    shap_by_class[[i]]$motif <- NULL
  }
  # normalize shap by cell
  shap_by_class <- lapply(shap_by_class, function(model_class)
    apply(model_class, 2, normalize_shap))
  shap_by_class <- lapply(shap_by_class, as.data.frame)
  ## motif name as a column again
  for(i in 1:nclasses){
    shap_by_class[[i]]$motif <- rownames(shap_by_class[[i]])
  }
  # cells as rows again
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

# get normalised shap - normalised by cell
shap <- lapply(shap, function(celltype) norm_shap_all_cells(x = celltype, nclasses = 7))

# add a column with the cell type of every cells set
for(i in 1:length(shap)){
  shap[[i]]$CellType <- CellTypes[i]
}

shap <- do.call("rbind", shap)#
dim(shap)
# [1] 110614245        10
dim(shap[complete.cases(shap),]) #
# [1] 110614245        10

model_classes <- c("Cardiomyocyte", "Endothelial", "Fibroblast"
                   , "Lymphocyte", "Macrophage", "Nervous", "Smooth_Muscle")
colnames(shap)[1:7] <-  model_classes

# calculate mean shap across cells and model classes
shap$motif <- NULL

shap_sum.pos <- aggregate(. ~ CellID + CellType, shap
                          , function(x) sum(x[x > 0], na.rm = TRUE))#78729     9
shap_sum.abs <- aggregate(. ~ CellID + CellType, shap
                          , function(x) sum(abs(x), na.rm = TRUE))# 78729     9

shap_sum.pos <- melt(shap_sum.pos, varnames = c("CellID", "CellType"))
shap_sum.abs <- melt(shap_sum.abs, varnames = c("CellID", "CellType"))

colnames(shap_sum.pos)[3:4] <- c("ModelClass", "sum.pos.shap")
colnames(shap_sum.abs)[3:4] <- c("ModelClass", "sum.abs.shap")

norm_shap <- merge(shap_sum.pos, shap_sum.abs, by = c("CellID", "CellType", "ModelClass"))# 551103      5
norm_shap$norm_shap <- with(norm_shap, sum.pos.shap/sum.abs.shap) # 551103      6

norm_shap <- tidyr::spread(norm_shap[,c("CellID", "CellType", "ModelClass", "norm_shap")]
                           , ModelClass, norm_shap)# 78729     9

# zscore of normalised shap - zscore per model class
norm_shap.z <- cbind(norm_shap[,c("CellID", "CellType")]
                     , scale(x = norm_shap[,3:ncol(norm_shap)], center = T, scale = T))

modelClasses <- c("Cardiomyocyte", "Endothelial", "Fibroblast"
                  , "Lymphocyte", "Macrophage", "Nervous", "Smooth_Muscle")

norm_shap.z$maxclass <- apply(norm_shap.z[,modelClasses], 1, function(x) modelClasses[which.max(x)])

# saving normalized and aggregated SHAP values

saveRDS(object = norm_shap.z, file = "human_by_ysedMarkerP_padj0.05_matched_byCellNorm_z.norm.shap.rds")

# norm_shap.z <- readRDS("human_by_ysedMarkerP_matched_z.norm.shap.rds")# 78729    10

# reading the number of marker peaks per cell
n_markers <- readRDS(file = "human_heart_n_markers_byCell.rds")# 78729     1
# keep cells with more than 20 marker peaks
n_markers <- n_markers[n_markers$n_markers > 20,,drop = F] # 55886     1

norm_shap.z <- norm_shap.z[norm_shap.z$CellID %in% rownames(n_markers), ]# 55886    10

# confusion matrix
conf.mat <- as.data.frame(table(norm_shap.z[,c("CellType", "maxclass")]))
colnames(conf.mat)[3] <- "Freq_predicted"

n.actual <- as.data.frame(table(norm_shap.z[,c("CellType")]))
colnames(n.actual) <- c("CellType", "n")

conf.mat <- merge(conf.mat, n.actual, by = "CellType", all.x = T)
conf.mat$pred_ratio <- with(conf.mat, Freq_predicted/n) # 56  5

# conf.mat <- tidyr::spread(conf.mat[,c("CellType", "maxclass", "pred_ratio")]
#                            , maxclass, pred_ratio)#

pdf("human_min20MarkP_by_ysedMarker_padj0.05_matched_bycellnormSHAP_confMat.pdf")
ggplot(conf.mat, aes(x = maxclass, y = CellType, fill = pred_ratio)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels 90 degrees
dev.off()

## prediction stats
norm_shap.z$CellType.1 <- norm_shap.z$CellType
norm_shap.z$CellType.1 <- sub("EC", "Endothelial", norm_shap.z$CellType.1)
norm_shap.z$CellType.1 <- sub("FB", "Fibroblast", norm_shap.z$CellType.1)
norm_shap.z$CellType.1 <- sub("SM", "Smooth_Muscle", norm_shap.z$CellType.1)
norm_shap.z$CellType.1 <- sub("aCM", "Cardiomyocyte", norm_shap.z$CellType.1)
norm_shap.z$CellType.1 <- sub("vCM", "Cardiomyocyte", norm_shap.z$CellType.1)
norm_shap.z$CellType.1 <- sub("LC", "Lymphocyte", norm_shap.z$CellType.1)
norm_shap.z$CellType.1 <- sub("MAC", "Macrophage", norm_shap.z$CellType.1)
norm_shap.z$CellType.1 <- sub("NER", "Nervous", norm_shap.z$CellType.1)

unique(norm_shap.z$CellType.1)
# [1] "Cardiomyocyte" "Macrophage"    "Endothelial"   "Nervous"      
# [5] "Fibroblast"    "Smooth_Muscle" "AD"  


human_stats <- calculate_metrics(norm_shap.z)
write.csv(x = human_stats, file = "human_min20MarkP_byCellNorm_z.norm.shap_STATS.csv"
          , quote = F)

human_stats <- melt(human_stats, varnames = "CellType")

my_theme <-  theme(panel.border = element_blank(), panel.grid.major = element_blank()
                   , panel.grid.minor = element_blank()#legend.position="none"
                   , axis.line = element_line(colour = "black", linewidth = 1)
                   , legend.title = element_blank()
                   , legend.key=element_blank()
                   , legend.position = "bottom"
                   , legend.key.width = unit(1.5,"cm")
                   , panel.background = element_blank()
                   , text = element_text(size=12)
                   , legend.text=element_text(size=8)
                   , axis.text.x=element_text(colour="black")
                   , axis.text.y=element_text(colour="black")
                   , axis.ticks = element_line(colour = "black"))

human_stats$CellType <- sub("AD", "Adipocyte", human_stats$CellType)
my_colors <- setNames(c("#F0A0FF", "#0075DC", "#4C005C", "#191919"
                        , "#005C31", "#94FFB5", "#8F7C00", "#C20088")
                      , c("Fibroblast", "Cardiomyocyte", "Macrophage"
                          , "Endothelial", "Lymphocyte", "Smooth_Muscle"
                          , "Nervous", "Adipocyte"))

# remove adipocyte
human_stats <- human_stats[human_stats$CellType != "Adipocyte",]

pdf("human_min20MarkP_by_ysedMarkerP_padj0.05_matched_byCellNorm_z.norm.shap_STATS_bw.pdf")
ggplot(human_stats, mapping = aes(variable, value, color = CellType)) +
  geom_beeswarm(dodge.width=.8,cex=2) + my_theme +
  stat_summary(
    fun = "mean", geom = "point", aes(group = variable),
    color = "black", size = 10, shape=95) +  scale_color_manual(values = my_colors) +
  scale_x_discrete(expand = c(0.1, 0.1))
dev.off()
