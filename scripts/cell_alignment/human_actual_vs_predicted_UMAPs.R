# Script used to plot human heart cells on the cross-species integrated UMAP space. The UMAPs show the human cells coloured by 1) actuall cell type and 2) predicted cell types based on SHAP
suppressMessages({
  library(reshape2)
  library(umap)
  library(spatstat.core)
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
})


# with normalization by cell
norm_shap.z <- readRDS("human_by_ysedMarkerP_matched_byCellNorm_z.norm.shap.rds")

# Update cell type names
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
# [1] "Cardiomyocyte" "Macrophage"    "Endothelial"   "Nervous"      
# [5] "Fibroblast"    "Smooth_Muscle" "Lymphocyte"    "AD"          

## READ INTEGRATED UMAP (SHAP normalised by cell)
sp.integrated <- readRDS("integrated_10Kcells_min20MarkP_normSHAP_Seu.rds")

# unique(sp.integrated@meta.data$CellType)
# [1] "Fibroblast"    "Cardiomyocyte" "Pericytes"     "Endothelial"  
# [5] "Macrophage"    "Lymphocytes"   "Dendritic"     "Mesothelial"  
# [9] "HSPC"          "Smooth_Muscle" "Nervous"       "ATAC_unknown1"
# [13] "Adipocyte"    

# Extract UMAP coordinates with metadata
umap_df <- sp.integrated@reductions$umap@cell.embeddings %>%
as.data.frame() %>% cbind(sp.integrated@meta.data[,c("CellType", "CellType_0", "species")])#

## Subset to human cells used in cross-species integration
umap_df <- umap_df[umap_df$species == "human",]
umap_df <- merge(umap_df, norm_shap.z[, c("CellID", "maxclass")]
                 , by.x = 0, by.y = "CellID", all.x = T)

# save colours as in previous UMAPs
my_umap_cols <- setNames(c("#F0A0FF", "#0075DC", "#4C005C", "#191919"
                           , "#005C31", "#94FFB5", "#8F7C00", "#C20088")
                         , c("Fibroblast", "Cardiomyocyte", "Macrophage"
                             , "Endothelial", "Lymphocyte", "Smooth_Muscle"
                             , "Nervous", "Adipocyte"))

# as.data.frame(table(umap_df$CellType))
#            Var1 Freq
# 1     Adipocyte   52
# 2 Cardiomyocyte 3923
# 3   Endothelial  619
# 4    Fibroblast 3541
# 5    Macrophage  882
# 6       Nervous  121
# 7 Smooth_Muscle  862

# Actualcell types
# unique(umap_df$CellType)
# [1] "Fibroblast"    "Cardiomyocyte" "Macrophage"    "Smooth_Muscle"
# [5] "Adipocyte"     "Endothelial"   "Nervous"      
# Predicted cell types based on SHAP
# unique(umap_df$maxclass)
# [1] "Lymphocyte"    "Cardiomyocyte" "Fibroblast"    "Endothelial"  
# [5] "Nervous"       "Macrophage"    "Smooth_Muscle"

pdf("integrated_byCellNorm_shap_SeuratUMAP_min20MarkerP_by_z.norm.shap.pdf")
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = maxclass)) +
  geom_point(size = 1) +
  theme_classic() + theme(legend.position = "bottom") +
  scale_color_manual(values = my_umap_cols)
dev.off()

pdf("integrated_byCellNorm_shap_SeuratUMAP_min20MarkerP_by_groundtruth.pdf")
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = CellType)) + geom_point(size = 1) +
  theme_classic() + theme(legend.position = "bottom") +
  scale_color_manual(values = my_umap_cols)
dev.off()
