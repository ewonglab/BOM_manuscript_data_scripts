# Script used for cross-species cell integration using SHAP scores
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
  library("gprofiler2", lib = "/g/data/zk16/software/Rpackages_paola/R_4.0.0")
  library("readxl")
  library(parallel)
  library(proxy)
  library(vegan)# This is vegan 2.5-7
  library("fda")
  # BiocManager::install("CCA", lib = "/g/data/zk16/software/Rpackages_paola/R_4.0.0")
  library("CCA")
})

# For UMAP colors
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

set.seed(123)

setwd("/g/data/zk16/cc3704/jack/bom/YSed")

# reading shap for mouse and human cells
shap.m <- readRDS("ysed_by_ysedMarkerP_matched_shap.rds")
shap.h <- readRDS("HumanHeart_by_ysedMarkerP_matched_shap.rds")

common_vars <- intersect(colnames(shap.m), colnames(shap.h))

shap.m$species <- "mouse"
shap.h$species <- "human"

# ---------------------------------- #
### SUBSET CELLS FOR EACH SPECIES

as.data.frame(table(shap.m$CellType))
#                Var1 Freq
# 1     ATAC_unknown1   10
# 2            B_cell   58
# 3                CM   13
# 4  Cardiac_Neuronal   36
# 5     Cardiomyocyte 4442
# 6       Coronary_EC 1306
# 7         Dendritic   23
# 8    Endocardial_EC  453
# 9        Fibroblast 3926
# 10             HSPC  128
# 11     Lymphatic_EC   96
# 12               M1  316
# 13               M2  309
# 14             M_Ex   87
# 15      Mesothelial   67
# 16               NK   40
# 17        Pericytes  721
# 18    Smooth_Muscle   66
# 19           T_cell   35

# remove Mex and CM from mouse (too few cells)
shap.m <- shap.m[!shap.m$CellType %in% c("CM", "M_Ex"), ]

as.data.frame(table(shap.h$CellType))
# Var1  Freq
# 1   AD   678
# 2   EC  7960
# 3   FB 23706
# 4   LC   664
# 5  MAC  6741
# 6  NER   754
# 7   SM  7298
# 8  aCM  8632
# 9  vCM 22296

### READ THE NUMBER OF ACCESSIBLE MARKER PEAKS PER CELL REMOVE THE CELLS WITH VERY LOW
# NUMBERS (LOW CONFIDENCE)
n_markers_mouse <- readRDS(file = "ysed_heart_n_markers_byCell.rds")
n_markers_human <- readRDS(file = "human_heart_n_markers_byCell.rds")

# Keeping cells with more than 20 marker peaks
shap.m <- shap.m[shap.m$CellID %in%
                   rownames(n_markers_mouse[n_markers_mouse$n_markers > 20,,drop = F]), ]
shap.h <- shap.h[shap.h$CellID %in%
                   rownames(n_markers_human[n_markers_human$n_markers > 20,,drop = F]), ]

# Function used to subsample cells from cell types with more than n_celltype cells
downsample_cells <- function(x, n_celltype, n_total){
  celltypes <- names(which(table(x$CellType) > n_celltype))
  # number of cells that should be sampled considerring cell types with no more than
  # n_celltype won't be sampled
  n_sample <- n_total - nrow(x[!x$CellType %in% celltypes, ])
  print(paste("A total of", n_sample, "cell will be sampled from cell types with more than"
              , n_celltype, "cells"))
  set.seed(123)
  cell_sample <- sample(x = x[x$CellType %in% celltypes, "CellID"]
                        , size = n_sample, replace = F)
  result <- x[(!x$CellType %in% celltypes) | (x$CellID %in% cell_sample), ]
  print("Number of cells in subsampled dataset")
  print(as.data.frame(table(result$CellType)))
  return(result)
}

# Subsampling for a total of 10,000 cells per species
shap.m_dnsamp <- downsample_cells(x = shap.m, n_celltype = 1000, n_total = 10000)
shap.h_dnsamp <- downsample_cells(x = shap.h, n_celltype = 1000, n_total = 10000)

# Removing motifs without variance
motifs <- common_vars[!common_vars %in% c("CellID", "CellType")]#

sd.h <- apply(shap.h_dnsamp[,motifs], 2, sd)
sd.m <- apply(shap.m_dnsamp[,motifs], 2, sd)

# table(sd.h == 0)
# FALSE  TRUE
# 6802  3005
# table(sd.m == 0)
# FALSE  TRUE
# 6802  3005

shap.h_dnsamp <- shap.h_dnsamp[, c("CellID", "CellType", "species", names(sd.h[sd.h != 0]))]
shap.m_dnsamp <- shap.m_dnsamp[, c("CellID", "CellType", "species", names(sd.m[sd.m != 0]))]

# Subset human and mouse data for motifs only
shap.m_motifs <- shap.m_dnsamp[, colnames(shap.m_dnsamp)[!colnames(shap.m_dnsamp) %in%
                                                           c("CellID", "CellType", "species")]]
shap.h_motifs <- shap.h_dnsamp[, colnames(shap.h_dnsamp)[!colnames(shap.h_dnsamp) %in%
                                                           c("CellID", "CellType", "species")]]

rownames(shap.h_motifs) <- shap.h_dnsamp$CellID
rownames(shap.m_motifs) <- shap.m_dnsamp$CellID

human.se <-  CreateSeuratObject(counts = t(shap.h_motifs), project = "human")
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
mouse.se <-  CreateSeuratObject(counts = t(shap.m_motifs), project = "mouse")
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

# Function to normalize each column by the sum of absolute values in that column
normalize_shap <- function(cell) {
  result <- cell / sum(abs(cell))
  return(result)
}

normalized_shap.h <- apply(human.se@assays$RNA@data, 2, normalize_shap)
normalized_shap.m <- apply(mouse.se@assays$RNA@data, 2, normalize_shap)

# Replace the original data in the Seurat object with the normalized values
human.se@assays$RNA@data <- normalized_shap.h
mouse.se@assays$RNA@data <- normalized_shap.m

# ---------------------------------- #
# HUMAN

#Find variable features
human.se <- FindVariableFeatures(human.se, selection.method = "vst", nfeatures = 2000)
all.motifs <- rownames(human.se)

# scale
human.se <- ScaleData(human.se, features = all.motifs)    

# Run PCA
human.se <- RunPCA(human.se, features = VariableFeatures(object = human.se))
# pdf("human_motifs_elbow.pdf")
# ElbowPlot(human.se, ndims = 50)
# dev.off()

# Run First UMAP
human.se <- FindNeighbors(human.se, dims = 1:40) #note dimensions

human.se <- FindClusters(human.se, resolution = 0.3)
human.se <- RunUMAP(human.se, dims = 1:30)

rownames(shap.h_dnsamp) <- shap.h_dnsamp$CellID
rownames(shap.m_dnsamp) <- shap.m_dnsamp$CellID

# same order of cells as in seurat objects
shap.h_dnsamp <- shap.h_dnsamp[colnames(human.se),]
shap.m_dnsamp <- shap.m_dnsamp[colnames(mouse.se),]

# sanity check
all(table(rownames(human.se@meta.data) == rownames(shap.h_dnsamp)))# T
all(table(rownames(mouse.se@meta.data) == rownames(shap.m_dnsamp)))# T

## ADD ANNOTATIONS TO METADATA - cell type and species
human.se@meta.data$CellType <- shap.h_dnsamp$CellType
human.se@meta.data$species <- shap.h_dnsamp$species

mouse.se@meta.data$CellType <- shap.m_dnsamp$CellType
mouse.se@meta.data$species <- shap.m_dnsamp$species

my_umap_cols <- alphabet(length(unique(human.se@meta.data$CellType)))
names(my_umap_cols) <- unique(human.se@meta.data$CellType)

pdf("human_normSHAP_10Kcells_SeuratUMAP_min20MarkerP.pdf")
DimPlot(human.se, reduction = "umap", group.by = "CellType", cols = my_umap_cols)
dev.off()

# ---------------------------------- #
# MOUSE

# Find variable features
mouse.se <- FindVariableFeatures(mouse.se, selection.method = "vst", nfeatures = 2000)
all.motifs <- rownames(mouse.se)
# scale
mouse.se <- ScaleData(mouse.se, features = all.motifs)    

#Run PCA
mouse.se <- RunPCA(mouse.se, features = VariableFeatures(object = mouse.se))

# pdf("mouse_motifs_elbow.pdf")
# ElbowPlot(mouse.se, ndims = 50)
# dev.off()


# Run First UMAP
mouse.se <- FindNeighbors(mouse.se, dims = 1:40) #note dimensions
mouse.se <- FindClusters(mouse.se, resolution = 0.3)
mouse.se <- RunUMAP(mouse.se, dims = 1:30)

my_umap_cols <- alphabet(length(unique(mouse.se@meta.data$CellType)))
names(my_umap_cols) <- unique(mouse.se@meta.data$CellType)

pdf("mouse_normSHAP_10Kcells_SeuratUMAP_min20MarkerP.pdf")
DimPlot(mouse.se, reduction = "umap", group.by = "CellType", cols = my_umap_cols)
dev.off()


# ---------------------------------- #
# INTEGRATE SPECIES

sp.anchors <- FindIntegrationAnchors(object.list = list(mouse.se, human.se), dims = 1:30)# default dims 1:30

sp.integrated <- IntegrateData(anchorset = sp.anchors, dims = 1:30)

library(ggplot2)
library(cowplot)
library(patchwork)

# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(sp.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
sp.integrated <- ScaleData(sp.integrated, verbose = FALSE)
sp.integrated <- RunPCA(sp.integrated, npcs = 30, verbose = FALSE)
sp.integrated <- RunUMAP(sp.integrated, reduction = "pca", dims = 1:30)

unique(sp.integrated@meta.data$CellType)
# [1] "Fibroblast"       "Cardiomyocyte"    "Pericytes"        "Coronary_EC"    
# [5] "M2"               "Endocardial_EC"   "T_cell"           "Lymphatic_EC"    
# [9] "Dendritic"        "Mesothelial"      "HSPC"             "B_cell"          
# [13] "M1"               "Smooth_Muscle"    "NK"               "Cardiac_Neuronal"
# [17] "ATAC_unknown1"    "FB"               "vCM"              "MAC"            
# [21] "SM"               "AD"               "EC"               "aCM"            
# [25] "NER"    

# match names for similar cell types
sp.integrated@meta.data$CellType_0 <- sp.integrated@meta.data$CellType
sp.integrated@meta.data$CellType <- sub("FB", "Fibroblast", sp.integrated@meta.data$CellType)
sp.integrated@meta.data$CellType <- sub("MAC", "Macrophage", sp.integrated@meta.data$CellType)
sp.integrated@meta.data$CellType <- sub("NER", "Nervous", sp.integrated@meta.data$CellType)
sp.integrated@meta.data$CellType <- sub("Cardiac_Neuronal", "Nervous", sp.integrated@meta.data$CellType)
sp.integrated@meta.data$CellType <- sub("SM", "Smooth_Muscle", sp.integrated@meta.data$CellType)

sp.integrated@meta.data$CellType <- sub("M1", "Macrophage"
                                        , sp.integrated@meta.data$CellType)
sp.integrated@meta.data$CellType <- sub("M2", "Macrophage"
                                        , sp.integrated@meta.data$CellType)

# change to cardiomyocytes:
sp.integrated@meta.data$CellType <- sub("vCM", "Cardiomyocyte"
                                        , sp.integrated@meta.data$CellType)
sp.integrated@meta.data$CellType <- sub("aCM", "Cardiomyocyte"
                                        , sp.integrated@meta.data$CellType)

# change to lymphocytes
sp.integrated@meta.data$CellType <- sub("LC", "Lymphocytes"
                                        , sp.integrated@meta.data$CellType)
sp.integrated@meta.data$CellType <- sub("NK", "Lymphocytes"
                                        , sp.integrated@meta.data$CellType)
sp.integrated@meta.data$CellType <- sub("B_cell", "Lymphocytes"
                                        , sp.integrated@meta.data$CellType)
sp.integrated@meta.data$CellType <- sub("T_cell", "Lymphocytes"
                                        , sp.integrated@meta.data$CellType)

# adjust adipocyte label
sp.integrated@meta.data$CellType <- sub("AD", "Adipocyte"
                                        , sp.integrated@meta.data$CellType)

# combine endothelial cells in mouse
sp.integrated@meta.data$CellType <- sub("Coronary_EC", "Endothelial"
                                        , sp.integrated@meta.data$CellType)
sp.integrated@meta.data$CellType <- sub("Endocardial_EC", "Endothelial"
                                        , sp.integrated@meta.data$CellType)
sp.integrated@meta.data$CellType <- sub("Lymphatic_EC", "Endothelial"
                                        , sp.integrated@meta.data$CellType)
sp.integrated@meta.data$CellType <- sub("EC", "Endothelial"
                                        , sp.integrated@meta.data$CellType)


# extract UMAP coordinates with metadata
umap_df <- sp.integrated@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% cbind(sp.integrated@meta.data[,c("CellType", "CellType_0", "species")])# 20000 4

my_umap_cols <- setNames(c("#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919"
                           , "#005C31", "#2BCE48", "#FFCC99", "#808080", "#94FFB5"
                           , "#8F7C00", "#9DCC00", "#C20088")
                         , c("Fibroblast", "Cardiomyocyte", "Pericytes", "Macrophage"
                             , "Endothelial", "Lymphocytes", "Dendritic", "Mesothelial"
                             , "HSPC", "Smooth_Muscle", "Nervous", "ATAC_unknown1", "Adipocyte"))

pdf("integrated_normSHAP_min20MarkP_SeuratUMAP.pdf")
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = CellType)) + geom_point(size = 1) +
  theme_classic() + theme(legend.position = "bottom") +
  scale_color_manual(values = my_umap_cols)
dev.off()

my_umap_cols <- alphabet(length(unique(umap_df$CellType_0)))
names(my_umap_cols) <- unique(umap_df$CellType_0)

pdf("integrated_normSHAP_min20MarkP_SeuratUMAP.1.pdf")
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = CellType_0)) + geom_point(size = 1) +
  theme_classic() + theme(legend.position = "bottom") +
  scale_color_manual(values = my_umap_cols)
dev.off()

pdf("integrated_normSHAP_min20MarkP_SeuratUMAP_bySp.pdf")
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = species)) + geom_point(size = 1) +
  theme_classic() + theme(legend.position = "bottom") #+
# scale_color_manual(values = my_umap_cols)
dev.off()

saveRDS(sp.integrated, "integrated_10Kcells_min20MarkP_normSHAP_Seu.rds")
