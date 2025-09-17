library(ggplot2)
library("ggpubr")
library("mclust")
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(reshape2)
library("data.table")

training_motif_mean <- function(x){
  train_1 <- x[x$binary_celltype==1,]
  train_0 <- x[x$binary_celltype==0,]
  train_1$binary_celltype <- NULL
  train_0$binary_celltype <- NULL
  #column sums
  train_1_colSums <- as.data.frame(apply(train_1, 2, mean))
  train_0_colSums <- as.data.frame(apply(train_0, 2, mean))
  train_colSums <- merge(train_1_colSums, train_0_colSums, by=0)
  colnames(train_colSums) <- c("motif", "sum_train_1", "sum_train_0")
  train_colSums$direction <- ifelse(train_colSums$sum_train_1 > train_colSums$sum_train_0
                                    , "enriched_in_class_1", "enriched_in_class_0")
  rownames(train_colSums) <- train_colSums$motif
  return(train_colSums)  
}

top_shap <- function(shap.df, train.df){
  abs.sum.shap_per_motif <- apply(shap.df, 2, function(x) sum(abs(x)))
  names(abs.sum.shap_per_motif) <- colnames(train.df)[-ncol(train.df)]
  abs.sum.shap_per_motif <- abs.sum.shap_per_motif[order(-abs.sum.shap_per_motif)]
  top_motifs <- names(head(abs.sum.shap_per_motif, 100))  
  return(top_motifs)}

std1 <- function(x){
  return ((x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

shap_f <- "endothelium_nc_vs_other_mouseCisBP_qval0.5_532T_shap.txt"
train_f <- "endothelium_nc_vs_other_mouseCisBP_qval0.5_532T_trainSet.csv"
data_path <- "/data/mouse/e8.25/"
data_f <- paste0(data_path, "endothelium_nc_vs_other_mouseCisBP_qval0.5.txt")

n <- 20

shap.df <- read.table(file = shap_f, header = F, stringsAsFactors = F)
train_counts.df <- read.csv(file = train_f, header = T, stringsAsFactors = F, row.names = 1)

top_motifs <- top_shap(shap.df, train_counts.df)
top_motifs <- head(top_motifs, n)

counts <- read.table(file = data_f, header = T, stringsAsFactors = F, sep = '\t')

# keep the motifs that are in the top n
counts <- counts[,c(top_motifs, "binary_celltype")] #21631    20

motifs.std <- counts
motifs.std$peakid <- rownames(motifs.std)
motifs.std$peakid <- paste(motifs.std$peakid, motifs.std$binary_celltype, sep="_")
motifs.std$binary_celltype <- NULL
#reshape counts
motifs.std <- reshape2::melt(motifs.std)#432620      3
#Using peakid as id variables

colnames(motifs.std) <-c("peakid","motif", "count")
motifs.std <- as.data.table(motifs.std)
motifs.std[, count_std := std1(count), by = "motif"]

motifs.std$count <- NULL
motifs.std <- tidyr::spread(motifs.std, motif, count_std)
motifs.std <- as.data.frame(motifs.std)#3528   38
rownames(motifs.std) <- motifs.std$peakid
motifs.std$peakid <- NULL

annot <- data.frame(class = sub(".*_", "", rownames(motifs.std)))
annot$class <- factor(annot$class, levels = c("1", "0"))
rownames(annot) <- rownames(motifs.std)

motifs.std$class <- sub(".*_", "", rownames(motifs.std))

motifs.std <- motifs.std[with(motifs.std, order(class)),]
motifs.std$class <- NULL
annot <- data.frame(class = sub(".*_", "", rownames(motifs.std)))
annot$class <- factor(annot$class, levels = c("1", "0"))
rownames(annot) <- rownames(motifs.std)

my_palette3 <- colorRampPalette(c("Light Blue", "red"))

pdf("endothelial_top_SHAP_normCount_top20.3.pdf")
pheatmap(
  mat               = motifs.std,
  color             = my_palette3(10),
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  # annotation_col    = mat_col,
  # annotation_colors = ann_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  annotation_row = annot,
  cluster_rows = F, cluster_cols = T)
dev.off()
