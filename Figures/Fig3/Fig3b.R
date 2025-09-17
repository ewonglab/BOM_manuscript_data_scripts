library("yardstick")
library("cvAUC")
library(pROC)
library(ggplot2)
library(reshape2)

f1 <- function(pred_tab){
  TP <- nrow(pred_tab[pred_tab$actual == 1 & pred_tab$predicted == 1, ])
  TN <- nrow(pred_tab[pred_tab$actual == 0 & pred_tab$predicted == 0, ])
  FP <- nrow(pred_tab[pred_tab$actual == 0 & pred_tab$predicted == 1, ])
  FN <- nrow(pred_tab[pred_tab$actual == 1 & pred_tab$predicted == 0, ])
  recall <- TP/(TP + FN)
  precision <- TP/(TP + FP)
  f1 <- 2*((precision*recall)/(precision+recall))
  f1 <- round(f1,4)
  # f1 <- round(recall, 4)
  return(f1)
}

acc <- function(pred_tab){
  acc <- nrow(pred_tab[pred_tab$predicted==pred_tab$actual,])/nrow(pred_tab)
  acc <- round(acc, 4)
  return(acc)
}

roc_auc <- function(pred_tab){
  roc.auc <- AUC(pred_tab$raw, pred_tab$actual)
  roc.auc <- round(roc.auc,4)
  return(roc.auc)
}

prauc <- function(pred_tab){
  pred_tab$actual <- factor(pred_tab$actual, levels=c(0,1))
  pr_auc_val <- pr_auc(pred_tab, truth=actual, raw, event_level="second", estimator="binary")
  pr_auc_val <- round(pr_auc_val$.estimate,4)
  return(pr_auc_val)
}

## cell type model predicting on test set
pred_files <- list.files(path = "."
                         , pattern = "(.*)_nc_vs_other_vert5.0_qval0.5_(.*)T_pred.txt")

cell_lines <- sub(pattern = "_nc_vs_other_vert5.0_qval0.5_.*", "", pred_files)

pred_tabs <- lapply(pred_files, read.table, header=T, stringsAsFactors = F)
names(pred_tabs) <- cell_lines

## calculate statistics

pred_f1 <- unlist(lapply(pred_tabs, f1))
pred_acc <- unlist(lapply(pred_tabs, acc))
pred_rocauc <- unlist(lapply(pred_tabs, roc_auc))
pred_prauc <- unlist(lapply(pred_tabs, prauc))

pred_f1 <- as.data.frame(pred_f1)
pred_acc <- as.data.frame(pred_acc)
pred_rocauc <- as.data.frame(pred_rocauc)
pred_prauc <- as.data.frame(pred_prauc)

all_stats <- merge(pred_f1, pred_acc, by=0)
all_stats <- merge(all_stats, pred_rocauc, by.x="Row.names", by.y=0)
all_stats <- merge(all_stats, pred_prauc, by.x="Row.names", by.y=0)

## read predictions on different cell types
other_stats <- read.table(file = "prediction_other_celltypes_stats", header = T
                          , stringsAsFactors = F, sep = '\t')
# bind all stats
colnames(all_stats) <- sub("Row.names", "model", colnames(all_stats))
all_stats$target <- all_stats$model
colnames(all_stats) <- sub("pred_", "", colnames(all_stats))
colnames(all_stats) <- sub("f1", "F1", colnames(all_stats))
#reorder
all_stats <- all_stats[,c("target", "model", "F1", "acc", "rocauc", "prauc")]
all_stats <- rbind(all_stats, other_stats)

F1.m <- tidyr::spread(all_stats[,c("model", "target", "F1")], target, F1)
accuracy.m <- tidyr::spread(all_stats[,c("model", "target", "acc")], target, acc)
ROCAUC.m <- tidyr::spread(all_stats[,c("model", "target", "rocauc")], target, rocauc)
PRAUC.m <- tidyr::spread(all_stats[,c("model", "target", "prauc")], target, prauc)
 
rownames(F1.m) <- F1.m$model
rownames(accuracy.m) <- accuracy.m$model
rownames(ROCAUC.m) <- ROCAUC.m$model
rownames(PRAUC.m) <- PRAUC.m$model

F1.m$model <- NULL
accuracy.m$model <- NULL
ROCAUC.m$model <- NULL
PRAUC.m$model <- NULL

model_order <- c("B", "Early.Eryth", "CD8.EM", "NK", "CD14.Mono.1", "cDC"
                 , "CMP.LMPP", "GMP", "Early.Baso", "Pre.B", "Unk_26"
                 , "CLP.1", "CD14.Mono.2", "pDC", "GMP.Neut", "CLP.2"
                 , "HSC", "CD8.CM", "Plasma", "CD4.M", "Late.Eryth"
                 , "CD4.N1", "CD4.N2", "CD8.N")

celltype_order <- c("B", "Early.Eryth", "CD8.EM", "NK", "CD14.Mono.1", "cDC"
                    , "CMP.LMPP", "GMP", "Early.Baso", "Pre.B", "Unk_26"
                    , "CLP.1", "CD14.Mono.2", "pDC", "GMP.Neut", "CLP.2"
                    , "HSC", "CD8.CM", "Plasma", "CD4.M", "Late.Eryth"
                    , "CD4.N1", "CD4.N2", "CD8.N")


F1.m <- F1.m[model_order,celltype_order]
accuracy.m <- accuracy.m[model_order,celltype_order]
ROCAUC.m <- ROCAUC.m[model_order,celltype_order]
PRAUC.m <- PRAUC.m[model_order,celltype_order]

my_breaks <- seq(0, 1, length.out = 101) # the range of statistic values is
my_color <- colorRampPalette(c("blue", "red"))(length(my_breaks) - 1)

pdf("binary_F1_otherCelltype_heatmap.pdf")
pheatmap(F1.m, color = my_color, breaks = my_breaks
         , cluster_rows = F, cluster_cols = F)
dev.off()
