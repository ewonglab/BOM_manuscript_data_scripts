library("yardstick")
library("cvAUC")
library(pROC)
library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)

pred_tab <- read.table(file = args[1], header =T, stringsAsFactors = F)

TP <- nrow(pred_tab[pred_tab$actual == 1 & pred_tab$predicted == 1, ])
TN <- nrow(pred_tab[pred_tab$actual == 0 & pred_tab$predicted == 0, ])
FP <- nrow(pred_tab[pred_tab$actual == 0 & pred_tab$predicted == 1, ])
FN <- nrow(pred_tab[pred_tab$actual == 1 & pred_tab$predicted == 0, ])
recall <- TP/(TP + FN)
precision <- TP/(TP + FP)
f1 <- 2*((precision*recall)/(precision+recall))

roc_auc <- round(AUC(pred_tab$raw, pred_tab$actual),4)
acc <- round(nrow(pred_tab[pred_tab$predicted==pred_tab$actual,])/nrow(pred_tab),4)
f1 <- round(f1,4)
recall <- round(recall, 4)
precision <- round(precision, 4)

pred_tab$actual <- factor(pred_tab$actual, levels=c(0,1))
pr_auc_val <- pr_auc(pred_tab, truth=actual, raw, event_level="second", estimator="binary")
pr_auc_val <- round(pr_auc_val$.estimate,4)

celltype <- basename(args[1])
celltype <- sub("_vs.*", "", celltype)
#output dataframe
out_df <- data.frame(celltype=celltype, ROC_AUC=roc_auc, RECALL=recall, PRECISION=precision
                     , Accuracy=acc, F1=f1, PRAUC=pr_auc_val)
out_df <- melt(out_df, "celltype")

write.table(x = out_df, file = args[2], col.names = F, row.names = F, sep = '\t', append=T, quote = F)
