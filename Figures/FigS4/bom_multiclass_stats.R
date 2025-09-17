library("yardstick")
library("cvAUC")
library(pROC)
library(ggplot2)
library("ggpubr")


celltype_rocauc <- function(celltype){
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
  # return((df))
  # df$actual <- as.factor(df$actual)
  return(AUC(df[,celltype], df$actual))
}

celltype_pr_auc <- function(celltype){
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
  colnames(df)[1] <- "CELLTYPE"
  df$actual <- factor(df$actual, levels=c(0,1))
  pr.auc <- pr_auc(df, truth=actual, CELLTYPE, event_level="second", estimator="binary")
  return(pr.auc$.estimate)
}


celltype_recall.BOM <- function(celltype){
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
  # colnames(df)[1] <- "CELLTYPE"
  df$predicted <- ifelse(df[,celltype] > 0.5, 1, 0)
  TP <- nrow(df[df$actual == 1 & df$predicted == 1, ])
  TN <- nrow(df[df$actual == 0 & df$predicted == 0, ])
  FP <- nrow(df[df$actual == 0 & df$predicted == 1, ])
  FN <- nrow(df[df$actual == 1 & df$predicted == 0, ])
  recall <- TP/(TP + FN)
  recall <- round(recall,4)
  return(recall)
}

celltype_prec.BOM <- function(celltype){
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
  # colnames(df)[1] <- "CELLTYPE"
  df$predicted <- ifelse(df[,celltype] > 0.5, 1, 0)
  TP <- nrow(df[df$actual == 1 & df$predicted == 1, ])
  TN <- nrow(df[df$actual == 0 & df$predicted == 0, ])
  FP <- nrow(df[df$actual == 0 & df$predicted == 1, ])
  FN <- nrow(df[df$actual == 1 & df$predicted == 0, ])
  precision <- TP/(TP + FP)
  precision <- round(precision,4)
  return(precision)
}

celltype_acc.BOM <- function(celltype){
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
  # colnames(df)[1] <- "CELLTYPE"
  df$predicted <- ifelse(df[,celltype] > 0.5, 1, 0)
  acc <- round(nrow(df[df$predicted==df$actual,])/nrow(df),4)
  return(acc)
}

celltype_MCC.BOM <- function(celltype){
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
  # colnames(df)[1] <- "CELLTYPE"
  df$predicted <- ifelse(df[,celltype] > 0.5, 1, 0)
  TP <- nrow(df[df$actual == 1 & df$predicted == 1, ])
  TN <- nrow(df[df$actual == 0 & df$predicted == 0, ])
  FP <- nrow(df[df$actual == 0 & df$predicted == 1, ])
  FN <- nrow(df[df$actual == 1 & df$predicted == 0, ])
  mcc_manual <- ((as.numeric(TP) * as.numeric(TN)) - (as.numeric(FP) * as.numeric(FN))) /
    sqrt((as.numeric(TP) + as.numeric(FP)) * (as.numeric(TP) + as.numeric(FN)) * (as.numeric(TN) + as.numeric(FP)) * (as.numeric(TN) + as.numeric(FN)))
  mcc_manual <- round(mcc_manual,4)
  return(mcc_manual)
}

celltype_F1.BOM <- function(celltype){
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
  # colnames(df)[1] <- "CELLTYPE"
  df$predicted <- ifelse(df[,celltype] > 0.5, 1, 0)
  TP <- nrow(df[df$actual == 1 & df$predicted == 1, ])
  TN <- nrow(df[df$actual == 0 & df$predicted == 0, ])
  FP <- nrow(df[df$actual == 0 & df$predicted == 1, ])
  FN <- nrow(df[df$actual == 1 & df$predicted == 0, ])
  recall <- TP/(TP + FN)
  precision <- TP/(TP + FP)
  f1 <- 2*((precision*recall)/(precision+recall))
  return(f1)
}


setwd("/g/data/zk16/cc3704/mouse_enh_grammar/xgb/e8.25")
bom <- read.table(file = "test_multiclass_500bp_4597T_pred.txt", header = T
                  , stringsAsFactors = F)


bom_recall <- unlist(lapply(colnames(bom)[1:17], celltype_recall.BOM))
bom_precision <- unlist(lapply(colnames(bom)[1:17], celltype_prec.BOM))
bom_acc <- unlist(lapply(colnames(bom)[1:17], celltype_acc.BOM))
bom_F1 <- unlist(lapply(colnames(bom)[1:17], celltype_F1.BOM))
bom_rocauc <- unlist(lapply(colnames(bom)[1:17], celltype_rocauc))
bom_prauc <- unlist(lapply(colnames(bom)[1:17], celltype_pr_auc))
bom_mcc <- unlist(lapply(colnames(bom)[1:17], celltype_MCC.BOM))

names(bom_recall) <- colnames(bom)[1:17]
names(bom_precision) <- colnames(bom)[1:17]
names(bom_acc) <- colnames(bom)[1:17]
names(bom_F1) <- colnames(bom)[1:17]
names(bom_rocauc) <- colnames(bom)[1:17]
names(bom_prauc) <- colnames(bom)[1:17]
names(bom_mcc) <- colnames(bom)[1:17]

bom_recall <- as.data.frame(bom_recall)
bom_precision <- as.data.frame(bom_precision)
bom_acc <- as.data.frame(bom_acc)
bom_F1 <- as.data.frame(bom_F1)
bom_rocauc <- as.data.frame(bom_rocauc)
bom_prauc <- as.data.frame(bom_prauc)
bom_mcc <- as.data.frame(bom_mcc)


# make a single table
multiclass_stats <- merge(bom_recall, bom_precision, by=0)
multiclass_stats <- merge(multiclass_stats, bom_acc, by.x="Row.names", by.y=0)
multiclass_stats <- merge(multiclass_stats, bom_F1, by.x="Row.names", by.y=0)
multiclass_stats <- merge(multiclass_stats, bom_rocauc, by.x="Row.names", by.y=0)
multiclass_stats <- merge(multiclass_stats, bom_prauc, by.x="Row.names", by.y=0)
multiclass_stats <- merge(multiclass_stats, bom_mcc, by.x="Row.names", by.y=0)

rownames(multiclass_stats) <- multiclass_stats$Row.names
multiclass_stats$Row.names <- NULL
multiclass_stats <- apply(multiclass_stats, 2, round, 4)

write.csv(x = multiclass_stats, file = "E8.25_multiclass_stats_500bp_mlogloss.csv"
          , quote = F, row.names = T)
