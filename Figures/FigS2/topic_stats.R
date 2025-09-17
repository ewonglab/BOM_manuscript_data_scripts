library("yardstick"))
library("cvAUC")
library(pROC)
library(ggplot2)
library("ggpubr")
library(reshape2)

celltype_rocauc <- function(celltype){
  #make a binary variable to represent actual celltype
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

celltype_acc <- function(celltype){
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
  # colnames(df)[1] <- "CELLTYPE"
  df$predicted <- ifelse(df[,celltype] > 0.5, 1, 0)
  acc <- round(nrow(df[df$predicted==df$actual,])/nrow(df),4)
  return(acc)
}

celltype_F1 <- function(celltype){
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

get_recall <- function(celltype){
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
  # colnames(df)[1] <- "CELLTYPE"
  df$predicted <- ifelse(df[,celltype] > 0.5, 1, 0)
  TP <- nrow(df[df$actual == 1 & df$predicted == 1, ])
  TN <- nrow(df[df$actual == 0 & df$predicted == 0, ])
  FP <- nrow(df[df$actual == 0 & df$predicted == 1, ])
  FN <- nrow(df[df$actual == 1 & df$predicted == 0, ])
  recall <- TP/(TP + FN)
  return(recall)
}

get_precision <- function(celltype){
  #make a binary variable to represent actual celltype
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
  # colnames(df)[1] <- "CELLTYPE"
  df$predicted <- ifelse(df[,celltype] > 0.5, 1, 0)
  TP <- nrow(df[df$actual == 1 & df$predicted == 1, ])
  TN <- nrow(df[df$actual == 0 & df$predicted == 0, ])
  FP <- nrow(df[df$actual == 0 & df$predicted == 1, ])
  FN <- nrow(df[df$actual == 1 & df$predicted == 0, ])
  precision <- TP/(TP + FP)
  return(precision)
}

get_MCC <- function(celltype){
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
  return(mcc_manual)
}


combined_stats <- function(rocauc_vals, prauc_vals
                           , acc_vals, f1_vals, prec_vals
                           , rec_vals, mcc_vals){
 
  rocauc_vals <- as.data.frame(rocauc_vals)
  prauc_vals <- as.data.frame(prauc_vals)
  acc_vals <- as.data.frame(acc_vals)
  f1_vals <- as.data.frame(f1_vals)
  prec_vals <- as.data.frame(prec_vals)
  rec_vals <- as.data.frame(rec_vals)
  mcc_vals <- as.data.frame(mcc_vals)
 
  out_df <- merge(rocauc_vals, prauc_vals, by=0)
  out_df <- merge(out_df, acc_vals, by.x="Row.names", by.y=0)
  out_df <- merge(out_df, f1_vals, by.x="Row.names", by.y=0)
  out_df <- merge(out_df, prec_vals, by.x="Row.names", by.y=0)
  out_df <- merge(out_df, rec_vals, by.x="Row.names", by.y=0)
  out_df <- merge(out_df, mcc_vals, by.x="Row.names", by.y=0)
 
  rownames(out_df) <- out_df$Row.names
  out_df$Row.names <- NULL
  return(out_df)}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func,
                    varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

my_theme <-  theme(panel.border = element_blank(), panel.grid.major = element_blank()
                   , panel.grid.minor = element_blank()#legend.position="none"
                   , axis.line = element_line(colour = "black", size = 1)
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

bom <- read.table(file = "topics_vert5.0_q0.5_9999T_pred.txt"
                  , header =T, stringsAsFactors = F)
bom <- bom[,c(paste0("V", 1:93), "label", "max_prob", "celltype")]
topic_annot <- unique(bom[,c("label", "celltype")])
topic_annot <- topic_annot[with(topic_annot, order(label)),]

n_topics <- 93
colnames(bom)[1:n_topics] <- topic_annot$celltype
bom$label <- bom$celltype


bom_ROCAUC <- unlist(lapply(colnames(bom)[1:n_topics], celltype_rocauc))
bom_PRAUC <- unlist(lapply(colnames(bom)[1:n_topics], celltype_pr_auc))
bom_acc <- unlist(lapply(colnames(bom)[1:n_topics], celltype_acc))
bom_F1 <- unlist(lapply(colnames(bom)[1:n_topics], celltype_F1))
bom_recall <- unlist(lapply(colnames(bom)[1:n_topics], get_recall))
bom_precision <- unlist(lapply(colnames(bom)[1:n_topics], get_precision))
bom_MCC <- unlist(lapply(colnames(bom)[1:n_topics], get_MCC))

names(bom_ROCAUC) <- colnames(bom)[1:n_topics]
names(bom_PRAUC) <- colnames(bom)[1:n_topics]
names(bom_acc) <- colnames(bom)[1:n_topics]
names(bom_F1) <- colnames(bom)[1:n_topics]
names(bom_recall) <- colnames(bom)[1:n_topics]
names(bom_precision) <- colnames(bom)[1:n_topics]
names(bom_MCC) <- colnames(bom)[1:n_topics]

BOM_stats_df <- combined_stats(bom_ROCAUC, bom_PRAUC, bom_acc
                               , bom_F1, bom_precision
                               , bom_recall, bom_MCC)
BOM_stats_df$Topic <- rownames(BOM_stats_df)

write.csv(x = BOM_stats_df[,c("Topic", "rocauc_vals", "prauc_vals"
                              , "acc_vals", "f1_vals", "prec_vals", "rec_vals", "mcc_vals")]
          , file = "BOM_multiclass_Topics_stats_mlogloss", quote = F)
