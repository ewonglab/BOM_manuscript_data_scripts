library("yardstick")
library("cvAUC")
library(pROC)
library(ggplot2)
library("ggpubr")

celltype_rocauc <- function(celltype){
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
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

celltype_rocauc_CNN <- function(celltype, scores, labels){
  df <- cbind(scores[,celltype], labels[,celltype])
  colnames(df) <- c(celltype, "actual")
  df <- as.data.frame(df)
  df$predicted <- ifelse(df[,celltype] > 0.5, 1, 0)
  # print(head(df))
  # df$actual <- as.factor(df$actual)
  return(AUC(df[,celltype], df$actual))
}

celltype_prauc_CNN <- function(celltype, scores, labels){
  df <- cbind(scores[,celltype], labels[,celltype])
  colnames(df) <- c(celltype, "actual")
  df <- as.data.frame(df)
  df$predicted <- ifelse(df[,celltype] > 0.5, 1, 0)
  # print(head(df))
  colnames(df)[1] <- "CELLTYPE"
  df$actual <- factor(df$actual, levels=c(0,1))
  pr.auc <- pr_auc(df, truth=actual, CELLTYPE, event_level="second", estimator="binary")
  return(pr.auc$.estimate)
}

celltype_acc_CNN <- function(celltype, scores, labels){
  df <- cbind(scores[,celltype], labels[,celltype])
  colnames(df) <- c(celltype, "actual")
  df <- as.data.frame(df)
  df$predicted <- ifelse(df[,celltype] > 0.5, 1, 0)
  acc <- round(nrow(df[df$predicted==df$actual,])/nrow(df),4)
  return(acc)
}

combined_df <- function(BOM, DeepSTARR, DeepMEL, Basset){
  BOM <- as.data.frame(BOM)
  DeepSTARR <- as.data.frame(DeepSTARR)
  DeepMEL <- as.data.frame(DeepMEL)
  Basset <- as.data.frame(Basset)
 
  out_df <- merge(BOM, DeepSTARR, by=0)
  out_df <- merge(out_df, DeepMEL, by.x="Row.names", by.y=0)
  out_df <- merge(out_df, Basset, by.x="Row.names", by.y=0)
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

bom <- read.table(file = "test_multiclass_500bp_1687T_pred.txt", header = T
                  , stringsAsFactors = F)#2489   20

basset_label <- read.table(file = "basset_label.txt", header = F)
basset_pred <- read.table(file = "basset_pred.txt", header = F)

deepSTARR_label <- read.table(file = "deepSTARR_label.txt", header = F)
deepSTARR_pred <- read.table(file = "deepSTARR_pred.txt", header = F)

deepMEL_label <- read.table(file = "deepMEL_label.txt", header = F)
deepMEL_pred <- read.table(file = "deepMEL_pred.txt", header = F)

celltypes <- c("allantois", "cardiom" , "endothelium", "erythroid"
               , "exe_endo", "forebrain", "gut", "mesenchyme", "mid_hindbrain"
               , "mixed_meso", "neuralcrest", "NMP", "paraxial_meso"
               , "pharyngeal_meso", "somitic_meso", "spinalcord", "surface_ecto")

colnames(basset_label) <- celltypes
colnames(basset_pred) <- celltypes

colnames(deepSTARR_label) <- celltypes
colnames(deepSTARR_pred) <- celltypes

colnames(deepMEL_label) <- celltypes
colnames(deepMEL_pred) <- celltypes


# BOM STATISTICS

bom_ROCAUC <- unlist(lapply(colnames(bom)[1:17], celltype_rocauc))
bom_PRAUC <- unlist(lapply(colnames(bom)[1:17], celltype_pr_auc))
bom_acc <- unlist(lapply(colnames(bom)[1:17], celltype_acc))

names(bom_ROCAUC) <- colnames(bom)[1:17]
names(bom_PRAUC) <- colnames(bom)[1:17]
names(bom_acc) <- colnames(bom)[1:17]

# deepMel
deepSTARR_ROCAUC <- unlist(lapply(colnames(deepSTARR_pred)
                                  , function(x) celltype_rocauc_CNN(celltype = x
                                                                    , scores = deepSTARR_pred
                                                                    , labels = deepSTARR_label)))
deepSTARR_PRAUC <- unlist(lapply(colnames(deepSTARR_pred)
                                     , function(x) celltype_prauc_CNN(celltype = x
                                                                     , scores = deepSTARR_pred
                                                                     , labels = deepSTARR_label)))
deepSTARR_acc <- unlist(lapply(colnames(deepSTARR_pred)
                                     , function(x) celltype_acc_CNN(celltype = x
                                                                     , scores = deepSTARR_pred
                                                                     , labels = deepSTARR_label)))

#deepstarr
deepMEL_ROCAUC <- unlist(lapply(colnames(deepMEL_pred)
                                , function(x) celltype_rocauc_CNN(celltype = x
                                                                  , scores = deepMEL_pred
                                                                  , labels = deepMEL_label)))
deepMEL_PRAUC <- unlist(lapply(colnames(deepMEL_pred)
                                   , function(x) celltype_prauc_CNN(celltype = x
                                                                   , scores = deepMEL_pred
                                                                   , labels = deepMEL_label)))
deepMEL_acc <- unlist(lapply(colnames(deepMEL_pred)
                                   , function(x) celltype_acc_CNN(celltype = x
                                                                   , scores = deepMEL_pred
                                                                   , labels = deepMEL_label)))

# basset
basset_ROCAUC <- unlist(lapply(colnames(basset_pred)
                               , function(x) celltype_rocauc_CNN(celltype = x
                                                                 , scores = basset_pred
                                                                 , labels = basset_label)))
basset_PRAUC <- unlist(lapply(colnames(basset_pred)
                                  , function(x) celltype_prauc_CNN(celltype = x
                                                                  , scores = basset_pred
                                                                  , labels = basset_label)))
basset_acc <- unlist(lapply(colnames(basset_pred)
                                  , function(x) celltype_acc_CNN(celltype = x
                                                                  , scores = basset_pred
                                                                  , labels = basset_label)))

names(deepSTARR_ROCAUC) <- colnames(deepSTARR_pred)
names(deepSTARR_PRAUC) <- colnames(deepSTARR_pred)
names(deepSTARR_acc) <- colnames(deepSTARR_pred)

names(deepMEL_ROCAUC) <- colnames(deepSTARR_pred)
names(deepMEL_PRAUC) <- colnames(deepSTARR_pred)
names(deepMEL_acc) <- colnames(deepSTARR_pred)

names(basset_ROCAUC) <- colnames(deepSTARR_pred)
names(basset_PRAUC) <- colnames(deepSTARR_pred)
names(basset_acc) <- colnames(deepSTARR_pred)

ROCAUC_df <- combined_df(bom_ROCAUC, deepSTARR_ROCAUC, deepMEL_ROCAUC, basset_ROCAUC)
PRAUC_df <- combined_df(bom_PRAUC, deepSTARR_PRAUC, deepMEL_PRAUC, basset_PRAUC)
acc_df <- combined_df(bom_acc, deepSTARR_acc, deepMEL_acc, basset_acc)

ROCAUC_df$celltype <- rownames(ROCAUC_df)
PRAUC_df$celltype <- rownames(PRAUC_df)
acc_df$celltype <- rownames(acc_df)

library(reshape2)
ROCAUC_df <- melt(ROCAUC_df, "celltype")
PRAUC_df <- melt(PRAUC_df, "celltype")
acc_df <- melt(acc_df, "celltype")

ROCAUC_summary <- data_summary(ROCAUC_df, "value", "variable")
PRAUC_summary <- data_summary(PRAUC_df, "value", "variable")
acc_summary <- data_summary(acc_df, "value", "variable")

ROCAUC_summary$stat <- "ROCAUC"
PRAUC_summary$stat <- "PRAUC"
acc_summary$stat <- "accuracy"

all_stats <- rbind(ROCAUC_summary, PRAUC_summary, acc_summary)
all_stats$variable <- factor(all_stats$variable, levels = c("Basset", "DeepMEL", "DeepSTARR", "BOM"))

p <- ggplot(all_stats, aes(x = variable, y = value, fill = variable))+
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#CC99FF")) +
  my_theme + facet_wrap(~stat)

pdf("multinomial_stats_500bp_mean_with_errorbars.pdf")
p
dev.off()
