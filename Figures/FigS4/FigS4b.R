
# scatterplots for multiclass models - BOM versus best CNN models from Xuan (DeepMEL, BASSET and DeepSTARR) - Precision, recall (and MCC)

library("yardstick")
library("cvAUC")
library(pROC)
library(ggplot2)
library("ggpubr")

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

bom <- read.table(file = "test_multiclass_500bp_4597T_pred.txt", header = T
                  , stringsAsFactors = F)#2489   20

best_models_path <- "/best_models/"

basset_label <- read.table(file = paste0(best_models_path, "basset_test_labels_BEST.txt")
                           , header = F)
basset_pred <- read.table(file = paste0(best_models_path, "basset_test_scores_BEST.txt")
                          , header = F)

deepSTARR_label <- read.table(file = paste0(best_models_path, "deepstar_test_labels_BEST.txt")
                              , header = F)
deepSTARR_pred <- read.table(file = paste0(best_models_path, "deepstar_test_scores_BEST.txt")
                             , header = F)

deepMEL_label <- read.table(file = paste0(best_models_path, "deepmel_test_labels_BEST.txt")
                            , header = F)
deepMEL_pred <- read.table(file = paste0(best_models_path, "deepmel_test_scores_BEST.txt")
                           , header = F)

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

celltype_recall.BOM <- function(celltype){
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
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
  df$predicted <- ifelse(df[,celltype] > 0.5, 1, 0)
  TP <- nrow(df[df$actual == 1 & df$predicted == 1, ])
  TN <- nrow(df[df$actual == 0 & df$predicted == 0, ])
  FP <- nrow(df[df$actual == 0 & df$predicted == 1, ])
  FN <- nrow(df[df$actual == 1 & df$predicted == 0, ])
  precision <- TP/(TP + FP)
  precision <- round(precision,4)
  return(precision)
}

celltype_MCC.BOM <- function(celltype){
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
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

bom_recall <- unlist(lapply(colnames(bom)[1:17], celltype_recall.BOM))
bom_precision <- unlist(lapply(colnames(bom)[1:17], celltype_prec.BOM))
bom_MCC <- unlist(lapply(colnames(bom)[1:17], celltype_MCC.BOM))

names(bom_recall) <- colnames(bom)[1:17]
names(bom_precision) <- colnames(bom)[1:17]
names(bom_MCC) <- colnames(bom)[1:17]

celltype_recall.CNN <- function(celltype, scores, labels){
  df <- cbind(scores[,celltype], labels[,celltype])
  colnames(df) <- c(celltype, "actual")
  df <- as.data.frame(df)
  df$predicted <- ifelse(df[,celltype] > 0.5, 1, 0)
  TP <- nrow(df[df$actual == 1 & df$predicted == 1, ])
  TN <- nrow(df[df$actual == 0 & df$predicted == 0, ])
  FP <- nrow(df[df$actual == 0 & df$predicted == 1, ])
  FN <- nrow(df[df$actual == 1 & df$predicted == 0, ])
  recall <- TP/(TP + FN)
  recall <- round(recall,4)
  return(recall)
}

celltype_prec.CNN <- function(celltype, scores, labels){
  df <- cbind(scores[,celltype], labels[,celltype])
  colnames(df) <- c(celltype, "actual")
  df <- as.data.frame(df)
  df$predicted <- ifelse(df[,celltype] > 0.5, 1, 0)
  TP <- nrow(df[df$actual == 1 & df$predicted == 1, ])
  TN <- nrow(df[df$actual == 0 & df$predicted == 0, ])
  FP <- nrow(df[df$actual == 0 & df$predicted == 1, ])
  FN <- nrow(df[df$actual == 1 & df$predicted == 0, ])
  # print(celltype)
  # print(TP)
  # print(FP)
  precision <- TP/(TP + FP)
  precision <- round(precision,4)
  return(precision)
}

celltype_MCC.CNN <- function(celltype, scores, labels){
  df <- cbind(scores[,celltype], labels[,celltype])
  colnames(df) <- c(celltype, "actual")
  df <- as.data.frame(df)
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

# deepMel
deepSTARR_recall <- unlist(lapply(colnames(deepSTARR_pred)
                                  , function(x) celltype_recall.CNN(celltype = x
                                                                    , scores = deepSTARR_pred
                                                                    , labels = deepSTARR_label)))
deepSTARR_precision <- unlist(lapply(colnames(deepSTARR_pred)
                                     , function(x) celltype_prec.CNN(celltype = x
                                                                     , scores = deepSTARR_pred
                                                                     , labels = deepSTARR_label)))
deepSTARR_MCC <- unlist(lapply(colnames(deepSTARR_pred)
                                     , function(x) celltype_MCC.CNN(celltype = x
                                                                     , scores = deepSTARR_pred
                                                                     , labels = deepSTARR_label)))

#deepstarr
deepMEL_recall <- unlist(lapply(colnames(deepMEL_pred)
                                , function(x) celltype_recall.CNN(celltype = x
                                                                  , scores = deepMEL_pred
                                                                  , labels = deepMEL_label)))
deepMEL_precision <- unlist(lapply(colnames(deepMEL_pred)
                                   , function(x) celltype_prec.CNN(celltype = x
                                                                   , scores = deepMEL_pred
                                                                   , labels = deepMEL_label)))
deepMEL_MCC <- unlist(lapply(colnames(deepMEL_pred)
                                   , function(x) celltype_MCC.CNN(celltype = x
                                                                   , scores = deepMEL_pred
                                                                   , labels = deepMEL_label)))

# basset
basset_recall <- unlist(lapply(colnames(basset_pred)
                               , function(x) celltype_recall.CNN(celltype = x
                                                                 , scores = basset_pred
                                                                 , labels = basset_label)))
basset_precision <- unlist(lapply(colnames(basset_pred)
                                  , function(x) celltype_prec.CNN(celltype = x
                                                                  , scores = basset_pred
                                                                  , labels = basset_label)))
basset_MCC <- unlist(lapply(colnames(basset_pred)
                                  , function(x) celltype_MCC.CNN(celltype = x
                                                                  , scores = basset_pred
                                                                  , labels = basset_label)))

# same colnames
# colnames(deepSTARR_pred)
# [1] "allantois"       "cardiom"         "endothelium"     "erythroid"      
# [5] "exe_endo"        "forebrain"       "gut"             "mesenchyme"    
# [9] "mid_hindbrain"   "mixed_meso"      "neuralcrest"     "NMP"            
# [13] "paraxial_meso"   "pharyngeal_meso" "somitic_meso"    "spinalcord"    
# [17] "surface_ecto"  
# > colnames(basset_pred)
# [1] "allantois"       "cardiom"         "endothelium"     "erythroid"      
# [5] "exe_endo"        "forebrain"       "gut"             "mesenchyme"    
# [9] "mid_hindbrain"   "mixed_meso"      "neuralcrest"     "NMP"            
# [13] "paraxial_meso"   "pharyngeal_meso" "somitic_meso"    "spinalcord"    
# [17] "surface_ecto"  
# > colnames(deepMEL_pred)
# [1] "allantois"       "cardiom"         "endothelium"     "erythroid"      
# [5] "exe_endo"        "forebrain"       "gut"             "mesenchyme"    
# [9] "mid_hindbrain"   "mixed_meso"      "neuralcrest"     "NMP"            
# [13] "paraxial_meso"   "pharyngeal_meso" "somitic_meso"    "spinalcord"    
# [17] "surface_ecto"

names(deepSTARR_recall) <- colnames(deepSTARR_pred)
names(deepSTARR_precision) <- colnames(deepSTARR_pred)
names(deepMEL_recall) <- colnames(deepSTARR_pred)
names(deepMEL_precision) <- colnames(deepSTARR_pred)
names(basset_recall) <- colnames(deepSTARR_pred)
names(basset_precision) <- colnames(deepSTARR_pred)

names(deepSTARR_MCC) <- colnames(deepSTARR_pred)
names(deepMEL_MCC) <- colnames(deepSTARR_pred)
names(basset_MCC) <- colnames(deepSTARR_pred)

combined_df <- function(BOM, DeepSTARR, DeepMEL, Basset){
  BOM <- as.data.frame(BOM)
  DeepSTARR <- as.data.frame(DeepSTARR)
  DeepMEL <- as.data.frame(DeepMEL)
  Basset <- as.data.frame(Basset)
 
  # print(setdiff(rownames(BOM), rownames(DeepSTARR)))
  out_df <- merge(BOM, DeepSTARR, by=0)
  # print(dim(out_df))
  out_df <- merge(out_df, DeepMEL, by.x="Row.names", by.y=0)
  # print(dim(out_df))
  out_df <- merge(out_df, Basset, by.x="Row.names", by.y=0)
  # print(dim(out_df))
  rownames(out_df) <- out_df$Row.names
  out_df$Row.names <- NULL
  return(out_df)}

recall_df <- combined_df(bom_recall, deepSTARR_recall, deepMEL_recall, basset_recall)
precision_df <- combined_df(bom_precision, deepSTARR_precision
                            , deepMEL_precision, basset_precision)
MCC_df <- combined_df(bom_MCC, deepSTARR_MCC, deepMEL_MCC, basset_MCC)

recall_df$celltype <- rownames(recall_df)
precision_df$celltype <- rownames(precision_df)
MCC_df$celltype <- rownames(MCC_df)

dot_colors <- data.frame(celltype = c("allantois", "mixed_meso", "cardiom"
                                      , "neuralcrest", "endothelium", "NMP"
                                      , "erythroid", "paraxial_meso"
                                      , "exe_endo", "pharyngeal_meso"
                                      , "forebrain", "somitic_meso"
                                      , "gut", "spinalcord", "mesenchyme"
                                      , "surface_ecto", "mid_hindbrain")
                         , colors = c("#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442"
                                                    , "#56B4E9", "#E69F00", "#009292", "#E41A1C", "#F781BF", "#66C2A5"
                                                    , "#FC8D62", "#A6D854", "#5A4A7F", "#F4A460", "#4682B4", "#E6B8B7"))

dot_colors <- dot_colors[with(dot_colors, order(celltype)),]


p1 <- ggplot(recall_df, aes(x=DeepSTARR, y=BOM, color=celltype)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray") + my_theme +
  scale_x_continuous(limits = c(0, 1)) + #geom_text(hjust=0, vjust=0, size = 4) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values=dot_colors$colors, aesthetics = c("colour", "fill"))

p2 <- ggplot(recall_df, aes(x=DeepMEL, y=BOM, color=celltype)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray") + my_theme + #geom_text(hjust=0, vjust=0, size = 4) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values=dot_colors$colors, aesthetics = c("colour", "fill"))

p3 <- ggplot(recall_df, aes(x=Basset, y=BOM, color=celltype)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray") + my_theme + #geom_text(hjust=0, vjust=0, size = 4) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values=dot_colors$colors, aesthetics = c("colour", "fill"))

pdf("multiclass_recall_BEST.pdf", onefile = F)
ggarrange(plotlist = list(p1, p2, p3), ncol=3, nrow=1, common.legend =T)
dev.off()

p1 <- ggplot(precision_df, aes(x=DeepSTARR, y=BOM, color=celltype)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray") + my_theme +
  scale_x_continuous(limits = c(0, 1)) + #geom_text(hjust=0, vjust=0, size = 4) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values=dot_colors$colors, aesthetics = c("colour", "fill"))

p2 <- ggplot(precision_df, aes(x=DeepMEL, y=BOM, color=celltype)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray") + my_theme + #geom_text(hjust=0, vjust=0, size = 4) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values=dot_colors$colors, aesthetics = c("colour", "fill"))

p3 <- ggplot(precision_df, aes(x=Basset, y=BOM, color=celltype)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray") + my_theme + #geom_text(hjust=0, vjust=0, size = 4) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values=dot_colors$colors, aesthetics = c("colour", "fill"))


pdf("multiclass_precision_BEST.pdf", onefile = F)
ggarrange(plotlist = list(p1, p2, p3), ncol=3, nrow=1, common.legend =T)
dev.off()


p1 <- ggplot(MCC_df, aes(x=DeepSTARR, y=BOM, color=celltype)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray") + my_theme +
  scale_x_continuous(limits = c(0, 1)) + #geom_text(hjust=0, vjust=0, size = 4) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values=dot_colors$colors, aesthetics = c("colour", "fill"))

p2 <- ggplot(MCC_df, aes(x=DeepMEL, y=BOM, color=celltype)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray") + my_theme + #geom_text(hjust=0, vjust=0, size = 4) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values=dot_colors$colors, aesthetics = c("colour", "fill"))

p3 <- ggplot(MCC_df, aes(x=Basset, y=BOM, color=celltype)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray") + my_theme + #geom_text(hjust=0, vjust=0, size = 4) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values=dot_colors$colors, aesthetics = c("colour", "fill"))

pdf("multiclass_MCC_BEST.pdf", onefile = F)
ggarrange(plotlist = list(p1, p2, p3), ncol=3, nrow=1, common.legend =T)
dev.off()
                            
write.table(x = recall_df, file = "recall_per_celltype_multiclass_BEST"
            , quote = F, sep ='\t')
write.table(x = precision_df, file = "precision_per_celltype_multiclass_BEST"
            , quote = F, sep ='\t')
write.table(x = MCC_df, file = "MCC_per_celltype_multiclass_BEST"
            , quote = F, sep ='\t')

