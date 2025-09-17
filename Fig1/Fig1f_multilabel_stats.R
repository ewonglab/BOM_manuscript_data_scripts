
# Prediction statistics of multilabel model of mouse E8.25 CREs

library("yardstick")
library("cvAUC")
library(pROC)
library(ggplot2)
library(tidyr)
library(data.table)


predictions <- readRDS("E8.25_multilab_vert5.0_qval0.5_merror_pred.rds")
predictions.df <- do.call("rbind", predictions)#55422     4

# reshape
predicted_raw.df <- tidyr::spread(predictions.df[,c("id", "model", "raw")], model, raw)#
predicted_classes.df <- tidyr::spread(predictions.df[,c("id", "predicted_class", "model")]
                                      , model, predicted_class)


# reading labels
label_filename <-  'E8.25_all_distal_nc_labels'
labels <- read.table(file = label_filename, header = T, stringsAsFactors = F
                       , sep ='\t')# 15394    19
labels_test_set <- labels[rownames(labels) %in% predictions.df$id, ]

as.data.frame(apply(labels_test_set, 2, function(x) length(x[x == 1])))
# apply(labels_test_set, 2, function(x) length(x[x == 1]))
# Allantois                                                                165
# Forebrain                                                                273
# Neural_crest                                                             209
# Somitic_mesoderm                                                         335
# Cardiomyocytes                                                           257
# Gut                                                                      206
# NMP                                                                      344
# Spinal_cord                                                              427
# Endothelium                                                              419
# Mesenchyme                                                               192
# Notochord                                                                 39
# Surface_ectoderm                                                         155
# Erythroid                                                                242
# Mid_Hindbrain                                                            295
# Paraxial_mesoderm                                                        192
# ExE_endoderm                                                              56
# Mixed_mesoderm                                                            89
# Pharyngeal_mesoderm                                                      171


# make sure data frames have same order
enh_ids <- predicted_classes.df$id
enh_ids <- enh_ids[order(enh_ids)]

rownames(predicted_classes.df) <- predicted_classes.df$id
predicted_classes.df$id <- NULL
rownames(predicted_raw.df) <- predicted_raw.df$id
predicted_raw.df$id <- NULL

predicted_classes.df <- predicted_classes.df[enh_ids, ]# [1] 3079   18
predicted_raw.df <- predicted_raw.df[enh_ids, ]# [1] 3079   18

labels <- labels[enh_ids, ]#3079   18

calculate_stats <- function(celltype){
  pred_tab <- data.frame(predicted=predicted_classes.df[,celltype]
                         , raw=predicted_raw.df[,celltype]
                         , actual = labels[,celltype])
 
  TP <- nrow(pred_tab[pred_tab$actual == 1 & pred_tab$predicted == 1, ])
  TN <- nrow(pred_tab[pred_tab$actual == 0 & pred_tab$predicted == 0, ])
  FP <- nrow(pred_tab[pred_tab$actual == 0 & pred_tab$predicted == 1, ])
  FN <- nrow(pred_tab[pred_tab$actual == 1 & pred_tab$predicted == 0, ])
  recall <- TP/(TP + FN)
  precision <- TP/(TP + FP)
  f1 <- 2*((precision*recall)/(precision+recall))
 
  ROCAUC <- round(AUC(pred_tab$raw, pred_tab$actual),4)
  acc <- round(nrow(pred_tab[pred_tab$predicted==pred_tab$actual,])/nrow(pred_tab),4)
 
  f1 <- round(f1,4)
  recall <- round(recall, 4)
  precision <- round(precision, 4)
 
  pred_tab$actual <- factor(pred_tab$actual, levels=c(0,1))
  pr.auc <- pr_auc(pred_tab, truth=actual, raw, event_level="second", estimator="binary")
  pr_auc_val <- round(pr.auc$.estimate,4)
 
  mcc_manual <- ((as.numeric(TP) * as.numeric(TN)) - (as.numeric(FP) * as.numeric(FN))) /
    sqrt((as.numeric(TP) + as.numeric(FP)) * (as.numeric(TP) + as.numeric(FN)) * (as.numeric(TN) + as.numeric(FP)) * (as.numeric(TN) + as.numeric(FN)))
  mcc_manual <- round(mcc_manual, 4)
 
  return(data.frame(celltype = celltype, ROCAUC = ROCAUC, accuracy = acc
                    , PRAUC = pr_auc_val, F1 = f1, Precision = precision
                    , Recall = recall, MCC = mcc_manual))
}


prediction_stats <- lapply(colnames(labels), calculate_stats)
prediction_stats <- do.call("rbind", prediction_stats)

write.table(x = prediction_stats, file = "E8.25_multilab_vert5.0_qval0.5_merror_pred_stats.txt", quote = F
            , sep = '\t', row.names = F)

rownames(prediction_stats) <- prediction_stats$celltype
prediction_stats$celltype <- NULL
apply(prediction_stats, 2, range, na.rm = T)
