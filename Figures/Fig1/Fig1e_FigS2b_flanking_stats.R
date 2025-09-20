suppressMessages({
  library(data.table)
  library(tidyr)
  library(dplyr)
  library("xgboost")
})

setwd("/E8.25_flanking")

counts_matrix <- function(fimoOutDir = "/E8.25_flanking/fimo_out"
                          , input_dir, qval_thresh = 0.5){
    motif_counts <- fread(file.path(fimoOutDir, paste0(input_dir, "/fimo.tsv")))

    motif_counts <- as.data.frame(motif_counts)

    #subset by q-value (column 9)
    motif_counts <- motif_counts[motif_counts[,9] <= qval_thresh,]

    # counts
    motif_counts <- as.data.frame(table(motif_counts[,c("motif_id", "sequence_name")]))

    # reshape
    motif_counts <- tidyr::spread(motif_counts, motif_id, Freq)#

    # dim(ups_counts)
    # dim(dwn_counts)
    rownames(motif_counts) <- motif_counts$sequence_name
    
    motif_counts$sequence_name <- NULL
    return(motif_counts)
}

card_ups <- counts_matrix(input_dir = "cardiom_windows_ups_50str")
card_dwn <- counts_matrix(input_dir = "cardiom_windows_dwn_50str")
endo_ups <- counts_matrix(input_dir = "Endothelium_windows_ups_50str")
endo_dwn <- counts_matrix(input_dir = "Endothelium_windows_dwn_50str")
ery_ups <- counts_matrix(input_dir = "Erythroid_windows_ups_50str")
ery_dwn <- counts_matrix(input_dir = "Erythroid_windows_dwn_50str")
ExEendo_ups <- counts_matrix(input_dir = "ExEendoderm_windows_ups_50str")
ExEendo_dwn <- counts_matrix(input_dir = "ExEendoderm_windows_dwn_50str")
forebrain_ups <- counts_matrix(input_dir = "Forebrain_windows_ups_50str")
forebrain_dwn <- counts_matrix(input_dir = "Forebrain_windows_dwn_50str")
Gut_ups <- counts_matrix(input_dir = "Gut_windows_ups_50str")
Gut_dwn <- counts_matrix(input_dir = "Gut_windows_dwn_50str")

write.table(x = card_ups, file = "cardiom_ups_flanking_counts_q0.5.txt", quote = F, sep = '\t')
write.table(x = card_dwn, file = "cardiom_dwn_flanking_counts_q0.5.txt", quote = F, sep = '\t')

write.table(x = endo_ups, file = "Endothelium_ups_flanking_counts_q0.5.txt", quote = F, sep = '\t')
write.table(x = endo_dwn, file = "Endothelium_dwn_flanking_counts_q0.5.txt", quote = F, sep = '\t')

write.table(x = ery_ups, file = "Erythroid_ups_flanking_counts_q0.5.txt", quote = F, sep = '\t')
write.table(x = ery_dwn, file = "Erythroid_dwn_flanking_counts_q0.5.txt", quote = F, sep = '\t')

write.table(x = ExEendo_ups, file = "ExEendoderm_ups_flanking_counts_q0.5.txt", quote = F, sep = '\t')
write.table(x = ExEendo_dwn, file = "ExEendoderm_dwn_flanking_counts_q0.5.txt", quote = F, sep = '\t')

write.table(x = forebrain_ups, file = "Forebrain_ups_flanking_counts_q0.5.txt", quote = F, sep = '\t')
write.table(x = forebrain_dwn, file = "Forebrain_dwn_flanking_counts_q0.5.txt", quote = F, sep = '\t')

write.table(x = Gut_ups, file = "Gut_ups_flanking_counts_q0.5.txt", quote = F, sep = '\t')
write.table(x = Gut_dwn, file = "Gut_dwn_flanking_counts_q0.5.txt", quote = F, sep = '\t')

add.missing.vars_xgb <- function(xgb.model , data){
  missing.vars <- setdiff(xgb.model$feature_names, colnames(data))
  if(length(missing.vars)>0){
    for(i in 1:length(missing.vars)){
      data[,missing.vars[i]] <- 0}}
  return(data)}

modelPath <- "/xgb/e8.25/bin_500bp/"
cardiom.xgb <- readRDS(file.path(modelPath, "cardiom_500bp_vert5.0_qval0.5.rds"))
endothelium.xgb <- readRDS(file.path(modelPath, "endothelium_500bp_vert5.0_qval0.5.rds"))
erythroid.xgb <- readRDS(file.path(modelPath, "erythroid_500bp_vert5.0_qval0.5.rds"))
exeEndo.xgb <- readRDS(file.path(modelPath, "exe_endo_500bp_vert5.0_qval0.5.rds"))
forebrain.xgb <- readRDS(file.path(modelPath, "forebrain_500bp_vert5.0_qval0.5.rds"))
gut.xgb <- readRDS(file.path(modelPath, "gut_500bp_vert5.0_qval0.5.rds"))

card_ups <- read.table(file = "cardiom_ups_flanking_counts_q0.5.txt", header = TRUE, stringsAsFactors = FALSE
                       , sep = '\t')
card_dwn <- read.table(file = "cardiom_dwn_flanking_counts_q0.5.txt", header = TRUE, stringsAsFactors = FALSE
                       , sep = '\t')

endo_ups <- read.table(file = "Endothelium_ups_flanking_counts_q0.5.txt", header = TRUE, stringsAsFactors = FALSE
                       , sep = '\t')
endo_dwn <- read.table(file = "Endothelium_dwn_flanking_counts_q0.5.txt", header = TRUE, stringsAsFactors = FALSE
                       , sep = '\t')

ery_ups <- read.table(file = "Erythroid_ups_flanking_counts_q0.5.txt", header = TRUE, stringsAsFactors = FALSE
                      , sep = '\t')
ery_dwn <- read.table(file = "Erythroid_dwn_flanking_counts_q0.5.txt", header = TRUE, stringsAsFactors = FALSE
                      , sep = '\t')

ExEendo_ups <- read.table(file = "ExEendoderm_ups_flanking_counts_q0.5.txt", header = TRUE
                          , stringsAsFactors = FALSE, sep = '\t')
ExEendo_dwn <- read.table(file = "ExEendoderm_dwn_flanking_counts_q0.5.txt", header = TRUE
                          , stringsAsFactors = FALSE, sep = '\t')

forebrain_ups <- read.table(file = "Forebrain_ups_flanking_counts_q0.5.txt", header = TRUE
                            , stringsAsFactors = FALSE, sep = '\t')
forebrain_dwn <- read.table(file = "Forebrain_dwn_flanking_counts_q0.5.txt", header = TRUE
                            , stringsAsFactors = FALSE, sep = '\t')

Gut_ups <- read.table(file = "Gut_ups_flanking_counts_q0.5.txt", header = TRUE, stringsAsFactors = FALSE
                      , sep = '\t')
Gut_dwn <- read.table(file = "Gut_dwn_flanking_counts_q0.5.txt", header = TRUE, stringsAsFactors = FALSE
                      , sep = '\t')


# cardiom
card_ups <- add.missing.vars_xgb(cardiom.xgb, card_ups)
card_ups <- card_ups[,cardiom.xgb$feature_names]  

card_dwn <- add.missing.vars_xgb(cardiom.xgb, card_dwn)
card_dwn <- card_dwn[,cardiom.xgb$feature_names]

# endothelium
endo_ups <- add.missing.vars_xgb(endothelium.xgb, endo_ups)
endo_ups <- endo_ups[,endothelium.xgb$feature_names]  

endo_dwn <- add.missing.vars_xgb(endothelium.xgb, endo_dwn)
endo_dwn <- endo_dwn[,endothelium.xgb$feature_names]

# erythroid
ery_ups <- add.missing.vars_xgb(erythroid.xgb, ery_ups)
ery_ups <- ery_ups[,erythroid.xgb$feature_names]  

ery_dwn <- add.missing.vars_xgb(erythroid.xgb, ery_dwn)
ery_dwn <- ery_dwn[,erythroid.xgb$feature_names]

# exe endo
ExEendo_ups <- add.missing.vars_xgb(exeEndo.xgb, ExEendo_ups)
ExEendo_ups <- ExEendo_ups[,exeEndo.xgb$feature_names]  

ExEendo_dwn <- add.missing.vars_xgb(exeEndo.xgb, ExEendo_dwn)
ExEendo_dwn <- ExEendo_dwn[,exeEndo.xgb$feature_names]

# forebrain
forebrain_ups <- add.missing.vars_xgb(forebrain.xgb, forebrain_ups)
forebrain_ups <- forebrain_ups[,forebrain.xgb$feature_names]  

forebrain_dwn <- add.missing.vars_xgb(forebrain.xgb, forebrain_dwn)
forebrain_dwn <- forebrain_dwn[,forebrain.xgb$feature_names]

# gut
Gut_ups <- add.missing.vars_xgb(gut.xgb, Gut_ups)
Gut_ups <- Gut_ups[,gut.xgb$feature_names]  

Gut_dwn <- add.missing.vars_xgb(gut.xgb, Gut_dwn)
Gut_dwn <- Gut_dwn[,gut.xgb$feature_names]


predict_windows <- function(model, ups_counts, dwn_counts){
    set.seed(123)
    y_pred <- predict(model, data.matrix(ups_counts), type = "response")
    predicted.class <- y_pred > 0.5 
    predicted.class <- gsub("TRUE", 1, predicted.class)
    predicted.class <- gsub("FALSE", 0, predicted.class)
    ups_predictions.df <- data.frame(predicted = predicted.class
                                 , raw = y_pred
                                 , stringsAsFactors = FALSE)
    rownames(ups_predictions.df) <- rownames(ups_counts)
    
    set.seed(123)
    y_pred <- predict(model, data.matrix(dwn_counts), type = "response")
    predicted.class <- y_pred > 0.5 
    predicted.class <- gsub("TRUE", 1, predicted.class)
    predicted.class <- gsub("FALSE", 0, predicted.class)
    dws_predictions.df <- data.frame(predicted = predicted.class
                                 , raw = y_pred
                                 , stringsAsFactors = FALSE)
    rownames(dws_predictions.df) <- rownames(dwn_counts)
    ups_predictions.df$type <- "upstream"
    dws_predictions.df$type <- "downstream"
    all_pred <- rbind(ups_predictions.df, dws_predictions.df)
    all_pred$actual <- 0 # add actual class
    return(all_pred)
} 

card_windows_pred <- predict_windows(model = cardiom.xgb, ups_counts = card_ups
                                     , dwn_counts = card_dwn)
endothelium_windows_pred <- predict_windows(model = endothelium.xgb, ups_counts = endo_ups
                                            , dwn_counts = endo_dwn)
erythroid_windows_pred <- predict_windows(model = erythroid.xgb, ups_counts = ery_ups
                                          , dwn_counts = ery_dwn)
exeEndo_windows_pred <- predict_windows(model = exeEndo.xgb, ups_counts = ExEendo_ups
                                        , dwn_counts = ExEendo_dwn)
forebrain_windows_pred <- predict_windows(model = forebrain.xgb, ups_counts = forebrain_ups
                                          , dwn_counts = forebrain_dwn)
gut_windows_pred <- predict_windows(model = gut.xgb, ups_counts = Gut_ups
                                    , dwn_counts = Gut_dwn)

suppressMessages({
    library("yardstick")
    library("cvAUC")
    library(pROC)
    library(ggplot2)
    library("mltools")
})


cardiom_test <- read.table(file = file.path(modelPath, "cardiom_500bp_vert5.0_qval0.5_578T_pred.txt")
                           , header = TRUE, stringsAsFactors = FALSE)
endoth_test <- read.table(file = file.path(modelPath, "endothelium_500bp_vert5.0_qval0.5_12T_pred.txt")
                           , header = TRUE, stringsAsFactors = FALSE)
ery_test <- read.table(file = file.path(modelPath, "erythroid_500bp_vert5.0_qval0.5_20T_pred.txt")
                           , header = TRUE, stringsAsFactors = FALSE)
exeEndo_test <- read.table(file = file.path(modelPath, "exe_endo_500bp_vert5.0_qval0.5_24T_pred.txt")
                           , header = TRUE, stringsAsFactors = FALSE)
forebr_test <- read.table(file = file.path(modelPath, "forebrain_500bp_vert5.0_qval0.5_330T_pred.txt")
                           , header = TRUE, stringsAsFactors = FALSE)
gut_test <- read.table(file = file.path(modelPath, "gut_500bp_vert5.0_qval0.5_450T_pred.txt")
                           , header = TRUE, stringsAsFactors = FALSE)


subset_positive <- function(x){
    x <- x[x$actual == 1, ]
    # print(dim(x))
    x$type <- "CRE"
    x <- x[, c("predicted", "raw", "type", "actual")]
    return(x)
}
cardiom_test <- subset_positive(cardiom_test)
endoth_test <- subset_positive(endoth_test)
ery_test <- subset_positive(ery_test)
exeEndo_test <- subset_positive(exeEndo_test)
forebr_test <- subset_positive(forebr_test)
gut_test <- subset_positive(gut_test)

all_pred_cardiom <- rbind(card_windows_pred, cardiom_test)
all_pred_endoth <- rbind(endothelium_windows_pred, endoth_test)
all_pred_ery <- rbind(erythroid_windows_pred, ery_test)
all_pred_exeEndo <- rbind(exeEndo_windows_pred, exeEndo_test)
all_pred_forebr <- rbind(forebrain_windows_pred, forebr_test)
all_pred_gut <- rbind(gut_windows_pred, gut_test)

pred_stats <- function(pred_tab){
    pred_tab$actual <- as.integer(pred_tab$actual)
    pred_tab$predicted <- as.integer(pred_tab$predicted)
    TP <- nrow(pred_tab[pred_tab$actual == 1 & pred_tab$predicted == 1, ])
    TN <- nrow(pred_tab[pred_tab$actual == 0 & pred_tab$predicted == 0, ])
    FP <- nrow(pred_tab[pred_tab$actual == 0 & pred_tab$predicted == 1, ])
    FN <- nrow(pred_tab[pred_tab$actual == 1 & pred_tab$predicted == 0, ])
    recall <- TP/(TP + FN)
    precision <- TP/(TP + FP)
    specificity <- TN / (TN + FP)
    FPR <- FP/(FP + TN)
    f1 <- 2*((precision*recall)/(precision+recall))
    mc <- mcc(pred_tab$predicted, pred_tab$actual)
    roc_auc <- AUC(pred_tab$raw, pred_tab$actual) # 2 actual classes needed
    acc <- nrow(pred_tab[pred_tab$predicted == pred_tab$actual,])/nrow(pred_tab)
    pred_tab$actual <- factor(pred_tab$actual, levels = c(0,1))
    pr_auc_val <- pr_auc(pred_tab, truth = actual, raw, event_level = "second", estimator = "binary")
    STATS <- data.frame(recall = round(recall, 3), precision = round(precision, 3)
                        , f1 = round(f1, 3), mc = round(mc, 3), specificity = round(specificity, 3) #, roc_auc = round(roc_auc, 3)
                        , acc = round(acc, 3), FPR = round(FPR, 3), roc_auc = round(roc_auc, 3)
                        , pr_auc_val = round(pr_auc_val$.estimate, 3) 
                        , TP = TP, TN = TN, FP = FP, FN = FN)#, pr_auc = round(pr_auc_val, 3))
    return(STATS)
}

card_stats <- pred_stats(all_pred_cardiom)
endoth_stats <- pred_stats(all_pred_endoth)
ery_stats <- pred_stats(all_pred_ery)
exeEndo_stats <- pred_stats(all_pred_exeEndo)
forebr_stats <- pred_stats(all_pred_forebr)
gut_stats <- pred_stats(all_pred_gut)

# add cell types to all predictions
all_pred_cardiom$Model <- "cardiomyocyte"
all_pred_endoth$Model <- "endothelium"
all_pred_ery$Model <- "erythroid"
all_pred_exeEndo$Model <- "exe endoderm"
all_pred_forebr$Model <- "forebrain"
all_pred_gut$Model <- "gut"

all_pred_cardiom$Region_ID <- rownames(all_pred_cardiom)
all_pred_endoth$Region_ID <- rownames(all_pred_endoth)
all_pred_ery$Region_ID <- rownames(all_pred_ery)
all_pred_exeEndo$Region_ID <- rownames(all_pred_exeEndo)
all_pred_forebr$Region_ID <- rownames(all_pred_forebr)
all_pred_gut$Region_ID <- rownames(all_pred_gut)


all_flanking_pred <- rbind(all_pred_cardiom, all_pred_endoth
                           , all_pred_ery, all_pred_exeEndo
                           , all_pred_forebr, all_pred_gut)

rownames(all_flanking_pred) <- 1:nrow(all_flanking_pred)
all_flanking_pred <- all_flanking_pred[, c("Model", "Region_ID", "type", "actual", "raw", "predicted")]
colnames(all_flanking_pred) <- sub("raw", "prob", colnames(all_flanking_pred))
write.csv(x = all_flanking_pred, file = "all_flanking_predictions.csv", quote = TRUE)

card_stats$celltype <- "cardiom"
endoth_stats$celltype <- "endothelium"
ery_stats$celltype <- "erythroid"
exeEndo_stats$celltype <- "exeEndo"
forebr_stats$celltype <- "forebrain"
gut_stats$celltype <- "gut"

all_stats <- rbind(card_stats, endoth_stats, ery_stats, exeEndo_stats, forebr_stats, gut_stats)
write.csv(x = all_stats, file = "flanking_regions_stats.csv", quote = FALSE)

breaksList = seq(0, 1, by = 0.1)
pheatmap(all_stats[, c("recall", "precision", "f1", "mc", "specificity", "FPR", "acc")]
         , scale = "none", cluster_rows = FALSE
         , cluster_cols = FALSE
         # , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))
         , color = colorRampPalette(c("white", "red"))(length(breaksList))
         , angle_col = 90, breaks = breaksList, display_numbers = TRUE, fontsize = 14)
