
suppressMessages({
  library("rsample")
  library(xgboost)
  library(readr)
  library(stringr)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(RColorBrewer)
  library(gridExtra)
  library(grid)
})

args <- commandArgs(trailingOnly=TRUE)

data_path <- "/data/human/blood/granja/"
model_celltype <- args[1]
target_celltype <- args[2]

xgb_model <- paste0(model_celltype, "_nc_vs_other_vert5.0_qval0.5_valE.rds")
model_data <- paste0(model_celltype, "_nc_vs_other_vert5.0_qval0.5.txt")
test_data <- paste0(target_celltype, "_nc_vs_other_vert5.0_qval0.5.txt")

add.missing.vars_xgb <- function(xgb.model , data){
  missing.vars <- setdiff(xgb.model$feature_names, colnames(data))
  if(length(missing.vars)>0){
    for(i in 1:length(missing.vars)){
      data[,missing.vars[i]] <- 0}}
  return(data)}

xgb <- readRDS(xgb_model)

counts.tab <- read.table(file = paste0(data_path, model_data), header =T, stringsAsFactors = F, sep = '\t')
counts.tab$celltype <- NULL

if("sequence_name" %in% colnames(counts.tab)){
  rownames(counts.tab) <- counts.tab$sequence_name
  counts.tab$sequence_name <- NULL
}

counts.tab.NAs <- sapply(counts.tab, function(x) sum(is.na(x)))
print(unique(counts.tab.NAs))#0
counts.tab$binary_celltype <- as.numeric(counts.tab$binary_celltype)

set.seed(123)
motifs_split <- initial_split(counts.tab, prop = .6)
motifs_train <- training(motifs_split)
motifs_test.0 <- testing(motifs_split)

set.seed(123)
motifs_split2 <- initial_split(motifs_test.0, prop = .5)
motifs_val <- training(motifs_split2)
motifs_test.0 <- testing(motifs_split2)

motifs_test <- read.table(file = paste0(data_path, test_data), header =T, stringsAsFactors = F, sep = '\t')
motifs_test <- motifs_test[!rownames(motifs_test) %in% c(rownames(motifs_train), rownames(motifs_val)),]

motifs_test <- add.missing.vars_xgb(xgb, motifs_test)
test_labels <- motifs_test$celltype
test_bin_labels <- motifs_test$binary_celltype
motifs_test <- motifs_test[,xgb$feature_names]  
set.seed(123)
y_pred <- predict(xgb, data.matrix(motifs_test), type="response")
predicted.class <- y_pred > 0.5 
predicted.class <- gsub("TRUE", 1, predicted.class)
predicted.class <- gsub("FALSE", 0, predicted.class)
actual.vs.predicted <- data.frame(actual_ct = test_labels
                                  , actual = test_bin_labels
                                  , predicted = predicted.class, raw = y_pred
                                  , stringsAsFactors = F)
rownames(actual.vs.predicted) <- rownames(motifs_test)

out_file <- paste0(target_celltype, "_by_", model_celltype,"_model_pred.txt")
print(paste("saving output file:", out_file))
write.table(x = actual.vs.predicted , file = out_file, quote = F)
