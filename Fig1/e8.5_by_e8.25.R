library("rsample")
library(xgboost)
library(readr)
library(stringr)
library(caret)
library(pROC)
library(ggplot2)      # model visualization
library(RColorBrewer)
library(gridExtra)
library(grid)

add.missing.vars_xgb <- function(xgb.model , data){
  missing.vars <- setdiff(xgb.model$feature_names, colnames(data))
  if(length(missing.vars)>0){
    for(i in 1:length(missing.vars)){
      data[,missing.vars[i]] <- 0}}
  return(data)}

args <- commandArgs(trailingOnly=TRUE)

target_celltype <- args[1]
test_set <- args[2]
xgb <- readRDS(paste0(target_celltype, "_500bp_vert5.0_qval0.5.rds"))

counts_path <- "/data/mouse/e8.5/"
counts.tab <- read.table(file = paste0(counts_path, test_set)
                         , header =T, stringsAsFactors = F, sep = '\t')
counts.tab$celltype <- NULL

if("sequence_name" %in% colnames(counts.tab)){
  rownames(counts.tab) <- counts.tab$sequence_name
  counts.tab$sequence_name <- NULL
}

counts.tab.NAs <- sapply(counts.tab, function(x) sum(is.na(x)))
motifs_test <- counts.tab

print("Test dataset:")
print(dim(motifs_test))
test_labels <- motifs_test$binary_celltype
motifs_test <- add.missing.vars_xgb(xgb, motifs_test)
motifs_test <- motifs_test[,xgb$feature_names]  
set.seed(123)
y_pred <- predict(xgb, data.matrix(motifs_test), type="response")
predicted.class <- y_pred > 0.5 
predicted.class <- gsub("TRUE", 1, predicted.class)
predicted.class <- gsub("FALSE", 0, predicted.class)
actual.vs.predicted <- data.frame(actual=test_labels, predicted = predicted.class, raw = y_pred
                                  , stringsAsFactors = F)
rownames(actual.vs.predicted) <- rownames(motifs_test)

write.table(x = actual.vs.predicted
            , file = paste0(gsub(".txt", "_", basename(test_set))
                            , xgb$best_iteration,"T_pred.txt"), quote = F)
