suppressMessages({
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

wd <- args[1]
input_data <- args[2]
xgb_model <- args[3]

add.missing.vars_xgb <- function(xgb.model , data){
  missing.vars <- setdiff(xgb.model$feature_names, colnames(data))
  if(length(missing.vars)>0){
    for(i in 1:length(missing.vars)){
      data[,missing.vars[i]] <- 0}}
  return(data)}

setwd(wd)

xgb <- readRDS(xgb_model)
xgb.save(xgb, gsub(".rds", ".bin", xgb_model))

counts.tab <- read.table(file = input_data, header =T, stringsAsFactors = F, sep = '\t')
counts.tab$celltype <- NULL

if("sequence_name" %in% colnames(counts.tab)){
  rownames(counts.tab) <- counts.tab$sequence_name
  counts.tab$sequence_name <- NULL
}

counts.tab.NAs <- sapply(counts.tab, function(x) sum(is.na(x)))
print(unique(counts.tab.NAs))
counts.tab$binary_celltype <-as.numeric(counts.tab$binary_celltype)

set.seed(123)
motifs_split <- initial_split(counts.tab, prop = .6)
motifs_train <- training(motifs_split)
motifs_test <- testing(motifs_split)

set.seed(123)
motifs_split2 <- initial_split(motifs_test, prop = .5)
motifs_val <- training(motifs_split2)
motifs_test <- testing(motifs_split2)

print("Training dataset:")
print(dim(motifs_train))
print("Validation dataset:")
print(dim(motifs_val))
print("Test dataset:")
print(dim(motifs_test))

motifs_train.sd <- apply(motifs_train, 2, sd)
motifs_train <- motifs_train[, names(which(motifs_train.sd != 0))]


write.csv(x = motifs_train
          , file = paste0(gsub(".txt", "_", basename(input_data)), xgb$best_iteration, "T_trainSet.csv")
          , quote = F)motifs_test <- add.missing.vars_xgb(xgb, motifs_test)
test_labels <- motifs_test$binary_celltype
motifs_test <- motifs_test[,xgb$feature_names]  
                         
set.seed(123)
y_pred <- predict(xgb, data.matrix(motifs_test), type="response")
predicted.class <- y_pred > 0.5 
predicted.class <- gsub("TRUE", 1, predicted.class)
predicted.class <- gsub("FALSE", 0, predicted.class)
actual.vs.predicted <- data.frame(actual = test_labels
                                  , predicted = predicted.class, raw = y_pred
                                  , stringsAsFactors = F)
rownames(actual.vs.predicted) <- rownames(motifs_test)
write.table(x = actual.vs.predicted
            , file = paste0(gsub(".txt", "_", basename(input_data))
                            , xgb$best_iteration,"T_pred.txt"), quote = F)
