library("rsample")
library(xgboost)# version 1.6.0.1
library(readr)
library(stringr)
library(pROC)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)

args <- commandArgs(trailingOnly=TRUE)

input_data <- args[1]
out_file <- args[2]
esr_value <- as.integer(args[3])

counts.tab <- read.table(file = input_data, header = T, stringsAsFactors = F, sep = '\t')
counts.tab$celltype <- NULL

if("sequence_name" %in% colnames(counts.tab)){
  rownames(counts.tab) <- counts.tab$sequence_name
  counts.tab$sequence_name <- NULL
}

counts.tab.NAs <- sapply(counts.tab, function(x) sum(is.na(x)))
print(unique(counts.tab.NAs))#0
counts.tab$binary_celltype <-as.numeric(counts.tab$binary_celltype)

set.seed(123)
motifs_split <- initial_split(counts.tab, prop = .6)
motifs_train <- training(motifs_split)
motifs_test <- testing(motifs_split)

set.seed(123)
motifs_split2 <- initial_split(motifs_test, prop = .5)
motifs_val <- training(motifs_split2)
motifs_test <- testing(motifs_split2)

motifs_train.sd <- apply(motifs_train, 2, sd)
motifs_train <- motifs_train[, names(which(motifs_train.sd != 0))]
print("Training dataset after filter:")
print(dim(motifs_train))

motifs_val <- motifs_val[,colnames(motifs_train)]

dtrain <- xgb.DMatrix(label = as.numeric(motifs_train$binary_celltype)
                      , data = as.matrix(motifs_train[, colnames(motifs_train)[colnames(motifs_train)!="binary_celltype"]]))
dvalid <- xgb.DMatrix(label = as.numeric(motifs_val$binary_celltype)
                      , data = as.matrix(motifs_val[, colnames(motifs_val)[colnames(motifs_val)!="binary_celltype"]]))


set.seed(123) 
xgb <- xgb.train(data = dtrain, nrounds = 10000, eta = 0.01, max_depth = 6, subsample = 0.5
                 , colsample_bytree = 0.5, objective = "binary:logistic"
                 , watchlist = list(train = dtrain, validation = dvalid)
                 , early_stopping_rounds = esr_value, nthread = 4, eval_metric = "error"
                 , maximize = F)
saveRDS(object = xgb, file = out_file)
