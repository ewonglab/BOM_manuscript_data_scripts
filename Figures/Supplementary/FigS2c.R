suppressMessages({
    library(xgboost)
    library("rsample")
    library(stringr)
    library(pROC)
    library(ggplot2)
    library(RColorBrewer)
    library(gridExtra)
    library(grid)
})

add.missing.vars_xgb <- function(xgb.model , data){
  missing.vars <- setdiff(xgb.model$feature_names, colnames(data))
  if(length(missing.vars)>0){
    for(i in 1:length(missing.vars)){
      data[,missing.vars[i]] <- 0}}
  return(data)}

xgb_with_nCREs <- function(counts.tab, n, sample_seed, esr_value = 100){
    
    ## sample n CREs from each class
    pos_seqs <- rownames(counts.tab[counts.tab$binary_celltype == 1, ])
    ned_ids <- rownames(counts.tab[counts.tab$binary_celltype == 0, ])
    
    message("Sample CREs")
    set.seed(sample_seed)
    pos_sample <- sample(x = pos_seqs, size = n, replace = FALSE)
    set.seed(sample_seed)
    neg_sample <- sample(x = ned_ids, size = n, replace = FALSE)
    
    counts.tab <- counts.tab[rownames(counts.tab) %in% c(pos_sample, neg_sample), ]
    message("nrow counts.tab")
    print(nrow(counts.tab))
    
    counts.tab$celltype <- NULL
    
    if("sequence_name" %in% colnames(counts.tab)){
        rownames(counts.tab) <- counts.tab$sequence_name
        counts.tab$sequence_name <- NULL
    }
    
    counts.tab$binary_celltype <- as.numeric(counts.tab$binary_celltype)

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

    motifs_val <- motifs_val[,colnames(motifs_train)]

    dtrain <- xgb.DMatrix(label = as.numeric(motifs_train$binary_celltype)
                          , data = as.matrix(motifs_train[, colnames(motifs_train)[colnames(motifs_train)!="binary_celltype"]]))
    dvalid <- xgb.DMatrix(label = as.numeric(motifs_val$binary_celltype)
                          , data = as.matrix(motifs_val[, colnames(motifs_val)[colnames(motifs_val)!="binary_celltype"]]))

    suppressMessages({
        set.seed(123) 
        xgb <- xgb.train(data = dtrain, nrounds = 10000, eta = 0.01, max_depth = 6, subsample = 0.5
                         , colsample_bytree = 0.5, objective = "binary:logistic"
                         , watchlist = list(train = dtrain, validation = dvalid)
                         , early_stopping_rounds = esr_value, nthread = 4, eval_metric = "error"
                         , maximize = FALSE, verbose = 0)
    })
                             
    motifs_test <- add.missing.vars_xgb(xgb, motifs_test)
    test_labels <- motifs_test$binary_celltype
    motifs_test <- motifs_test[,xgb$feature_names]  
    set.seed(123)
    y_pred <- predict(xgb, data.matrix(motifs_test), type = "response")
    predicted.class <- y_pred > 0.5 
    predicted.class <- gsub("TRUE", 1, predicted.class)
    predicted.class <- gsub("FALSE", 0, predicted.class)
    actual.vs.predicted <- data.frame(actual = test_labels
                                      , predicted = predicted.class, raw = y_pred
                                      , stringsAsFactors = F)
    rownames(actual.vs.predicted) <- rownames(motifs_test)
    return(list(pos_seqs = pos_sample, neg_seqs = neg_sample, sampling_seed = sample_seed
                , model = xgb, predictions = actual.vs.predicted))

}
                             
data_path <- "/data/mouse/e8.25"
celltype <- "endothelium"
endothelium <- read.table(file = file.path(data_path, paste0(celltype, "_500bp_vert5.0_qval0.5.txt"))
                          , header = TRUE, stringsAsFactors = FALSE, sep = '\t')
# get number of pos and neg seqs
npos <- table(endothelium$binary_celltype)["1"]
nneg <- table(endothelium$binary_celltype)["0"]
n_min <- min(npos, nneg)

n_values <- seq(from = 30, to = n_min, by = 50)

set.seed(123)
SEEDS <- sample.int(n = 1000000, size = length(n_values), replace = FALSE)
names(SEEDS) <- n_values

n_models <- lapply(n_values, function(n) xgb_with_nCREs(counts.tab = endothelium, n = n
                                                        , sample_seed = SEEDS[[as.character(n)]]))
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
                   
n_stats <- lapply(n_models, function(n) pred_stats(n$predictions))# add n by class
for(i in 1:length(n_stats)){
    n_stats[[i]]$n_by_class <- n_values[i]
    n_stats[[i]] <- n_stats[[i]][, c("n_by_class"
                                     , colnames(n_stats[[i]])[colnames(n_stats[[i]]) != "n_by_class"])]
}
                  
                  
                  
n_stats.df <- do.call("rbind", n_stats)suppressMessages(library("yardstick"))
                  
suppressMessages({
  library("cvAUC")
  library(pROC)
  library(ggplot2)
})

my_theme <-  theme(panel.border = element_blank(), panel.grid.major = element_blank()
                   , panel.grid.minor = element_blank()
                   , axis.line = element_line(colour = "black", linewidth = 1)
                   , legend.title = element_blank()
                   , legend.key=element_blank()
                   , legend.position = "bottom"
                   # , legend.key.size = unit(0.5, "lines")
                   , legend.key.width = unit(0.8,"cm")
                   , panel.background = element_blank()
                   , text = element_text(size=15)
                   , axis.text.x=element_text(colour="black")
                   , axis.text.y=element_text(colour="black")
                   , legend.text=element_text(size=6)
                   , axis.ticks = element_line(colour = "black"))

plot_roc <- function(pred_tab){
    roc_pred_tab <- roc(pred_tab$actual, pred_tab$raw, direction="<")
    rocs.list <- list(roc_pred_tab)
    names(rocs.list) <- c("pred_tab")
    # Plot ROC curves
    p <- ggroc(rocs.list, size = 1) + theme(panel.border = element_blank(), panel.grid.major = element_blank()
                                            , panel.grid.minor = element_blank()#legend.position="none" 
                                        , axis.line = element_line(colour = "black", size = 1)
                                        , legend.title = element_blank()
                                        , legend.key=element_blank()
                                        , legend.position = c(0.75, 0.25)
                                        , legend.key.width = unit(1.5,"cm")
                                        , panel.background = element_blank()
                                        , text = element_text(size=23)
                                        , axis.text.x=element_text(colour="black")
                                        , axis.text.y=element_text(colour="black")
                                        , legend.text=element_text(size=24)
                                        , axis.ticks = element_line(colour = "black")) +
  guides(linetype = guide_legend(override.aes = list(size = 3))) +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.8, color = "grey") + coord_equal() +
  scale_colour_manual(values=c("goldenrod1")
                      , aesthetics = c("colour", "fill")
                      , labels = c(args[2])) +
    labs(y= "Sensitivity", x = "Specificity")
    return(p)
}

celltype_rocs <- function(df){
    roc_pred_tab <- roc(df$actual, df$raw, direction = "<")
    return(roc_pred_tab)
}


n_rocs <- lapply(n_models, function(n) celltype_rocs(n$predictions))
names(n_rocs) <- n_values
roc_df <- purrr::map2_dfr(n_rocs, names(n_rocs), function(roc_obj, n_train) {
  tibble(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    n_train = as.numeric(n_train)
  )
})
                 
                 
ggplot(roc_df, aes(x = fpr, y = tpr, group = n_train, color = n_train)) +
  geom_line(size = 1.2) +
  scale_color_gradient(low = "grey80", high = "purple4", name = "Dataset Size by class") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  coord_equal() +
  labs(
    title = "ROC Curves Colored by Training Size",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  theme_classic()
