library(ggplot2)
library("ggpubr")
library("mclust")
library("data.table")

rank_motifs <- function(shap.df, train.df){
  # shap only for the target class
  # shap.df <- shap.df[which(train.df$celltype_numeric == model_class),]
  abs.sum.shap_per_motif <- apply(shap.df, 2, function(x) sum(abs(x)))
  names(abs.sum.shap_per_motif) <- colnames(train.df)[-ncol(train.df)]
  abs.sum.shap_per_motif <- abs.sum.shap_per_motif[order(-abs.sum.shap_per_motif)]
  ranked_motifs <- (abs.sum.shap_per_motif)  
  ranked_motifs <- as.data.frame(ranked_motifs)
  colnames(ranked_motifs)[1] <- "sub_abs_SHAP"
  ranked_motifs$motif <- rownames(ranked_motifs)
  rownames(ranked_motifs) <- 1:nrow(ranked_motifs)
  ranked_motifs$rank <- 1:nrow(ranked_motifs)
  return(ranked_motifs[,c("motif", "sub_abs_SHAP", "rank")])}

shap_f <- paste0("trimm600bp_4Celltypes_miltiXGB_SHAP_", 1:4, ".txt")
train_f <- "atha_CisBP2.0_4CellT_600bp_q0.5_1137T_trainSet.csv"

#READING SHAP VALUES
shap.df <- lapply(shap_f, fread)
shap.df <- lapply(shap.df, as.data.frame)
train_counts.df <- read.csv(file = train_f, header = T, stringsAsFactors = F, row.names = 1)

class1_ranked <- rank_motifs(shap.df[[1]], train_counts.df)
class2_ranked <- rank_motifs(shap.df[[2]], train_counts.df)
class3_ranked <- rank_motifs(shap.df[[3]], train_counts.df)
class4_ranked <- rank_motifs(shap.df[[4]], train_counts.df)

# add model class
class1_ranked$cellType <- "c.e_precursor"
class2_ranked$cellType <- "endodermis_2"
class3_ranked$cellType <- "endodermis_3"
class4_ranked$cellType <- "stele_1.xylem"

# combine in single data frame
allClass <- rbind(class1_ranked, class2_ranked, class3_ranked, class4_ranked)

# save:
write.table(x = allClass, file = "trimm600bp_4Celltypes_SHAPrank", quote = F, sep ='\t')
