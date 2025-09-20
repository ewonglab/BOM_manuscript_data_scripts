library(reshape2)
library(tidyr)
library(pheatmap)

setwd("/xgb/Atha/root")

# READ SHAP RANKING AND GET TOP 10 FOR EACH MODEL CLASS
shapRank <- read.table(file = "trimm600bp_4Celltypes_SHAPrank", header = T
                       , sep ='\t', stringsAsFactors = F)
top10 <- (shapRank[shapRank$rank %in% 1:10, "motif"])

shapRank <- shapRank[shapRank$motif %in% top10,]
shapRank <- tidyr::spread(shapRank[,c("motif", "cellType", "sub_abs_SHAP")]
                          , cellType, sub_abs_SHAP)
rownames(shapRank) <- shapRank$motif
shapRank$motif <- NULL

pdf("trimm600bp_4Celltypes_SHAP_top10.pdf")
pheatmap(shapRank, cluster_rows = T, cluster_cols = T, scale = "row"
         , color = rev(hcl.colors(50, "BluYl")))
#colorRampPalette(c("darkseagreen1", "#4B0082"))(100)
dev.off()
