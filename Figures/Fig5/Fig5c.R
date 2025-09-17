library(data.table)
pred_f <- list.files(path = ".", pattern = "cardiom_selected_near_Nkx2.5_GIMME_by_(.*).txt$")
pred <- lapply(pred_f, read.table, header = T, stringsAsFactors = F)

for(i in 1:length(pred)){
  pred[[i]]$model <- sub(".txt", "", sub(".*_by_", "", pred_f[i]))
  pred[[i]]$peakid <- rownames(pred[[i]])
  rownames(pred[[i]]) <- NULL
}

pred.df <- do.call("rbind", pred)
pred.df$actual <- NULL
pred.df$predicted <- NULL

library(tidyr)
pred.df <- tidyr::spread(pred.df, model, raw)
rownames(pred.df) <- pred.df$peakid
pred.df$peakid <- NULL


library(pheatmap)
library(RColorBrewer)
cols <- brewer.pal(10, "RdBu")
myColor <- colorRampPalette(rev(cols))(100)

myColor3 <- colorRampPalette(c("#f2edee","#D46228"))(100)
pdf("cardiom_selected_near_Nkx2.5_GIMME_prob.2.pdf")
pheatmap(mat = pred.df, color=(myColor3), cluster_cols = F, cluster_rows = F, scale = "none")
dev.off()

