library("yardstick")
library("cvAUC")
library(pROC)
library(ggplot2)

celltype_roc_curve <- function(celltype){
  #make a binary variable to represent actual celltype
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
  # print(head(df))
  colnames(df) <- c("raw", "label", "actual")
  ct_roc <- roc(df$actual, df$raw, direction="<")
  return(ct_roc)
}

setwd("/bom/YSed")

bom <- read.table(file = "human_by_ysed_matchedCT_ysed.posMarkerP_padj0.05_pred.txt"
                  , header = T, stringsAsFactors = F) #
celltype_annot <- unique(bom[,c("label", "celltype")])
celltype_annot <- celltype_annot[with(celltype_annot, order(label)),]

# change column names
n_celltypes <- 7
# replace label as the actual topic name
bom$label <- bom$celltype
bom_rocs <- (lapply(colnames(bom)[1:n_celltypes], celltype_roc_curve))
names(bom_rocs) <- colnames(bom)[1:n_celltypes]

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

my_colors <- rainbow(n_celltypes)

# Plot ROC curves
pdf("human_by_ysed_matchedCT_ysed.posMarkerP_padj0.05_ROCs.pdf")
ggroc(bom_rocs) + my_theme +
  guides(linetype = guide_legend(override.aes = list(size = 3))) +
  geom_abline(slope=1, intercept = 1, linetype = "solid", alpha=0.8
              , color = "grey") + coord_equal() +
  scale_colour_manual(values=my_colors, aesthetics = c("colour", "fill")
                      , labels = c(names(bom_rocs))) +
  labs(y = "Sensitivity", x = "Specificity")
dev.off()
