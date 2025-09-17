suppressMessages({
  library("yardstick")
  library("cvAUC")
  library(pROC)
  library(ggplot2)
})

celltype_roc_curve <- function(celltype){
  #make a binary variable to represent actual celltype
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
  # print(head(df))
  colnames(df) <- c("raw", "label", "actual")
  ct_roc <- roc(df$actual, df$raw, direction="<")
  return(ct_roc)
}

bom <- read.table(file = "topics_vert5.0_q0.5_9999T_pred.txt"
                  , header =T, stringsAsFactors = F)
bom <- bom[,c(paste0("V", 1:93), "label", "max_prob", "celltype")]
topic_annot <- unique(bom[,c("label", "celltype")])
topic_annot <- topic_annot[with(topic_annot, order(label)),]

n_topics <- 93
colnames(bom)[1:n_topics] <- topic_annot$celltype
bom$label <- bom$celltype
bom_rocs <- (lapply(colnames(bom)[1:n_topics], celltype_roc_curve))
names(bom_rocs) <- colnames(bom)[1:n_topics]

my_theme <-  theme(panel.border = element_blank(), panel.grid.major = element_blank()
                   , panel.grid.minor = element_blank()
                   , axis.line = element_line(colour = "black", size = 1)
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

my_colors <- rainbow(n_topics)

p <- ggroc(bom_rocs) + my_theme +
  guides(linetype = guide_legend(override.aes = list(size = 3))) +
  geom_abline(slope=1, intercept = 1, linetype = "solid", alpha=0.8
              , color = "grey") + coord_equal() +
  scale_colour_manual(values=my_colors, aesthetics = c("colour", "fill")
                      , labels = c(names(bom_rocs))) +
  labs(y= "Sensitivity", x = "Specificity")

pdf("topics_multiclass_pred_ROCs_mlogloss.pdf")
p
dev.off()
