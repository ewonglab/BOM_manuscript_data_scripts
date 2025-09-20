suppressMessages({
  library("yardstick")
  library("cvAUC")
  library(pROC)
  library(ggplot2)
})

plot_pre_curve <- function(x){
  p <- ggplot(x, aes(x = recall, y = precision)) +
    geom_path() + coord_equal() +  theme_classic() #+ ggtitle(my_title)
  return(p)}

celltype_PR_curve <- function(celltype){
  #make a binary variable to represent actual celltype
  df <- bom[, c(celltype, "label")]
  df$actual <- ifelse(df$label == celltype, 1, 0)
  # print(head(df))
  colnames(df) <- c("raw", "label", "actual")
  df$actual <- factor(df$actual, levels=c(0,1))
  bop_PRcurve <- pr_curve(df, truth=actual, raw,event_level="second")
  return(bop_PRcurve)
}

bom <- read.table(file = "topics_vert5.0_q0.5_9999T_pred.txt"
                  , header =T, stringsAsFactors = F)
# order columns
bom <- bom[,c(paste0("V", 1:93), "label", "max_prob", "celltype")]
# firt column is NOT topic1 etc
topic_annot <- unique(bom[,c("label", "celltype")])
topic_annot <- topic_annot[with(topic_annot, order(label)),]

# change column names
n_topics <- 93
colnames(bom)[1:n_topics] <- topic_annot$celltype
# replace label as the actual topic name
bom$label <- bom$celltype
bom_prcurves <- (lapply(colnames(bom)[1:n_topics], celltype_PR_curve))

bom_prcurves <- lapply(bom_prcurves, as.data.frame)
names(bom_prcurves) <- colnames(bom)[1:n_topics]
for(i in 1:length(bom_prcurves)){
  bom_prcurves[[i]]$Topic <- names(bom_prcurves)[i]
}

bom_prcurves <- do.call("rbind", bom_prcurves)#[1] 1424832       4

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

## setting line colours, line types and line transparency
my_colors <- rainbow(n_topics)


# Plot PR curves
p <- ggplot(bom_prcurves, aes(x = recall, y = precision, color = Topic)) +
  geom_path() + coord_equal() +  my_theme +
  scale_color_manual(values=my_colors) +
  scale_y_continuous(limits = c(0, 1))


pdf("topics_multiclass_pred_PRcurves_mlogloss.pdf")
p
dev.off()
