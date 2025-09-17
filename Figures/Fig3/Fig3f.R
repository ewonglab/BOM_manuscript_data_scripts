
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

bom <- read.table(file = "atha_CisBP2.0_4CellT_600bp_q0.5_1137T_pred.txt"
                  , header = T, stringsAsFactors = F) # 122   7

# order columns
bom <- bom[,c(paste0("V", 1:4), "label", "max_prob", "celltype")]

# firt column is NOT topic1 etc
celltype_annot <- unique(bom[,c("label", "celltype")])
celltype_annot <- celltype_annot[with(celltype_annot, order(label)),]
celltype_annot
#                     label      celltype
# 1:12844771-12845371     0 c.e_precursor
# 1:11620972-11621572     1  endodermis_2
# 4:5607263-5607863       2  endodermis_3
# 1:13361345-13361945     3 stele_1.xylem

# change column names
n_celltypes <- nrow(celltype_annot)
colnames(bom)[1:n_celltypes] <- celltype_annot$celltype
# change label names
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

## setting line colours
# my_colors <- c("#ffef96", "#50394c", "#b2b2b2", "#f4e1d2")
my_colors <- c("#c9cba3", "#ffe1a8", "#e26d5c", "#723d46")

# Plot ROC curves
pdf("Ath600bp_mult4CT_CisBP_pred_ROCs.pdf")
ggroc(bom_rocs) + my_theme +
  guides(linetype = guide_legend(override.aes = list(size = 3))) +
  geom_abline(slope=1, intercept = 1, linetype = "solid", alpha=0.8
              , color = "grey") + coord_equal() +
  scale_colour_manual(values = my_colors, aesthetics = c("colour", "fill")
                      , labels = c(names(bom_rocs))) +
  labs(y= "Sensitivity", x = "Specificity")
dev.off()

ROC_auc <- sapply(1:length(bom_rocs), function(x) bom_rocs[[x]]$auc)
ROC_auc <- as.data.frame(ROC_auc)
ROC_auc$CellType <- names(bom_rocs)
ROC_auc$ROC_auc <- round(ROC_auc$ROC_auc, 3)
                
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
bom <- read.table(file = "atha_CisBP2.0_4CellT_600bp_q0.5_1137T_pred.txt"
                  , header = T, stringsAsFactors = F) 

# order columns
bom <- bom[,c(paste0("V", 1:4), "label", "max_prob", "celltype")]

celltype_annot <- unique(bom[,c("label", "celltype")])
celltype_annot <- celltype_annot[with(celltype_annot, order(label)),]
celltype_annot
#                     label      celltype
# 1:12844771-12845371     0 c.e_precursor
# 1:11620972-11621572     1  endodermis_2
# 4:5607263-5607863       2  endodermis_3
# 1:13361345-13361945     3 stele_1.xylem

# change column names
n_celltypes <- nrow(celltype_annot)
colnames(bom)[1:n_celltypes] <- celltype_annot$celltype
# change label names
bom$label <- bom$celltype

bom_prcurves <- (lapply(colnames(bom)[1:n_celltypes], celltype_PR_curve))
bom_prcurves <- lapply(bom_prcurves, as.data.frame)
names(bom_prcurves) <- colnames(bom)[1:n_celltypes]
for(i in 1:length(bom_prcurves)){
  bom_prcurves[[i]]$CellType <- names(bom_prcurves)[i]
}

bom_prcurves <- do.call("rbind", bom_prcurves)

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

## setting line colours
my_colors <- c("#c9cba3", "#ffe1a8", "#e26d5c", "#723d46")


# Plot PR curves

pdf("Ath600bp_mult4CT_CisBP_pred_PRcurves.pdf")
ggplot(bom_prcurves, aes(x = recall, y = precision, color = CellType)) +
  geom_path() + coord_equal() +  my_theme +
  scale_color_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_abline(slope = -1, intercept = 1, linetype = "solid", alpha = 0.8
              , color = "grey")
dev.off()
