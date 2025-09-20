suppressMessages({
  library("cvAUC")
  library(pROC)
  library(ggplot2)
})

args <- commandArgs(trailingOnly=TRUE)
celltype <- args[1]
pred_files <- list.files(path = "."
                         , pattern = paste0(celltype, "_nc_vs_other_vert5.0_qval(.*)_pred.txt"))
print("Reading predictions")
pred_q0.1 <- read.table(file = grep("qval0.1", pred_files, value =T)
                        , header = T, stringsAsFactors = F)
pred_q0.3 <- read.table(file = grep("qval0.3", pred_files, value =T)
                        , header = T, stringsAsFactors = F)
pred_q0.5 <- read.table(file = grep("qval0.5", pred_files, value =T)
                        , header = T, stringsAsFactors = F)

roc_q0.1 <- roc(pred_q0.1$actual, pred_q0.1$raw, direction="<")
roc_q0.3 <- roc(pred_q0.3$actual, pred_q0.3$raw, direction="<")
roc_q0.5 <- roc(pred_q0.5$actual, pred_q0.5$raw, direction="<")

rocs.list <- list(roc_q0.1, roc_q0.3, roc_q0.5)
names(rocs.list) <- c("qval<=0.1", "qval<=0.3", "qval<=0.5")
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
  scale_colour_manual(values=c("#E69F00", "#0072B2", "#009E73")
                      , aesthetics = c("colour", "fill")
                      , labels = names(rocs.list)) +
  labs(y= "Sensitivity", x = "Specificity") + ggtitle(paste(celltype, "vs other"))

pdf(paste0(celltype, "_qval_cutoffs_rocs.pdf"))
p
dev.off()
