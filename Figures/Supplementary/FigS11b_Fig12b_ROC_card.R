### CROSS SPECIES ROC CURVES - ERYTHROID/ERYTHROBLASTS
library("cvAUC")
library(pROC)
library(ggplot2)
library("ggpubr")

setwd("/xgb/cross_sp")
fetal_path <- "/xgb/fetal/"
e8.25_path <- "/xgb/e8.25/bin_500bp/"
human <- read.table(file = paste0(fetal_path, "Erythroblasts_MULTI_v_other_distal_nc_500bp_vert5.0_10T_pred.txt")
                    , header = T, stringsAsFactors = F)
mouse <- read.table(file = paste0(e8.25_path, "erythroid_500bp_vert5.0_qval0.5_20T_pred.txt")
                    , header = T, stringsAsFactors = F)
human_by_mouse <- read.table(file = "human_Erythroblasts_vs_other_byMouseE8.25Eryth_vert5.0_q0.5_500bp.txt"
                             , header = T, stringsAsFactors = F)
mouse_by_human <- read.table(file = "mouse_eryth_vs_other_byFetalErythroblasts_vert5.0_q0.5_500bp.txt"
                             , header = T, stringsAsFactors = F)

roc_human <- roc(human$actual, human$raw, direction="<")
# Setting levels: control = 0, case = 1
roc_mouse <- roc(mouse$actual, mouse$raw, direction="<")
# Setting levels: control = 0, case = 1
roc_human_by_mouse <- roc(human_by_mouse$actual, human_by_mouse$raw, direction="<")
# Setting levels: control = 0, case = 1
roc_mouse_by_human <- roc(mouse_by_human$actual, mouse_by_human$raw, direction="<")
# Setting levels: control = 0, case = 1

rocs.human <- list(roc_human, roc_human_by_mouse)
names(rocs.human) <- c("human", "human_by_mouse")

rocs.mouse <- list(roc_mouse, roc_mouse_by_human)
names(rocs.mouse) <- c("mouse", "mouse_by_human")

p_human <- ggroc(rocs.human, size = 1) + theme(panel.border = element_blank(), panel.grid.major = element_blank()
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
  scale_colour_manual(values=c("#005B5C", "#FF6978")
                      , aesthetics = c("colour", "fill")
                      , labels = c("human", "human_by_mouse")) +
  labs(y= "Sensitivity", x = "Specificity")

p_mouse <- ggroc(rocs.mouse, size = 1) + theme(panel.border = element_blank(), panel.grid.major = element_blank()
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
  scale_colour_manual(values=c("#0B3954", "#FF6B6B")
                      , aesthetics = c("colour", "fill")
                      , labels = c("mouse", "mouse_by_human")) +
  labs(y= "Sensitivity", x = "Specificity")

pdf("cross_species_erythr_rocs_500bp.pdf")
ggarrange(plotlist = list(p_human, p_mouse))
dev.off()
