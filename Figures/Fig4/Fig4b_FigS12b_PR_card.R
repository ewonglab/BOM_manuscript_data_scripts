library("yardstick")
library("cvAUC")
library(pROC)
library(ggplot2)
library(reshape2)
library(pheatmap)
library("ggpubr")

setwd("/xgb/cross_sp")
fetal_path <- "/xgb/fetal/"
e8.25_path <- "/xgb/e8.25/bin_500bp/"
human <- read.table(file = paste0(fetal_path, "Cardiomyocytes_HEART_v_other_distal_nc_500bp_vert5.0_304T_pred.txt")
                    , header = T, stringsAsFactors = F)
mouse <- read.table(file = paste0(e8.25_path, "cardiom_500bp_vert5.0_qval0.5_578T_pred.txt")
                    , header = T, stringsAsFactors = F)
human_by_mouse <- read.table(file = "human_cm_vs_other_byMouseE8.25_vert5.0_q0.5_500bp.txt"
                             , header = T, stringsAsFactors = F)
mouse_by_human <- read.table(file = "mouse_cm_vs_other_byHumanFetal_vert5.0_q0.5_500bp.txt"
                             , header = T, stringsAsFactors = F)

get_pr_curve <- function(df){
  df$actual <- factor(df$actual, levels=c(0,1))  
  curve <- pr_curve(df, truth=actual, raw, event_level="second")
  return(as.data.frame(curve))
}

pr_human <- get_pr_curve(human)
pr_mouse <- get_pr_curve(mouse)
pr_human_by_mouse <- get_pr_curve(human_by_mouse)
pr_mouse_by_human <- get_pr_curve(mouse_by_human)

## add model inforation
pr_human$model <- "human"
pr_mouse$model <- "mouse"
pr_human_by_mouse$model <- "human_by_mouse"
pr_mouse_by_human$model <- "mouse_by_human"

mouse_curves <- rbind(pr_mouse, pr_mouse_by_human)
human_curves <- rbind(pr_human, pr_human_by_mouse)

range(mouse_curves$recall)
# [1] 0 1
range(mouse_curves$precision)
# [1] 0.5014426 1.0000000
range(human_curves$recall)
# [1] 0 1
range(human_curves$precision)
# [1] 0.4944681 1.0000000

mouse_curves$model <- as.factor(mouse_curves$model)#
human_curves$model <- as.factor(human_curves$model)#

my_theme <-  theme(panel.border = element_blank(), panel.grid.major = element_blank()
                   , panel.grid.minor = element_blank()#legend.position="none"
                   , axis.line = element_line(colour = "black", size = 1)
                   , legend.title = element_blank()
                   , legend.key=element_blank()
                   , legend.position = c(0.75, 0.25)
                   , legend.key.width = unit(1,"cm")
                   , panel.background = element_blank()
                   , legend.text=element_text(size=15)
                   , text = element_text(size=20)
                   , axis.text.x=element_text(colour="black")
                   , axis.text.y=element_text(colour="black")
                   , axis.ticks = element_line(colour = "black"))

p_mouse <- ggplot(mouse_curves, aes(x = recall, y = precision, color = model)) +
  geom_path() + coord_equal() +  my_theme +
  scale_colour_manual(values=c("#0B3954", "#FF6B6B")
                      , aesthetics = c("colour", "fill")) + ylim(0.45, 1)

p_human <- ggplot(human_curves, aes(x = recall, y = precision, color = model)) +
  geom_path() + coord_equal() +  my_theme +
  scale_colour_manual(values=c("#0B3954", "#FF6B6B")
                      , aesthetics = c("colour", "fill")) + ylim(0.45, 1)

pdf("cross_species_cardiom_PR_curves_500bp.pdf")
ggarrange(plotlist = list(p_human, p_mouse))
dev.off()
