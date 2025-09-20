library(data.table)
library(ggplot2)
library(reshape2)
library("ggpubr")
setwd("/xgb/cross_sp")

cm_human_by_mouse <- read.table(file = "human_cm_vs_other_byMouseE8.25_vert5.0_q0.5_500bp.txt"
                                , header = T, stringsAsFactors = F)
cm_mouse_by_human <- read.table(file = "mouse_cm_vs_other_byHumanFetal_vert5.0_q0.5_500bp.txt"
                                , header = T, stringsAsFactors = F)
ery_human_by_mouse <- read.table(file = "human_Erythroblasts_vs_other_byMouseE8.25Eryth_vert5.0_q0.5_500bp.txt"
                                 , header = T, stringsAsFactors = F)
ery_mouse_by_human <- read.table(file = "mouse_eryth_vs_other_byFetalErythroblasts_vert5.0_q0.5_500bp.txt"
                                 , header = T, stringsAsFactors = F)

get_conf_mat <- function(x){
  conf <- table(x[,c("actual", "predicted")])
  return(conf)
}

cm_human_by_mouse.mat <- get_conf_mat(cm_human_by_mouse)
cm_mouse_by_human.mat <- get_conf_mat(cm_mouse_by_human)
ery_human_by_mouse.mat <- get_conf_mat(ery_human_by_mouse)
ery_mouse_by_human.mat <- get_conf_mat(ery_mouse_by_human)

plot_conf_mat <- function(df, actual_label){
  df <- as.data.frame(df)
  df$actual <- ifelse(df$actual ==1, actual_label, "other")
  df$predicted <- ifelse(df$predicted ==1, actual_label, "other")
  # df <- melt(df)
  p <- ggplot(data =  df, mapping = aes(x = predicted, y = actual)) +
    geom_tile(aes(fill = Freq), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", Freq)), colour = "white"
              , vjust = 1, size=2) +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_bw() + theme(legend.position = "bottom"
                       , axis.text.x = element_text(size = 5)
                       , axis.text.y = element_text(size = 5))  +
    ggtitle(actual_label)
  return(p)
}

cm_human_by_mouse.p <- plot_conf_mat(cm_human_by_mouse.mat, "human_cardiom_by_mouse")
cm_mouse_by_human.p <- plot_conf_mat(cm_mouse_by_human.mat, "mouse_cardiom_by_human")
ery_human_by_mouse.p <- plot_conf_mat(ery_human_by_mouse.mat, "human_eryth_by_mouse")
ery_mouse_by_human.p <- plot_conf_mat(ery_mouse_by_human.mat, "mouse_eryth_by_human")

conf_mat_plots <- list(cm_human_by_mouse.p, cm_mouse_by_human.p
                       , ery_human_by_mouse.p, ery_mouse_by_human.p)

pdf("fetal_and_e8.25_crossSpec_confMat_500bp.pdf", onefile = F)
ggarrange(plotlist=conf_mat_plots, ncol=2, nrow=2, common.legend = T)
dev.off()
