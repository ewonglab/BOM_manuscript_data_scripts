suppressMessages({
  library("yardstick")
  library("cvAUC")
  library(pROC)
  library(ggplot2)
})

setwd("/xgb/e8.25/non_overlapping")

pred_files <- c("allantois_vs_other_vert5.0_qval0.5_nonOv_38T_pred.txt"
                , "mixed_meso_vs_other_vert5.0_qval0.5_nonOv_21T_pred.txt"
                , "cardiom_vs_other_vert5.0_qval0.5_nonOv_61T_pred.txt"
                ,"neuralcrest_vs_other_vert5.0_qval0.5_nonOv_72T_pred.txt"
                , "endothelium_vs_other_vert5.0_qval0.5_nonOv_176T_pred.txt"
                , "NMP_vs_other_vert5.0_qval0.5_nonOv_213T_pred.txt"
                , "erythroid_vs_other_vert5.0_qval0.5_nonOv_12T_pred.txt"
                , "paraxial_meso_vs_other_vert5.0_qval0.5_nonOv_90T_pred.txt"
                , "exe_endo_vs_other_vert5.0_qval0.5_nonOv_18T_pred.txt"
                , "pharyngeal_meso_vs_other_vert5.0_qval0.5_nonOv_103T_pred.txt"
                , "forebrain_vs_other_vert5.0_qval0.5_nonOv_3T_pred.txt"
                , "somitic_meso_vs_other_vert5.0_qval0.5_nonOv_11T_pred.txt"
                , "gut_vs_other_vert5.0_qval0.5_nonOv_132T_pred.txt"
                , "spinalcord_vs_other_vert5.0_qval0.5_nonOv_145T_pred.txt"
                , "mesenchyme_vs_other_vert5.0_qval0.5_nonOv_168T_pred.txt"
                , "surface_ecto_vs_other_vert5.0_qval0.5_nonOv_11T_pred.txt"
                , 'mid_hindbrain_vs_other_vert5.0_qval0.5_nonOv_38T_pred.txt')

pred_li <- lapply(pred_files, read.table, header = T, stringsAsFactors = F)

celltypes <- sub("_vs_other_.*", "", pred_files)
names(pred_li) <- celltypes
bom_rocs <- lapply(pred_li, function(x) roc(x$actual, x$raw, direction="<"))

names(bom_rocs) <- celltypes
my_colours <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442"
                , "#56B4E9", "#E69F00", "#009292", "#E41A1C", "#F781BF", "#66C2A5"
                , "#FC8D62", "#A6D854", "#5A4A7F", "#F4A460", "#4682B4", "#E6B8B7")

colour.df <- data.frame(celltype=celltypes, celltype_colour = my_colours)

my_theme <-  theme(panel.border = element_blank(), panel.grid.major = element_blank()
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
                   , legend.text=element_text(size=10)
                   , axis.ticks = element_line(colour = "black"))

rownames(colour.df) <- colour.df$celltype
colour.df <- colour.df[names(bom_rocs), ]

p <- ggroc(bom_rocs, size = 1) + my_theme +
  guides(linetype = guide_legend(override.aes = list(size = 3))) +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.8, color = "grey") +
  coord_equal() +
  scale_color_manual(values=unique(colour.df$celltype_colour)) +
  labs(y= "Sensitivity", x = "Specificity")


pdf("mouseE8.25_nc_nonOverlap_BOM_binary_pred_ROCs.pdf")
p
dev.off()

pred_li <- lapply(pred_files, read.table, header = T, stringsAsFactors = F)

celltypes <- sub("_vs_other_.*", "", pred_files)
names(pred_li) <- celltypes

# actual class as factor
for(i in 1:length(pred_li)){
  pred_li[[i]]$actual <- factor(pred_li[[i]]$actual, levels=c(0,1))
}

pr_curves <- lapply(pred_li, pr_curve, truth=actual, raw,event_level="second")

plot_pre_curve <- function(x){
  p <- ggplot(x, aes(x = recall, y = precision)) +
    geom_path() + coord_equal() +  theme_classic() #+ ggtitle(my_title)
  return(p)}

pr_curves <- lapply(pr_curves, as.data.frame)
for(i in 1:length(pr_curves)){
  pr_curves[[i]]$xgb_model <- names(pr_curves)[i]
}

pr_curves <- do.call("rbind", pr_curves)#59449     4

my_theme <-  theme(panel.border = element_blank(), panel.grid.major = element_blank()
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
                   , legend.text=element_text(size=10)
                   , axis.ticks = element_line(colour = "black"))

my_colours <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442"
                , "#56B4E9", "#E69F00", "#009292", "#E41A1C", "#F781BF", "#66C2A5"
                , "#FC8D62", "#A6D854", "#5A4A7F", "#F4A460", "#4682B4", "#E6B8B7")

colour.df <- data.frame(celltype=celltypes, celltype_colour = my_colours)

pr_curves$ord <- 1:nrow(pr_curves)
pr_curves <- merge(pr_curves, colour.df, by.x="xgb_model", by.y="celltype", all.x=T)
pr_curves <- pr_curves[with(pr_curves, order(ord)),]
pr_curves$ord <- NULL

pr_curves$xgb_model <- factor(pr_curves$xgb_model
                              , levels = c("allantois", "mixed_meso", "cardiom"
                                           , "neuralcrest", "endothelium", "NMP"
                                           , "erythroid", "paraxial_meso"
                                           , "exe_endo", "pharyngeal_meso"
                                           , "forebrain", "somitic_meso"
                                           , "gut", "spinalcord", "mesenchyme"
                                           , "surface_ecto", "mid_hindbrain"))

p <- ggplot(pr_curves, aes(x = recall, y = precision, color = xgb_model)) +
  geom_path() + coord_equal() +  my_theme +
  scale_color_manual(values=unique(pr_curves$celltype_colour)) +
  scale_y_continuous(limits = c(0, 1))

pdf("E8.25_nc_nonOverlap_binary_PRcurves.pdf")
p
dev.off()


