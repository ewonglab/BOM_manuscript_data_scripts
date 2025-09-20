suppressMessages({
  library("yardstick")
  library("cvAUC")
  library(pROC)
  library(ggplot2)
})

setwd("xgb/e8.25/qvals")
pred_files <- c("allantois_vs_other_vert5.0_qval0.3_49T_pred.txt"
                , "mixed_meso_vs_other_vert5.0_qval0.3_303T_pred.txt"
                , "cardiom_vs_other_vert5.0_qval0.3_582T_pred.txt"
                , "neuralcrest_vs_other_vert5.0_qval0.3_7T_pred.txt"
                , "endothelium_vs_other_vert5.0_qval0.3_79T_pred.txt"
                , "NMP_vs_other_vert5.0_qval0.3_336T_pred.txt"
                , "erythroid_vs_other_vert5.0_qval0.3_245T_pred.txt"
                , "paraxial_meso_vs_other_vert5.0_qval0.3_13T_pred.txt"
                , "exe_endo_vs_other_vert5.0_qval0.3_5T_pred.txt"
                , "pharyngeal_meso_vs_other_vert5.0_qval0.3_34T_pred.txt"
                , "forebrain_vs_other_vert5.0_qval0.3_269T_pred.txt"
                , "somitic_meso_vs_other_vert5.0_qval0.3_114T_pred.txt"
                , "gut_vs_other_vert5.0_qval0.3_389T_pred.txt"
                , "spinalcord_vs_other_vert5.0_qval0.3_355T_pred.txt"
                , "mesenchyme_vs_other_vert5.0_qval0.3_109T_pred.txt"
                , "surface_ecto_vs_other_vert5.0_qval0.3_5T_pred.txt"
                , 'mid_hindbrain_vs_other_vert5.0_qval0.3_179T_pred.txt')

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
  # scale_colour_manual(values=roc_colours, aesthetics = c("colour", "fill")
  #                     , labels = c(names(bom_rocs))) +
  scale_color_manual(values=unique(colour.df$celltype_colour)) +
  labs(y= "Sensitivity", x = "Specificity")

pdf("mouseE8.25_nc_BOM_binary_qval0.3_pred_ROCs.pdf")
p
dev.off()

pred_files <- c("allantois_vs_other_vert5.0_qval0.1_3T_pred.txt"
                , "mixed_meso_vs_other_vert5.0_qval0.1_108T_pred.txt"
                , "cardiom_vs_other_vert5.0_qval0.1_51T_pred.txt"
                , "neuralcrest_vs_other_vert5.0_qval0.1_8T_pred.txt"
                , "endothelium_vs_other_vert5.0_qval0.1_14T_pred.txt"
                , "NMP_vs_other_vert5.0_qval0.1_107T_pred.txt"
                , "erythroid_vs_other_vert5.0_qval0.1_20T_pred.txt"
                , "paraxial_meso_vs_other_vert5.0_qval0.1_5T_pred.txt"
                , "exe_endo_vs_other_vert5.0_qval0.1_2T_pred.txt"
                , "pharyngeal_meso_vs_other_vert5.0_qval0.1_94T_pred.txt"
                , "forebrain_vs_other_vert5.0_qval0.1_65T_pred.txt"
                , "somitic_meso_vs_other_vert5.0_qval0.1_52T_pred.txt"
                , "gut_vs_other_vert5.0_qval0.1_13T_pred.txt"
                , "spinalcord_vs_other_vert5.0_qval0.1_18T_pred.txt"
                , "mesenchyme_vs_other_vert5.0_qval0.1_13T_pred.txt"
                , "surface_ecto_vs_other_vert5.0_qval0.1_25T_pred.txt"
                , 'mid_hindbrain_vs_other_vert5.0_qval0.1_271T_pred.txt')
         
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
  # scale_colour_manual(values=roc_colours, aesthetics = c("colour", "fill")
  #                     , labels = c(names(bom_rocs))) +
  scale_color_manual(values=unique(colour.df$celltype_colour)) +
  labs(y= "Sensitivity", x = "Specificity")

pdf("mouseE8.25_nc_BOM_binary_qval0.1_pred_ROCs.pdf")
p
dev.off()
