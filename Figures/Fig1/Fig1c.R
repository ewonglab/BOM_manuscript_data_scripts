suppressMessages({
  library("yardstick")
  library("cvAUC")
  library(pROC)
  library(ggplot2)
})

setwd("/xgb/e8.25/bin_500bp")
pred_files <- list.files(path = ".", pattern = "(.*)_500bp_vert5.0_qval0.5_(.*)T_pred.txt")

pred_li <- lapply(pred_files, read.table, header = T, stringsAsFactors = F)
celltypes <- sub("_500bp.*", "", pred_files)
names(pred_li) <- celltypes
bom_rocs <- lapply(pred_li, function(x) roc(x$actual, x$raw, direction="<"))

names(bom_rocs) <- celltypes

my_colours <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442"
                , "#56B4E9", "#E69F00", "#009292", "#E41A1C", "#F781BF", "#66C2A5"
                , "#FC8D62", "#A6D854", "#5A4A7F", "#F4A460", "#4682B4", "#E6B8B7")


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

# match colors to other figures

celltype_order <- c("allantois", "mixed_meso", "cardiom", "neuralcrest", "endothelium", "NMP"
                    , "erythroid", "paraxial_meso", "exe_endo", "pharyngeal_meso"
                    , "forebrain", "somitic_meso", "gut", "spinalcord", "mesenchyme"
                    , "surface_ecto", "mid_hindbrain")

bom_rocs <- bom_rocs[celltype_order]


p <- ggroc(bom_rocs, size = 1) + my_theme +
  guides(linetype = guide_legend(override.aes = list(size = 3))) +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.8, color = "grey") +
  coord_equal() +
  # scale_colour_manual(values=roc_colours, aesthetics = c("colour", "fill")
  #                     , labels = c(names(bom_rocs))) +
  scale_color_manual(values=unique(colour.df$celltype_colour)) +
  labs(y= "Sensitivity", x = "Specificity")

pdf("E8.25_BOM_binary_500bp_pred_ROCs.pdf")
p
dev.off()
