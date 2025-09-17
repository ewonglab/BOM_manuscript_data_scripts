# ROC and PR curves - e8.5 predicted by e8.25 model
suppressMessages({
  library("yardstick")
  library("cvAUC")
  library(pROC)
  library(ggplot2)
  })


pred_files <- list.files(path = ".", pattern = "(.*)_vs_other_vert5.0_qval0.5_(.*)_pred.txt$")

pred_li <- lapply(pred_files, read.table, header = T, stringsAsFactors = F)

celltypes <- sub("_vs_.*", "", pred_files)

names(pred_li) <- celltypes
bom_rocs <- lapply(pred_li, function(x) roc(x$actual, x$raw, direction="<"))

colour.df <- data.frame(celltype = c("Allantois", "Mixed_mesoderm", "Cardiomyocytes", "Neural_crest"
                                   , "Endothelium", "NMP", "Erythroid", "Paraxial_mesoderm"
                                   , "ExE_endoderm", "Pharyngeal_mesoderm", "forebrain", "Somitic_mesoderm"
                                   , "Gut", "Spinal_cord", "Mesenchyme", "Surface_ectoderm", "mid_hindbrain")
                        , celltype_colour = c("#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442"
                                              , "#56B4E9", "#E69F00", "#009292", "#E41A1C", "#F781BF"
                                              , "#66C2A5", "#FC8D62", "#A6D854", "#5A4A7F", "#F4A460"
                                              , "#4682B4", "#E6B8B7"))

colour.df <- colour.df[!colour.df$celltype %in% c("forebrain", "mid_hindbrain"),]

rownames(colour.df) <- colour.df$celltype
colour.df <- colour.df[names(bom_rocs), ]


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


p <- ggroc(bom_rocs, size = 1) + my_theme +
  guides(linetype = guide_legend(override.aes = list(size = 3))) +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.8, color = "grey") +
  coord_equal() +
  # scale_colour_manual(values=roc_colours, aesthetics = c("colour", "fill")
  #                     , labels = c(names(bom_rocs))) +
  scale_color_manual(values=(colour.df$celltype_colour)) +
  labs(y= "Sensitivity", x = "Specificity")

pdf("mouseE8.5_nc_500bp_BALANCED_by_E8.25_500bp_BOM_binary_pred_ROCs.pdf")
p
dev.off()

pred_files <- list.files(path = ".", pattern = "(.*)_vs_other_vert5.0_qval0.5_(.*)_pred.txt$")

pred_li <- lapply(pred_files, read.table, header = T, stringsAsFactors = F)

# pred$actual <- ifelse(pred$celltype %in% celltype, 1, 0)
# add model
celltypes <- sub("_vs_.*", "", pred_files)
names(pred_li) <- celltypes

# actual class as factor
for(i in 1:length(pred_li)){
  pred_li[[i]]$actual <- factor(pred_li[[i]]$actual, levels=c(0,1))
}

pr_curves <- lapply(pred_li, pr_curve, truth=actual, raw,event_level = "second")
pr_curves <- lapply(pr_curves, as.data.frame)

for(i in 1:length(pr_curves)){
  pr_curves[[i]]$xgb_model <- names(pr_curves)[i]
}

pr_curves <- do.call("rbind", pr_curves)#26573     4

# match previous color code
colour.df <- data.frame(celltype=c("Allantois", "Mixed_mesoderm", "Cardiomyocytes", "Neural_crest"
                                   , "Endothelium", "NMP", "Erythroid", "Paraxial_mesoderm"
                                   , "ExE_endoderm", "Pharyngeal_mesoderm", "forebrain", "Somitic_mesoderm"
                                   , "Gut", "Spinal_cord", "Mesenchyme", "Surface_ectoderm", "mid_hindbrain")
                        , celltype_colour = c("#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442"
                                              , "#56B4E9", "#E69F00", "#009292", "#E41A1C", "#F781BF"
                                              , "#66C2A5", "#FC8D62", "#A6D854", "#5A4A7F", "#F4A460"
                                              , "#4682B4", "#E6B8B7"))

colour.df <- colour.df[!colour.df$celltype %in% c("forebrain", "mid_hindbrain"),]

rownames(colour.df) <- colour.df$celltype
colour.df <- colour.df[unique(pr_curves$xgb_model), ]

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

p <- ggplot(pr_curves, aes(x = recall, y = precision, color = xgb_model)) +
  geom_path() + coord_equal() +  my_theme +
  scale_color_manual(values=unique(colour.df$celltype_colour))
# scale_colour_manual(values = my_colours, aesthetics = c("colour"))

pdf("mouseE8.5_nc_500bp_BALANCED_by_E8.25_500bp_BOM_binary_pred_PRcurvess.pdf")
p
dev.off()
