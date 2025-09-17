# zebrafish ROC and PR curves

suppressMessages({
  library("yardstick")
  library("cvAUC")
  library(pROC)
  library(ggplot2)
})

pred_files <- list.files(path = ".", pattern = "(.*)_v_other_nc_(.*)_pred.txt")
pred_li <- lapply(pred_files, read.table, header = T, stringsAsFactors = F)

tissues <- sub("_v_other_nc_.*", "", pred_files)
names(pred_li) <- tissues

bom_rocs <- lapply(pred_li, function(x) roc(x$actual, x$raw, direction="<"))
names(bom_rocs) <- tissues

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

my_colours <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"
                , "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#393b79")

p <- ggroc(bom_rocs, size = 1) + my_theme +
  guides(linetype = guide_legend(override.aes = list(size = 3))) +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.8, color = "grey") + coord_equal() +
  scale_colour_manual(values=my_colours, aesthetics = c("colour", "fill")
                      , labels = c(names(bom_rocs))) +
  labs(y= "Sensitivity", x = "Specificity")

pdf("zebrafish_nc_BOM_binary_pred_ROCs.pdf")
p
dev.off()

pred_files <- list.files(path = ".", pattern = "(.*)_v_other_nc_(.*)_pred.txt")
pred_li <- lapply(pred_files, read.table, header = T, stringsAsFactors = F)

tissues <- sub("_v_other_nc_.*", "", pred_files)
names(pred_li) <- tissues

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

my_colours <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"
                , "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#393b79")

p <- ggplot(pr_curves, aes(x = recall, y = precision, color = xgb_model)) +
  geom_path() + coord_equal() +  my_theme +
  scale_colour_manual(values = my_colours, aesthetics = c("colour"))

pdf("zebrafish_binary_PRcurves.pdf")
p
dev.off()
