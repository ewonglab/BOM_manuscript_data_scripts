# Beeswarm plot and heatmap for the prediction statistics of binary models (DNABERT, gkm-SVM and BOM) - 500bp enhancers
library(tidyr)
all_stats <- read.table(file = "stats_500bp_enh", header = F
                        , stringsAsFactors = F, sep ='\t')

all_stats.spr <- tidyr::spread(all_stats, V2, V3)
all_stats.spr <- all_stats.spr[with(all_stats.spr, order(V4)),]
write.csv(x = all_stats.spr, file = "stats_500bp_enh_byCelltype.csv", row.names = F, quote = F)

all_stats_summ <- aggregate(V3~., all_stats[,1:3], mean)
all_stats_summ.spr <- tidyr::spread(all_stats_summ, V2, V3)

write.csv(x = all_stats_summ.spr, file = "stats_500bp_enh_summary.csv", row.names = F, quote = F)


library(pheatmap)
library(RColorBrewer)
# heatmap_trial_2 <- read.csv("Final genes_log2CPM.csv")
# heatmap_trial_2 <- data.frame(heatmap_trial_2[,-1], row.names=heatmap_trial_2[,1])
# sc_1 <-t(scale(t(heatmap_trial_2), center = TRUE, scale = TRUE))

rownames(all_stats_summ.spr) <- all_stats_summ.spr$V1
all_stats_summ.spr$V1 <- NULL

pdf("binary_stats_summary.pdf")
pheatmap(all_stats_summ.spr, scale = "none"
         , cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = TRUE, show_colnames = TRUE,
         colorRampPalette(brewer.pal(9,"BuPu"))(100))
dev.off()


library(data.table)
library(ggplot2)
library(reshape2)
library(ggbeeswarm)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func,
                    varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

my_theme <-  theme(panel.border = element_blank(), panel.grid.major = element_blank()
                   , panel.grid.minor = element_blank()
                   , axis.line = element_line(colour = "black", size = 1)
                   , legend.title = element_blank()
                   , legend.key=element_blank()
                   , legend.position = "none"
                   , legend.key.width = unit(1.5,"cm")
                   , panel.background = element_blank()
                   , text = element_text(size=10)
                   , axis.text.x=element_text(colour="black"
                                              , angle = 90, vjust = 0.5, hjust=1)
                   , axis.text.y=element_text(colour="black")
                   # , legend.text=element_text(size=10)
                   , axis.ticks = element_line(colour = "black"))


my_colours <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442"
                , "#56B4E9", "#E69F00", "#009292", "#E41A1C", "#F781BF", "#66C2A5"
                , "#FC8D62", "#A6D854", "#5A4A7F", "#F4A460", "#4682B4", "#E6B8B7")

celltype_order <- c("allantois", "mixed_meso", "cardiom", "neuralcrest", "endothelium", "NMP"
                    , "erythroid", "paraxial_meso", "exe_endo", "pharyngeal_meso"
                    , "forebrain", "somitic_meso", "gut", "spinalcord", "mesenchyme"
                    , "surface_ecto", "mid_hindbrain")

colour.df <- data.frame(celltype=celltype_order, celltype_colour=my_colours)
rownames(colour.df) <- colour.df$celltype
colour.df <- colour.df[with(colour.df, order(celltype)),]

p <- ggplot(all_stats, mapping=aes(V1,V3,color=V4)) +
  geom_beeswarm(dodge.width=.8,cex=2) + facet_wrap(~V2) + my_theme +
  stat_summary(
    fun = "mean", geom = "point", aes(group = interaction(V1, V2)),
    color = "black", size = 5, shape=95) +  scale_color_manual(values=colour.df$celltype_colour) +
  scale_x_discrete(expand = c(0.1, 0.1))

pdf("E8.25_binary_methods_ggbeeswarm.pdf")
p
dev.off()

