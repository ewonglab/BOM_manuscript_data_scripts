library("yardstick")
library("cvAUC")
library(pROC)
library(ggplot2))
library("ggpubr")
library(reshape2)
library(ggbeeswarm)

my_theme <-  theme(panel.border = element_blank(), panel.grid.major = element_blank()
                   , panel.grid.minor = element_blank()#legend.position="none"
                   , axis.line = element_line(colour = "black", size = 1)
                   , legend.title = element_blank()
                   , legend.key=element_blank()
                   , legend.position = "bottom"
                   , legend.key.width = unit(1.5,"cm")
                   , panel.background = element_blank()
                   , text = element_text(size=12)
                   , legend.text=element_text(size=8)
                   , axis.text.x=element_text(colour="black")
                   , axis.text.y=element_text(colour="black")
                   , axis.ticks = element_line(colour = "black"))

BOM_stats_df <- read.csv(file = "BOM_multiclass_Topics_stats_mlogloss"
                         , header = T, stringsAsFactors = F, row.names = 1)

n_topics <- length(unique(BOM_stats_df$Topic))
my_colors <- rainbow(n_topics)

BOM_stats_df <- melt(BOM_stats_df, "Topic")


p <- ggplot(BOM_stats_df, mapping=aes(variable,value,color=Topic)) +
  geom_beeswarm(dodge.width=.8,cex=2) + my_theme +
  stat_summary(
    fun = "mean", geom = "point", aes(group = variable),
    color = "black", size = 10, shape=95) +  scale_color_manual(values=my_colors) +
  scale_x_discrete(expand = c(0.1, 0.1))

pdf("topics_multiclass_stats_ggbeeswarm_mlogloss.pdf")
p
dev.off()

BOM_stats_summ <- aggregate(value~variable, BOM_stats_df, mean)
BOM_stats_summ$value <- round(BOM_stats_summ$value, 3)
