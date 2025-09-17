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

stats_files <- list.files(path = ".", pattern = "BOM_binary_seed(.*)_stats")

seed_stats <- lapply(stats_files, fread)
names(seed_stats) <- sub(".*seed", "", stats_files)
names(seed_stats) <- sub("_stats", "", names(seed_stats))

# as data frame
seed_stats <- lapply(seed_stats, as.data.frame)

## add seed number as column
for(i in 1:length(seed_stats)){
  seed_stats[[i]]$random_seed <- names(seed_stats)[[i]]
}

# list to data frame
seed_stats <- do.call("rbind", seed_stats)#510   4
colnames(seed_stats) <- c("celltype", "stat", "value", "random_seed")

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

p <- ggplot(seed_stats, mapping=aes(celltype,value,color=celltype)) +
  geom_beeswarm(dodge.width=.8,cex=2) + facet_wrap(~stat) + my_theme +
  stat_summary(
    fun = "mean", geom = "point", aes(group = interaction(celltype, stat)),
    color = "black", size = 5, shape=95) +  scale_color_manual(values=colour.df$celltype_colour) +
  scale_x_discrete(expand = c(0.1, 0.1))

pdf("E8.25_binary_stats_randomSeeds_ggbeeswarm.pdf")
p
dev.off()
