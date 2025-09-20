# Summary of statistics for binary predictors - BOM, kgm-svm, Enformer and DNABERT including MCC
library(ggplot2)

binary_class_stats <- read.csv(file = "./data/mouse/binary_classifiers_stats.csv"
                               , header = T, stringsAsFactors = F)

dnabert <- binary_class_stats[binary_class_stats$Model == "DNABERT",]
enformer <- binary_class_stats[binary_class_stats$Model == "Enformer",]

library(reshape2)
dnabert <- melt(dnabert, c("Cell.type", "Model"))
enformer <- melt(enformer, c("Cell.type", "Model"))

gkm <- read.table(file = "./gkmsvm/e8.25_june_2023/t4_l11_k7_d3_model_stats_500bp_mcc"
                  , header = F, stringsAsFactors = F, sep = '\t')
bom <- read.table(file = "xgb/e8.25/bin_500bp/BOM_binary_500bp_stats_full"
                  , header = F, stringsAsFactors = F, sep = '\t')

gkm$V3 <- round(gkm$V3, 3)
bom$V3 <- round(bom$V3, 3)

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

dnabert_summary <- data_summary(dnabert, "value", "variable")

enformer_summary <- data_summary(enformer, "value", "variable")

colnames(gkm) <- c("celltype", "variable", "value")
colnames(bom) <- c("celltype", "variable", "value")

gkm_summary <- data_summary(gkm, "value", "variable")

bom_summary <- data_summary(bom, "value", "variable")

dnabert_summary$Model <- "DNABERT"
enformer_summary$Model <- "Enformer"
gkm_summary$Model <- "gkm"
bom_summary$Model <- "BOM"

all_stats <- rbind(dnabert_summary
                   , enformer_summary
                   , gkm_summary
                   , bom_summary)

all_stats$Model <- factor(all_stats$Model
                             , levels = c("gkm", "DNABERT"
                                          , "Enformer", "BOM"))
all_stats$variable <- as.character(all_stats$variable)
all_stats$variable <- gsub("PRAUC", "auPR", all_stats$variable)
all_stats$variable <- gsub("ROC_AUC", "auROC", all_stats$variable)

all_stats <- all_stats[all_stats$variable %in% c("Accuracy", "F1", "auPR", "auROC", "MCC"),]

method_colors <- c("#AAD9BC", "#8B0000", "#FF4500", "#B89BC9")

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

all_stats$Model <- factor(all_stats$Model
                          , levels = c("gkm", "DNABERT"
                                       , "Enformer", "BOM"))


p <- ggplot(all_stats, aes(x = Model, y = value, fill = Model))+
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = method_colors) +
  my_theme + facet_wrap(~variable)

pdf("binary_stats_mean_with_errorbars_revised.pdf")
p
dev.off()
