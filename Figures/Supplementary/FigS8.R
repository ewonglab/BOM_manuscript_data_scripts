library(tidyr)
library(ggplot2)

setwd("xgb/e8.25")
all_motifs_stats <- read.table(file = "E8.25_binary_stats_full", header = F
                               , stringsAsFactors = F, sep = '\t')#119   3

perc_stats <- read.table(file = "./sample_Ovlp/binary_stats_motif_perc"
                         , header = F, sep = '\t', stringsAsFactors = F)

perc_stats$celltype <- sub("_.*", "", perc_stats$V1)
perc_stats$percentage <- sub(".*_sample_", "", perc_stats$V1)
perc_stats$percentage <- sub("_pred.*", "", perc_stats$percentage)
perc_stats$percentage <- sub("_.*", "", perc_stats$percentage)

all_motifs_stats <- all_motifs_stats[all_motifs_stats$V1 %in% c("endothelium", "neuralcrest"), ]

colnames(all_motifs_stats)[1] <- c("celltype")
all_motifs_stats$percentage <- "1"

combined <- rbind(perc_stats[,c("celltype", "V2", "V3", "percentage")]
                  , all_motifs_stats[,c("celltype", "V2", "V3", "percentage")])
combined$percentage <- as.numeric(combined$percentage)

pdf("stats_by_motif_perc.pdf")
ggplot(data = combined, aes(x = percentage, y = V3, color=celltype))+
  geom_line() + theme_classic() + facet_wrap(.~V2, nrow = 2) +
  scale_x_continuous(breaks=seq(0.1, 1, 0.1))
dev.off()
