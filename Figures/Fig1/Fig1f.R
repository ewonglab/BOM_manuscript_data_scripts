# Multilabel prediction stats barplot with error bars

library(ggplot2)
library(reshape2)

stats <- read.table(file = "E8.25_multilab_vert5.0_qval0.5_merror_pred_stats.txt"
                    , header = T, stringsAsFactors = F, sep = '\t')

round(apply(stats[, 2:ncol(stats)], 2, mean, na.rm = T), 2)
#   ROCAUC  accuracy     PRAUC        F1 Precision    Recall       MCC
 #    0.83      0.94      0.44      0.30      0.66      0.20      0.32

stats <- melt(stats, "celltype")

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm = TRUE),
      sd = sd(x[[col]], na.rm = TRUE))
  }
  data_sum <- ddply(data, groupnames, .fun = summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

stats_summ <- data_summary(stats, "value", "variable")

pdf("E8.25_multilabel_stats_barplot.pdf")
ggplot(stats_summ, aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat = "identity", color = "black",
           position = position_dodge()) +
  geom_errorbar(aes(ymin = value-sd, ymax = value+sd), width = .2,
                position = position_dodge(.9))
dev.off()
