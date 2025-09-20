suppressMessages({
  library(tidyr)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

setwd("/xgb/e8.25")
bom <- read.table(file = "test_multiclass_500bp_4597T_pred.txt", header = T
                  , stringsAsFactors = F)
conf.df <- as.data.frame(table(bom[, c("predicted", "label")]))
conf.mat <- tidyr::spread(data = conf.df, key = predicted, value = Freq)

conf.df <- conf.df %>%
  group_by(label) %>%
  mutate(Proportion = Freq / sum(Freq))

celltypes <- unique(c(conf.df$label, conf.df$predicted))
celltypes[order(tolower(celltypes))]
conf.df$label <- factor(conf.df$label, levels = rev(celltypes[order(tolower(celltypes))]))
conf.df$predicted <- factor(conf.df$predicted, levels = (celltypes[order(tolower(celltypes))]))

ggplot(conf.df, aes(x = predicted, y = label, fill = Proportion)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Proportion") +
  theme_minimal() +
  theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12), 
      panel.grid = element_blank(),
      axis.title = element_blank()
  )
