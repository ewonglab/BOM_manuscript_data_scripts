library(ggplot2)
library(viridis)
library(ggplot2)
library(reshape2)
library(tidyr)
library("ggpubr")

minMatch <- c(0.95, 0.6)

human_cardi <- read.table(file = "Cardiomyocytes_HEART_distal_nc_500bp.bed"
                          , header =F, stringsAsFactors = F, sep ='\t')#581   3
mouse_cardi <- read.table(file = "cardiom_distal_nc_500bp.bed"
                          , header =F, stringsAsFactors = F, sep ='\t')#869   5

# add full coordinates
# id
human_cardi$id <- paste0("id", 1:nrow(human_cardi))
mouse_cardi$id <- paste0("id", 1:nrow(mouse_cardi))

# cross species mapping of enhancers
human_aligned_bed <- expand.grid("Cardiomyocytes_HEART_500bp_mm10_", minMatch, ".bed")
mouse_aligned_bed <- expand.grid("mouse_cardiom_500bp_hg19_", minMatch, ".bed")

human_aligned <- lapply(with(human_aligned_bed, paste0(Var1, Var2, Var3))
                        , read.table, header = F, stringsAsFactors = F)
names(human_aligned) <- basename(with(human_aligned_bed, paste0(Var1, Var2)))

mouse_aligned <- lapply(with(mouse_aligned_bed, paste0(Var1, Var2, Var3))
                        , read.table, header = F, stringsAsFactors = F)
names(mouse_aligned) <- basename(with(mouse_aligned_bed, paste0(Var1, Var2)))

for(i in 1:length(human_aligned)){
  human_aligned[[i]]$dataset <- names(human_aligned)[i]
}

for(i in 1:length(mouse_aligned)){
  mouse_aligned[[i]]$dataset <- names(mouse_aligned)[i]
}

human_aligned <- do.call("rbind", human_aligned)#454   5
mouse_aligned <- do.call("rbind", mouse_aligned)#932   5

# I am not overestimating the number of CREs that map to the other genome becuase below I am setting the
# variables align_0.95 and align_0.6 which depend on the presence of the original ids in the mapped regions

aligned <- rbind(human_aligned, mouse_aligned)
# all unique sequences in each set
align0.95 <- aligned[grepl("0.95", aligned$dataset),]
align0.6 <- aligned[grepl("0.6", aligned$dataset),]

unique(align0.95$dataset)
# [1] "Cardiomyocytes_HEART_500bp_mm10_0.95"
# [2] "mouse_cardiom_500bp_hg19_0.95"
unique(align0.6$dataset)
# [1] "Cardiomyocytes_HEART_500bp_mm10_0.6" "mouse_cardiom_500bp_hg19_0.6"

mouse_cardi$align_0.95 <- ifelse(mouse_cardi$id %in%
                                   align0.95[align0.95$dataset=="mouse_cardiom_500bp_hg19_0.95", "V4"]
                                 ,1, 0)
mouse_cardi$align_0.6 <- ifelse(mouse_cardi$id %in%
                                  align0.6[align0.6$dataset=="mouse_cardiom_500bp_hg19_0.6", "V4"]
                                ,1, 0)
human_cardi$align_0.95 <- ifelse(human_cardi$id %in%
                                   align0.95[align0.95$dataset=="Cardiomyocytes_HEART_500bp_mm10_0.95", "V4"]
                                 ,1, 0)
human_cardi$align_0.6 <- ifelse(human_cardi$id %in%
                                  align0.6[align0.6$dataset=="Cardiomyocytes_HEART_500bp_mm10_0.6", "V4"]
                                ,1, 0)

mouse_by_humanBOM <- read.table(file = "mouse_cm_vs_other_byHumanFetal_vert5.0_q0.5_500bp.txt"
                                , header = T, stringsAsFactors = F)
human_by_mouseBOM <- read.table(file = "human_cm_vs_other_byMouseE8.25_vert5.0_q0.5_500bp.txt"
                                , header = T, stringsAsFactors = F)
# ONLY ACTUAL CARDIOMYOCYTE ENHANCERS
mouse_by_humanBOM <- mouse_by_humanBOM[mouse_by_humanBOM$actual ==1,]#869   3
human_by_mouseBOM <- human_by_mouseBOM[human_by_mouseBOM$actual ==1,]#581   3

### add binary predictions
### here I need to change the ids of mouse_cardi and human_cardi and crip them to 500bp

mouse_cardi$peakid <- with(mouse_cardi, paste(V1, paste(V2, V3, sep = "-"), sep = ":"))
human_cardi$peakid <- with(human_cardi, paste(V1, paste(V2, V3, sep = "-"), sep = ":"))

mouse_cardi <- merge(mouse_cardi, mouse_by_humanBOM[,c("predicted", "raw")]
                     , by.x="peakid", by.y = 0)#869  10
human_cardi <- merge(human_cardi, human_by_mouseBOM[,c("predicted", "raw")]
                     , by.x="peakid", by.y = 0) #581  10
library(reshape2)

calc_perc <- function(df){
  enh_freq <- as.data.frame(table(df[,c("predicted", "align_0.95", "align_0.6")]))
  enh_freq <- enh_freq[enh_freq$Freq != 0,]
 
  enh_freq <- melt(enh_freq, c("predicted", "Freq"))
  enh_freq <- aggregate(Freq ~., enh_freq, sum)
  tmp <- aggregate(Freq ~ ., enh_freq[,c("value", "variable", "Freq")], sum)
  # print(tmp)
  colnames(tmp)[3] <- "total"
  enh_freq <- merge(enh_freq, tmp, by =c("value", "variable"))
  enh_freq$perc <- with(enh_freq, Freq/total)
  return(enh_freq)
}


mouse_freq <- calc_perc(mouse_cardi)

human_freq <- calc_perc(human_cardi)

mouse_freq$species <- "mouse"
human_freq$species <- "human"

mouse_human <- rbind(mouse_freq, human_freq)
mouse_human
# value   variable predicted Freq total      perc species
# 1      0  align_0.6         0  118   245 0.4816327   mouse
# 2      0  align_0.6         1  127   245 0.5183673   mouse
# 3      0 align_0.95         0  259   561 0.4616756   mouse
# 4      0 align_0.95         1  302   561 0.5383244   mouse
# 5      1  align_0.6         0  282   624 0.4519231   mouse
# 6      1  align_0.6         1  342   624 0.5480769   mouse
# 7      1 align_0.95         0  141   308 0.4577922   mouse
# 8      1 align_0.95         1  167   308 0.5422078   mouse
# 9      0  align_0.6         0   64   176 0.3636364   human
# 10     0  align_0.6         1  112   176 0.6363636   human
# 11     0 align_0.95         0  194   532 0.3646617   human
# 12     0 align_0.95         1  338   532 0.6353383   human
# 13     1  align_0.6         0  142   405 0.3506173   human
# 14     1  align_0.6         1  263   405 0.6493827   human
# 15     1 align_0.95         0   12    49 0.2448980   human
# 16     1 align_0.95         1   37    49 0.7551020   human


### MAKE PIECHARTS

human_0.95 <- subset(mouse_human, species=="human" & variable=="align_0.95")
human_0.6 <- subset(mouse_human, species=="human" & variable=="align_0.6")

mouse_0.95 <- subset(mouse_human, species=="mouse" & variable=="align_0.95")
mouse_0.6 <- subset(mouse_human, species=="mouse" & variable=="align_0.6")

make_piecharts <- function(x){
  #remove barplot totals and percentages
  x$total <- NULL
  x$perc <- NULL
 
  x$align <- paste0("align", x$value)
  x$variable <- NULL
  x$species <- NULL
  x$value <- NULL
  x <- tidyr::spread(x, align, Freq)
  x$align1_perc <- with(x, align1/(align1 + align0) * 100)
  x$align0_perc <- with(x, align0/(align1 + align0) * 100)
  x <- x[,c("predicted", "align1_perc", "align0_perc")]
 
  x$predicted <- paste0("class", x$predicted)
  x <- melt(x, "predicted")
  x$predicted <- factor(x$predicted, levels = c("class1","class0"))
  x$variable <- factor(x$variable, levels = c("align1_perc","align0_perc"))
  x$value_round <- with(x, round(value, 1))
  print(x)
  p <- ggplot(x, aes(x="", y=value, fill=variable)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) + facet_wrap(.~predicted) +
    theme_void() + scale_fill_manual(values = c("#FFA500", "#1F77B4"))# +
  # theme(legend.position="bottom")
  return(p)
}


human_0.95_pie <- make_piecharts(human_0.95)
# predicted    variable     value value_round
# 1    class0 align1_perc  5.825243         5.8
# 2    class1 align1_perc  9.866667         9.9
# 3    class0 align0_perc 94.174757        94.2
# 4    class1 align0_perc 90.133333        90.1
human_0.6_pie <- make_piecharts(human_0.6)
# predicted    variable    value value_round
# 1    class0 align1_perc 68.93204        68.9
# 2    class1 align1_perc 70.13333        70.1
# 3    class0 align0_perc 31.06796        31.1
# 4    class1 align0_perc 29.86667        29.9

pdf("Cardiomyocytes_HEART_alignment_to_mm10_500bpModels.1.pdf", onefile = F)
ggarrange(plotlist = list(human_0.95_pie, human_0.6_pie)
          , common.legend = T, nrow=2, ncol=1)
dev.off()

mouse_0.95_pie <- make_piecharts(mouse_0.95)
# predicted    variable    value value_round
# 1    class0 align1_perc 35.25000        35.2
# 2    class1 align1_perc 35.60768        35.6
# 3    class0 align0_perc 64.75000        64.8
# 4    class1 align0_perc 64.39232        64.4
mouse_0.6_pie <- make_piecharts(mouse_0.6)
# predicted    variable    value value_round
# 1    class0 align1_perc 70.50000        70.5
# 2    class1 align1_perc 72.92111        72.9
# 3    class0 align0_perc 29.50000        29.5
# 4    class1 align0_perc 27.07889        27.1

pdf("mouse_cardiom_alignment_to_hg19_500bpModels.1.pdf", onefile = F)
ggarrange(plotlist = list(mouse_0.95_pie
                          , mouse_0.6_pie)
          , common.legend = T, ncol=1, nrow=2)
dev.off()
