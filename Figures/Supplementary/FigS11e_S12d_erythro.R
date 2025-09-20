## ERITHROID/ERYTHROBLAST - ENHANCERS - liftOver with 500bp
library(ggplot2)
library(viridis)
library(ggplot2)
library(reshape2)
library(tidyr)
library("ggpubr")

setwd("/xgb/cross_sp")
minMatch <- c(0.95, 0.6)

human_ery <- read.table(file = "/data/human/fetal/spec_distal_nonexonic/Erythroblasts_MULTI_distal_nc_500bp.bed"
                        , header =F, stringsAsFactors = F, sep ='\t')
mouse_ery <- read.table(file = "/g/data/zk16/cc3704/mouse_data/gottgens_scATAC/erythroid_distal_nc_500bp.bed"
                        , header =F, stringsAsFactors = F, sep ='\t')

# add full coordinates
# id
human_ery$id <- paste0("id", 1:nrow(human_ery))
mouse_ery$id <- paste0("id", 1:nrow(mouse_ery))

# cross species mapping of enhancers
human_aligned_bed <- expand.grid("/data/human/fetal/spec_distal_nonexonic/Erythroblasts_MULTI_500bp_mm10_", minMatch, ".bed")
mouse_aligned_bed <- expand.grid("/data/mouse_erythroid_500bp_hg19_", minMatch, ".bed")

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

human_aligned <- do.call("rbind", human_aligned)#
mouse_aligned <- do.call("rbind", mouse_aligned)#

aligned <- rbind(human_aligned, mouse_aligned)
# all unique sequences in each set
align0.95 <- aligned[grepl("0.95", aligned$dataset),]
align0.6 <- aligned[grepl("0.6", aligned$dataset),]

unique(align0.95$dataset)
# [1] "Erythroblasts_MULTI_500bp_mm10_0.95" "mouse_erythroid_500bp_hg19_0.95"
unique(align0.6$dataset)
# [1] "Erythroblasts_MULTI_500bp_mm10_0.6" "mouse_erythroid_500bp_hg19_0.6"    

mouse_ery$align_0.95 <- ifelse(mouse_ery$id %in%
                                 align0.95[align0.95$dataset=="mouse_erythroid_500bp_hg19_0.95", "V4"]
                               ,1, 0)
mouse_ery$align_0.6 <- ifelse(mouse_ery$id %in%
                                align0.6[align0.6$dataset=="mouse_erythroid_500bp_hg19_0.6", "V4"]
                              ,1, 0)
human_ery$align_0.95 <- ifelse(human_ery$id %in%
                                 align0.95[align0.95$dataset=="Erythroblasts_MULTI_500bp_mm10_0.95", "V4"]
                               ,1, 0)
human_ery$align_0.6 <- ifelse(human_ery$id %in%
                                align0.6[align0.6$dataset=="Erythroblasts_MULTI_500bp_mm10_0.6", "V4"]
                              ,1, 0)


mouse_by_humanBOM <- read.table(file = "mouse_eryth_vs_other_byFetalErythroblasts_vert5.0_q0.5_500bp.txt"
                                , header = T, stringsAsFactors = F)
human_by_mouseBOM <- read.table(file = "human_Erythroblasts_vs_other_byMouseE8.25Eryth_vert5.0_q0.5_500bp.txt"
                                , header = T, stringsAsFactors = F)
# ONLY ACTUAL CARDIOMYOCYTE ENHANCERS
mouse_by_humanBOM <- mouse_by_humanBOM[mouse_by_humanBOM$actual ==1,]#1002    3
human_by_mouseBOM <- human_by_mouseBOM[human_by_mouseBOM$actual ==1,]#140   3

# add peak id
mouse_ery$peakid <- with(mouse_ery, paste(V1, paste(V2, V3, sep = "-"), sep = ":"))
human_ery$peakid <- with(human_ery, paste(V1, paste(V2, V3, sep = "-"), sep = ":"))

# crop_enh <- function(x, n){
#   x$width <- with(x, V3-V2)
#   x$centre <- with(x, V2 + round(width/2))
#   x$start <- round(x$centre - n/2)
#   x$end <- round(x$centre + n/2)
#   print(unique(with(x, end-start)))
#   x$peakid <- with(x, paste(V1, paste(start, end, sep = "-"), sep = ":"))
#   x$start <- NULL
#   x$end <- NULL
#   x$centre <- NULL
#   return(x)
# }
#
# # adding peak ids that correspond to 500bp peaks
# mouse_ery <- crop_enh(mouse_ery, 500)
# # [1] 500
# human_ery <- crop_enh(human_ery, 500)
# # [1] 500

mouse_ery <- merge(mouse_ery, mouse_by_humanBOM[,c("predicted", "raw")]
                   , by.x="peakid", by.y = 0)#[1] 1002   10
human_ery <- merge(human_ery, human_by_mouseBOM[,c("predicted", "raw")]
                   , by.x="peakid", by.y = 0) #[1] 140   9

# now I need to make tables with the percentages of enhancers alighning
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


mouse_freq <- calc_perc(mouse_ery)
human_freq <- calc_perc(human_ery)

mouse_freq$species <- "mouse"
human_freq$species <- "human"

mouse_human <- rbind(mouse_freq, human_freq)

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
  x <- tidyr::spread(x, align, Freq, fill = 0)
  # print(x)
  x$align1_perc <- with(x, align1/(align1 + align0) * 100)
  x$align0_perc <- with(x, align0/(align1 + align0) * 100)
  x <- x[,c("predicted", "align1_perc", "align0_perc")]
 
  x$predicted <- paste0("class", x$predicted)
  x <- melt(x, "predicted")
  # print(x)
  x$value <- ifelse(is.na(x$value), 0, x$value)
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
# 1    class0 align1_perc  5.769231         5.8
# 2    class1 align1_perc  3.409091         3.4
# 3    class0 align0_perc 94.230769        94.2
# 4    class1 align0_perc 96.590909        96.6

human_0.6_pie <- make_piecharts(human_0.6)
# predicted    variable value value_round
# 1    class0 align1_perc  50.0        50.0
# 2    class1 align1_perc  62.5        62.5
# 3    class0 align0_perc  50.0        50.0
# 4    class1 align0_perc  37.5        37.5

pdf("Erythroblasts_alignment_to_mm10_500bpModels.1.pdf", onefile = F)
ggarrange(plotlist = list(human_0.95_pie, human_0.6_pie)
          , common.legend = T, nrow=2, ncol=1)
dev.off()

mouse_0.95_pie <- make_piecharts(mouse_0.95)
# predicted    variable    value value_round
# 1    class0 align1_perc 16.61342        16.6
# 2    class1 align1_perc 22.64151        22.6
# 3    class0 align0_perc 83.38658        83.4
# 4    class1 align0_perc 77.35849        77.4

mouse_0.6_pie <- make_piecharts(mouse_0.6)
# predicted    variable    value value_round
# 1    class0 align1_perc 55.27157        55.3
# 2    class1 align1_perc 58.78084        58.8
# 3    class0 align0_perc 44.72843        44.7
# 4    class1 align0_perc 41.21916        41.2

pdf("mouse_erythroid_alignment_to_hg19_500bpModels.1.pdf", onefile = F)
ggarrange(plotlist = list(mouse_0.95_pie
                          , mouse_0.6_pie)
          , common.legend = T, ncol=1, nrow=2)
dev.off()

print_perc <- function(x,y){
  print(round((x/(x+y))*100, 1))
  print(round((y/(x+y))*100, 1))
}
