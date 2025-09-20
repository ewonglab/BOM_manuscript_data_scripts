library(ggplot2)
library(tidyr)
library(stringr)

setwd("bom/processing_time")

# loading times
mouse_annot <- read.csv(file = "mouse_motif_annot_times.csv"
                        , header = T, stringsAsFactors = F)
mouse_train <- read.csv(file = "mouse_train_time.csv", header = T
                        , stringsAsFactors = F)
mouse_shap <- read.csv(file = "mouse_shap_time.csv", header = T
                       , stringsAsFactors = F)

celllines_annot <- read.csv(file = "celllines_motif_annot.csv", header = T, stringsAsFactors = F)
celllines_train <- read.csv(file = "celllines_train_time.csv"
                            , header = T, stringsAsFactors = F)
celllines_shap <- read.csv(file = "celllines_shap_time.csv"
                           , header = T, stringsAsFactors = F)

# read number of training and validation CREs
n_train_CREs <- read.table(file = "NCREs_celllines_and_e8.25.txt"
                           , header = T, stringsAsFactors = F, sep = '\t')

mouse_annot <- unique(mouse_annot[,c("celltypes", "time", "number.of.CRES", "percentage")])# 85  4
mouse_train <- unique(mouse_train[,c("celltypes", "time", "percentage")])# 85  3
mouse_shap <- unique(mouse_shap[,c("celltypes", "cpu.time", "percentage")])# 85  3

# add CRE numbers
celllines_annot <- unique(celllines_annot[,c("celltypes", "time", "percentage")])# 30  3
celllines_train <- unique(celllines_train[, c("celltypes", "time", "number.of.CRES", "percentage")])# 30  4
celllines_shap <- unique(celllines_shap[,c("celltypes", "time", "percentage")])# 30  3

celllines_annot$celltypes <- sub("\"", "", celllines_annot$celltypes)
celllines_annot <- merge(celllines_annot
                         , celllines_train[,c("celltypes", "percentage", "number.of.CRES")]
                         , by = c("celltypes", "percentage"), all.x = T)
celllines_train$number.of.CRES <- NULL

n_train_CREs <- tidyr::spread(n_train_CREs, V2, V3)
colnames(n_train_CREs) <- sub("V1", "celltypes", colnames(n_train_CREs))
colnames(n_train_CREs) <- sub("data_percentage", "percentage", colnames(n_train_CREs))

mouse_train <- merge(mouse_train, n_train_CREs[,c("celltypes", "percentage", "train_neg"
                                                  , "train_pos", "val_neg", "val_pos")]# 85  7
                     , by = c("celltypes", "percentage"), all.x = T)
celllines_train <- merge(celllines_train, n_train_CREs[,c("celltypes", "percentage", "train_neg"
                                                          , "train_pos", "val_neg", "val_pos")]
                         , by = c("celltypes", "percentage"), all.x = T)# 30  7

mouse_shap <- merge(mouse_shap, n_train_CREs[,c("celltypes", "percentage"
                                                , "train_neg", "train_pos")]
                    , by = c("celltypes", "percentage"), all.x = T)# 85  5
celllines_shap <- merge(celllines_shap, n_train_CREs[,c("celltypes", "percentage"
                                                        , "train_neg", "train_pos")]
                        , by = c("celltypes", "percentage"), all.x = T) # 30  5

split_time <- function(x){
  colnames(x) <- sub("cpu.time", "time", colnames(x))
  x$time <- sub("resources_used.cput=", "", x$time)
  time_df <- str_split_fixed(x$time, ":", 3)
  time_df <- as.data.frame(time_df)
  colnames(time_df) <- c("hours", "minutes", "seconds")
  time_df$hours <- as.numeric(time_df$hours)
  time_df$minutes <- as.numeric(time_df$minutes)
  time_df$seconds <- as.numeric(time_df$seconds)
  x <- cbind(x, time_df)
  return(x)
}

mouse_annot <- split_time(mouse_annot)# 85  7
mouse_train <- split_time(mouse_train)# 85 10
mouse_shap <- split_time(mouse_shap)# 85  8

celllines_annot <- split_time(celllines_annot)# 30  7
celllines_train <- split_time(celllines_train)# 30 10
celllines_shap <- split_time(celllines_shap)# 30  8

mouse_annot$all_seconds <- with(mouse_annot, (hours * 60 * 60) + (minutes * 60) + seconds)
mouse_train$all_seconds <- with(mouse_train, (hours * 60 * 60) + (minutes * 60) + seconds)
mouse_shap$all_seconds <- with(mouse_shap, (hours * 60 * 60) + (minutes * 60) + seconds)

celllines_annot$all_seconds <- with(celllines_annot, (hours * 60 * 60) + (minutes * 60) + seconds)
celllines_train$all_seconds <- with(celllines_train, (hours * 60 * 60) + (minutes * 60) + seconds)
celllines_shap$all_seconds <- with(celllines_shap, (hours * 60 * 60) + (minutes * 60) + seconds)

# maybe I can plot mouse time in seconds and human in minutes

# combine time in minutes
mouse_annot$all_minutes <- with(mouse_annot, all_seconds/60)
mouse_train$all_minutes <- with(mouse_train, all_seconds/60)
mouse_shap$all_minutes <- with(mouse_shap, all_seconds/60)

celllines_annot$all_minutes <- with(celllines_annot, all_seconds/60)
celllines_train$all_minutes <- with(celllines_train, all_seconds/60)
celllines_shap$all_minutes <- with(celllines_shap, all_seconds/60)

dot_colors <- data.frame(celltype = c("allantois", "mixed_meso", "cardiom"
                                      , "neuralcrest", "endothelium", "NMP"
                                      , "erythroid", "paraxial_meso"
                                      , "exe_endo", "pharyngeal_meso"
                                      , "forebrain", "somitic_meso"
                                      , "gut", "spinalcord", "mesenchyme"
                                      , "surface_ecto", "mid_hindbrain")
                         , colors = c("#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442"
                                      , "#56B4E9", "#E69F00", "#009292", "#E41A1C", "#F781BF", "#66C2A5"
                                      , "#FC8D62", "#A6D854", "#5A4A7F", "#F4A460", "#4682B4", "#E6B8B7"))

dot_colors <- dot_colors[with(dot_colors, order(celltype)),]

### VERSION 2: PLOT ONLY 100% DATA
mouse_annot <- mouse_annot[mouse_annot$percentage == 100, ]
mouse_train <- mouse_train[mouse_train$percentage == 100, ]
mouse_shap <- mouse_shap[mouse_shap$percentage == 100, ]

celllines_annot <- celllines_annot[celllines_annot$percentage == 100, ]
celllines_train <- celllines_train[celllines_train$percentage == 100, ]
celllines_shap <- celllines_shap[celllines_shap$percentage == 100, ]

# MOUSE MOTIF ANNOTATION
pdf("mouse_annot_seconds_100perc.1.pdf")
ggplot(mouse_annot, aes(x = number.of.CRES, y = all_seconds, color = celltypes)) + #
  geom_point() + geom_smooth(method=lm, aes(group=1), color="darkred") +  theme_classic() +
  theme(text = element_text(size=20)) + scale_colour_manual(values=dot_colors$colors
                                                            , aesthetics = c("colour"))
dev.off()

mouse_train$n_train <- with(mouse_train, train_neg + train_pos)
mouse_train$n_train_and_val <- with(mouse_train, train_neg + train_pos + val_neg + val_pos)

pdf("mouse_train_seconds_trainCREs_100perc.1.pdf")
ggplot(mouse_train, aes(x = n_train, y = all_seconds, color = celltypes)) +
  geom_point() + theme_classic() + geom_smooth(method=lm, aes(group=1), color="darkred") +
  theme(text = element_text(size=20)) + scale_colour_manual(values=dot_colors$colors
                                                            , aesthetics = c("colour"))
dev.off()

mouse_shap$n_train <- with(mouse_shap, train_neg + train_pos)
pdf("mouse_shap_seconds_trainCREs_100perc.1.pdf")
ggplot(mouse_shap, aes(x = n_train, y = all_seconds, color = celltypes)) +
  geom_point() + theme_classic() + geom_smooth(method=lm, aes(group=1), color="darkred") +
  theme(text = element_text(size=20)) + scale_colour_manual(values=dot_colors$colors
                                                            , aesthetics = c("colour"))
dev.off()

pdf("celllines_annot_minutes_100perc.1.pdf")
ggplot(celllines_annot, aes(x = number.of.CRES, y = all_minutes, color = celltypes)) +
  geom_point() + theme_classic() + geom_smooth(method=lm, aes(group=1), color="darkred") +
  theme(text = element_text(size=20))
dev.off()

celllines_train$n_train <- with(celllines_train, train_neg + train_pos)
celllines_train$n_train_and_val <- with(celllines_train, train_neg + train_pos + val_neg + val_pos)

pdf("celllines_train_minutes_trainCREs_100perc.1.pdf")
ggplot(celllines_train, aes(x = n_train, y = all_minutes, color = celltypes)) +
  geom_point() + theme_classic() + geom_smooth(method=lm, aes(group=1), color="darkred") +
  theme(text = element_text(size=20))
dev.off()

celllines_shap$n_train <- with(celllines_shap, train_neg + train_pos)
pdf("celllines_shap_minutes_trainCREs_100perc.1.pdf")
ggplot(celllines_shap, aes(x = n_train, y = all_minutes, color = celltypes)) +
  geom_point() + theme_classic() + geom_smooth(method=lm, aes(group=1), color="darkred") +
  theme(text = element_text(size=20))
dev.off()

