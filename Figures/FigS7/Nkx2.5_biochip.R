
Mean biochip signal across q-value thresholds: Nkx2-5 motifs

# OVERLAP MM10 SUMMIT COORDINATES WITH MOTIFS
library("Hmisc", lib = "/g/data/zk16/software/Rpackages_paola/R_4.0.0")
library(GenomicRanges)
library(ggplot2)
library(viridis)

setwd("/g/data/zk16/cc3704/mouse_enh_grammar/data/mouse/e12.5")

## MOTIF: SRF
tf <- "Nkx2-5"
mm10_summits <- read.table(file = "chipseq_mouse_heart_fetal_mm10_0.bed"
                           , header = F, stringsAsFactors = F, sep = '\t')
mm10_summits$TF <- sub(".*_", "", mm10_summits$V4)
mm10_summits$id <- sub("_.*", "", mm10_summits$V4)
# read biochip signal for TF
mm10_summits <- mm10_summits[mm10_summits$TF == tf, ] # 25859     6

Nkx2.5_biochip <- read.table(file = "./biochip_mm9/Nkx2.5_combined_mean", header = T
                          , stringsAsFactors = F, sep ='\t') # 25864     4
# add mean signal to summits merging by mm19 id
mm10_summits <- merge(mm10_summits, Nkx2.5_biochip[,c("id", "mean")], by = "id")
table(is.na(mm10_summits$mean))
# FALSE
# 25859

# reading cardiomyocyte gimme
gimme <- read.table(file = "/g/data/zk16/cc3704/mouse_data/gottgens_scATAC/fimo/gimme/cardiom_Gimme_vertv5.0/fimo.tsv"
                    , header = T, stringsAsFactors = F, sep = '\t')
gimme <- gimme[gimme$q.value <= 0.5, ]

# read annotation and keep only Srf motifs
gimme_annot <- read.table(file = "/g/data/zk16/useful/gimmemotifs/gimme.vertebrate.v5.0.motif2factors.txt"
                          , header = T, stringsAsFactors = F, sep ='\t')

Nkx2.5_motifs <- unique(gimme_annot[gimme_annot$Factor %in% c("Nkx2-5", "NKX2-5", "Nkx2.5")
                                 , c("Factor", "Motif")])
(Nkx2.5_motifs)
# Factor                   Motif
# 1189 NKX2-5 GM.5.0.Homeodomain.0012
# 2947 NKX2-5 GM.5.0.Homeodomain.0033
# 4264 Nkx2.5 GM.5.0.Homeodomain.0054
# 4277 NKX2-5 GM.5.0.Homeodomain.0054
# 4278 Nkx2-5 GM.5.0.Homeodomain.0054
# 4804 Nkx2-5 GM.5.0.Homeodomain.0066
# 5352 Nkx2-5 GM.5.0.Homeodomain.0070
# 5841 Nkx2-5 GM.5.0.Homeodomain.0077
# 9147 NKX2-5 GM.5.0.Homeodomain.0168
# 9528 NKX2-5 GM.5.0.Homeodomain.0178
# 9674 NKX2-5 GM.5.0.Homeodomain.0183


gimme_annot.agr <- aggregate(Factor ~ Motif, gimme_annot, paste, collapse = ";")

gimme_annot.agr[gimme_annot.agr$Motif %in% Nkx2.5_motifs$Motif, ]
#                        Motif
# 840  GM.5.0.Homeodomain.0012
# 861  GM.5.0.Homeodomain.0033
# 882  GM.5.0.Homeodomain.0054
# 894  GM.5.0.Homeodomain.0066
# 898  GM.5.0.Homeodomain.0070
# 905  GM.5.0.Homeodomain.0077
# 996  GM.5.0.Homeodomain.0168
# 1006 GM.5.0.Homeodomain.0178
# 1011 GM.5.0.Homeodomain.0183
# Factor
# 840  NKX2-3;Nkx3-1;NKX2-4;NKX2-6;Nkx2-6;NKX2-3;NKX2-4;NKX2-6;NKX2-5;NKX3-2;Nkx3-1;NKX3-1;Nkx3-1;NKX3-1;Nkx3-1;Nkx3-2;NKX3-2;NKX3-1;Nkx3-1;NKX3-2;Nkx3-2;Bapx1
# 861                                                                              NKX2-3;NKX2-4;NKX2-6;NKX2-3;NKX2-4;NKX2-6;NKX2-5;NKX2-3;NKX2-6;NKX2-4;NKX2-5
# 882                      Nkx2.5;Nkx2.2;NKX2-3;NKX2-4;NKX2-3;NKX2-4;NKX2-3;NKX2-6;NKX2-4;NKX2-8;Nkx2-1;Nkx2-2;Nkx2-4;NKX2-5;Nkx2-5;NKX2-6;Nkx2-6;NKX2-8;Nkx2.1
# 894                                                                                                                 NKX2-3;NKX2-6;NKX2-3;NKX2-3;Nkx2-2;Nkx2-5
# 898                                                                               NKX2-3;NKX2-4;NKX2-6;Nkx2-5;HOXA6;HOXD3;AC012531.1;HOXB4;HOXC5;Hoxd4;Nkx2-5
# 905                                                                                                          NKX2-3;NKX2-4;NKX2-6;Nkx2-5;Nkx2-5;NKX2-1;Nkx2-1
# 996                                                                                                                                                    NKX2-5
# 1006                                                                                                                                                   NKX2-5
# 1011                                                                                                                       NKX2-3;NKX2-4;NKX2-6;NKX2-5;NKX2-5

gimme <- gimme[gimme$motif_id %in% Nkx2.5_motifs$Motif,] #
gimme$q.value_bin <- ""
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.01, 1, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.05 & gimme$q.value > 0.01, 2, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.1 & gimme$q.value > 0.05, 3, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.3 & gimme$q.value > 0.1, 4, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.5 & gimme$q.value > 0.3, 5, gimme$q.value_bin)

table(gimme$q.value_bin)
#   4   5
# 204 825

gimme$q.value_bin <- factor(x = gimme$q.value_bin, levels = c(1, 2, 3, 4, 5))

aggregate(q.value ~ q.value_bin, gimme, range)
# q.value_bin q.value.1 q.value.2
# 1           4     0.141     0.299
# 2           5     0.302     0.482


## define motifs absolute coordinates

gimme$CRE_chr <- sub(":.*", "", gimme$sequence_name)
gimme$CRE_start <- sub("-.*", "", sub(".*:", "", gimme$sequence_name))
gimme$CRE_end <- sub(".*-", "", sub(".*:", "", gimme$sequence_name))
gimme$CRE_start <- as.integer(gimme$CRE_start)
gimme$CRE_end <- as.integer(gimme$CRE_end)

#add absolute motif locations
gimme$motif_start <- with(gimme, CRE_start + start)
gimme$motif_end <- with(gimme, CRE_start + stop)


gimme_gr <- with(gimme, GRanges( CRE_chr , IRanges( motif_start+1, motif_end )))
# using mm10 coordinates
mm10_summits$V1 <- sub("chr", "", mm10_summits$V1)
mm10_summits_gr <- with(mm10_summits, GRanges( V1 , IRanges( V2+1, V3 )))

x <- as.data.frame(findOverlaps(gimme_gr, mm10_summits_gr)) #

gimme_biochip <- cbind(gimme[x$queryHits,], mm10_summits[x$subjectHits, "mean", drop = F])# 487  17
# all unique motif instances
dim(unique(gimme_biochip[,c("CRE_chr", "motif_start", "motif_end", "strand")]))# 487   4

pdf("Nkx2.5_qval_versus_biochip_signal.pdf")
ggplot(gimme_biochip, aes(x = q.value_bin, y = mean, fill = q.value_bin)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_classic() +
  theme(legend.position="none",
        plot.title = element_text(size=11)) + scale_x_discrete(drop = FALSE)
# xlab("")
dev.off()

pdf("Nkx2.5_qval_versus_biochip_signal.1.pdf")
ggplot(gimme_biochip, aes(x = q.value_bin, y = log10(mean), fill = q.value_bin)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_classic() +
  theme(legend.position="none",
        plot.title = element_text(size=11)) + scale_x_discrete(drop = FALSE)
# xlab("")
dev.off()

pdf("Nkx2.5_qval_versus_biochip_signal_violin.pdf")
ggplot(gimme_biochip, aes(x = q.value_bin, y = log10(mean), fill = q.value_bin)) +
  geom_violin(trim = FALSE, color="black") +
  geom_boxplot(width=0.1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_classic() +
  theme(legend.position="none",
        plot.title = element_text(size=11)) + scale_x_discrete(drop = FALSE)
# xlab("")
dev.off()

table(gimme_biochip$q.value_bin)
# 1   2   3   4   5
# 0   0   0  98 389

n_motifs <- data.frame(quitile=c(1, 2, 3, 4, 5), n_motifs = c(0, 0, 0, 98, 389))
n_motifs$quitile <- factor(n_motifs$quitile ,levels = c(1, 2, 3, 4, 5))


pdf("Nkx2.5_qval_number_motifs.pdf")
ggplot(data=n_motifs, aes(x=quitile, y=n_motifs)) +
  geom_bar(stat="identity") + theme_classic() +
  scale_x_discrete(drop = FALSE)
dev.off()

write.table(x = gimme_biochip, file = "Nkx2.5_gimme_biochip.txt"
            , sep ='\t', quote = F)

# Line plots

gimme_biochip.median <- aggregate(mean ~ q.value_bin, gimme_biochip, median)
gimme_biochip.sd <- aggregate(mean ~ q.value_bin, gimme_biochip, sd)
colnames(gimme_biochip.sd)[2] <- "sd"
# sttandard error formula
std <- function(x) sd(x)/sqrt(length(x))
gimme_biochip.se <- aggregate(mean ~ q.value_bin, gimme_biochip, std)
colnames(gimme_biochip.se)[2] <- "se"

# adding levels without motifs
gimme_biochip.median <- rbind(gimme_biochip.median, data.frame(q.value_bin = c(1, 2, 3)
                                                               , mean = c(0, 0, 0)))

gimme_biochip.summary <- merge(gimme_biochip.median, gimme_biochip.sd
                               , by = "q.value_bin", all.x = T)
gimme_biochip.summary <- merge(gimme_biochip.summary, gimme_biochip.se
                               , by = "q.value_bin", all.x = T)

# mean = median of mean signal
pdf("Nkx2.5_qval_versus_biochip_signal_line.pdf")
ggplot(gimme_biochip.summary, aes(x = q.value_bin, y = mean, group=1)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) + theme_classic()
dev.off()

pdf("Nkx2.5_qval_versus_biochip_signal_line.se.pdf")
ggplot(gimme_biochip.summary, aes(x = q.value_bin, y = mean, group=1)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(0.05)) + theme_classic()
dev.off()


gimme_biochip.summary
# q.value_bin     mean       sd        se
# 1           1 0.000000       NA        NA
# 2           2 0.000000       NA        NA
# 3           3 0.000000       NA        NA
# 4           4 7.885252 6.209817 0.6272862
# 5           5 7.899775 7.189226 0.3645082
