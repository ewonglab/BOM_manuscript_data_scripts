library("Hmisc")
library(GenomicRanges)
library(ggplot2)
library(viridis)

setwd("data/mouse/e12.5")

## MOTIF: tead1
tf <- "Tead"
mm10_summits <- read.table(file = "chipseq_mouse_heart_fetal_mm10_0.bed"
                           , header = F, stringsAsFactors = F, sep = '\t')
mm10_summits$TF <- sub(".*_", "", mm10_summits$V4)
mm10_summits$id <- sub("_.*", "", mm10_summits$V4)
# read biochip signal for TF
mm10_summits <- mm10_summits[mm10_summits$TF == tf, ] # 45800     6

tead1_biochip <- read.table(file = "./biochip_mm9/Tead_combined_mean", header = T
                             , stringsAsFactors = F, sep ='\t') # 45807     4
# add mean signal to summits merging by mm19 id
mm10_summits <- merge(mm10_summits, tead1_biochip[,c("id", "mean")], by = "id")
table(is.na(mm10_summits$mean))
# FALSE
# 45800

# reading cardiomyocyte gimme
gimme <- read.table(file = "/gimme/cardiom_Gimme_vertv5.0/fimo.tsv"
                    , header = T, stringsAsFactors = F, sep = '\t')
gimme <- gimme[gimme$q.value <= 0.5, ]

# read annotation and keep only Srf motifs
gimme_annot <- read.table(file = "/useful/gimmemotifs/gimme.vertebrate.v5.0.motif2factors.txt"
                          , header = T, stringsAsFactors = F, sep ='\t')
unique(grep("tead1", gimme_annot$Factor, ignore.case = T, value=T))
# [1] "TEAD1" "Tead1"

tead1_motifs <- unique(gimme_annot[gimme_annot$Factor %in% c("TEAD1", "Tead1")
                                    , c("Factor", "Motif")])
(tead1_motifs)
#       Factor           Motif
# 5095   TEAD1 GM.5.0.TEA.0001
# 5096   Tead1 GM.5.0.TEA.0001
# 6097   TEAD1 GM.5.0.TEA.0003
# 9520   TEAD1 GM.5.0.TEA.0006
# 9964   TEAD1 GM.5.0.TEA.0008
# 10041  TEAD1 GM.5.0.TEA.0009
# 10171  TEAD1 GM.5.0.TEA.0010


gimme_annot.agr <- aggregate(Factor ~ Motif, gimme_annot, paste, collapse = ";")

gimme_annot.agr[gimme_annot.agr$Motif %in% tead1_motifs$Motif, ]
# Motif
# 1572 GM.5.0.TEA.0001
# 1574 GM.5.0.TEA.0003
# 1577 GM.5.0.TEA.0006
# 1579 GM.5.0.TEA.0008
# 1580 GM.5.0.TEA.0009
# 1581 GM.5.0.TEA.0010
# Factor
# 1572 Tead2;Tead4;TEAD3;Tead3;TEAD4;Tead4;TEAD3;TEAD4;TEAD1;Tead1;TEAD3;TEAD4;TEAD4
# 1574                                                                         TEAD1
# 1577                                                 TEAD1;TEAD3;Tead3;TEAD1;TEAD3
# 1579                                                                         TEAD1
# 1580                                                             TEAD1;TEAD1;TEAD1
# 1581                                                                         TEAD1

gimme <- gimme[gimme$motif_id %in% tead1_motifs$Motif,] #
gimme$q.value_bin <- ""
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.01, 1, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.05 & gimme$q.value > 0.01, 2, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.1 & gimme$q.value > 0.05, 3, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.3 & gimme$q.value > 0.1, 4, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.5 & gimme$q.value > 0.3, 5, gimme$q.value_bin)

table(gimme$q.value_bin)
# 4   5
# 237 510

gimme$q.value_bin <- factor(x = gimme$q.value_bin, levels = c(1, 2, 3, 4, 5))

aggregate(q.value ~ q.value_bin, gimme, range)
# q.value_bin q.value.1 q.value.2
# 1           4     0.181     0.300
# 2           5     0.301     0.476


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

gimme_biochip <- cbind(gimme[x$queryHits,], mm10_summits[x$subjectHits, "mean", drop = F])# 422  17
# all unique motif instances
dim(unique(gimme_biochip[,c("CRE_chr", "motif_start", "motif_end", "strand")]))# 422   4

pdf("tead1_qval_versus_biochip_signal_violin.pdf")
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
# 0   0   0 157 265

n_motifs <- data.frame(quitile=c(1, 2, 3, 4, 5), n_motifs = c(0, 0, 0, 157, 265))
n_motifs$quitile <- factor(n_motifs$quitile ,levels = c(1, 2, 3, 4, 5))

pdf("tead1_qval_number_motifs.pdf")
ggplot(data=n_motifs, aes(x=quitile, y=n_motifs)) +
  geom_bar(stat="identity") + theme_classic() +
  scale_x_discrete(drop = FALSE)
dev.off()

write.table(x = gimme_biochip, file = "Tead1_gimme_biochip.txt"
            , sep ='\t', quote = F)
