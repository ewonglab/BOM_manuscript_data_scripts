library("Hmisc")
library(GenomicRanges)
library(ggplot2)
library(viridis)

setwd("data/mouse/e12.5")

## MOTIF: Tbx5
tf <- "Tbx5"
mm10_summits <- read.table(file = "chipseq_mouse_heart_fetal_mm10_0.bed"
                           , header = F, stringsAsFactors = F, sep = '\t')
mm10_summits$TF <- sub(".*_", "", mm10_summits$V4)
mm10_summits$id <- sub("_.*", "", mm10_summits$V4)
# read biochip signal for TF
mm10_summits <- mm10_summits[mm10_summits$TF == tf, ] # 30785     6

Tbx5_biochip <- read.table(file = "./biochip_mm9/Tbx5_combined_mean", header = T
                             , stringsAsFactors = F, sep ='\t') # 30789     4
# add mean signal to summits merging by mm19 id
mm10_summits <- merge(mm10_summits, Tbx5_biochip[,c("id", "mean")], by = "id")
table(is.na(mm10_summits$mean))
# FALSE
# 30785

# reading cardiomyocyte gimme
gimme <- read.table(file = "/g/data/zk16/cc3704/mouse_data/gottgens_scATAC/fimo/gimme/cardiom_Gimme_vertv5.0/fimo.tsv"
                    , header = T, stringsAsFactors = F, sep = '\t')
gimme <- gimme[gimme$q.value <= 0.5, ]

# read annotation and keep only Srf motifs
gimme_annot <- read.table(file = "/g/data/zk16/useful/gimmemotifs/gimme.vertebrate.v5.0.motif2factors.txt"
                          , header = T, stringsAsFactors = F, sep ='\t')
unique(grep("Tbx5", gimme_annot$Factor, ignore.case = T, value=T))
# "Tbx5" "TBX5"

Tbx5_motifs <- unique(gimme_annot[gimme_annot$Factor %in% c("Tbx5", "TBX5")
                                    , c("Factor", "Motif")])
(Tbx5_motifs)
#       Factor             Motif
# 2109    Tbx5 GM.5.0.T-box.0003
# 2110    TBX5 GM.5.0.T-box.0003
# 4085    Tbx5 GM.5.0.T-box.0005
# 4100    TBX5 GM.5.0.T-box.0005
# 4633    TBX5 GM.5.0.T-box.0008
# 6710    TBX5 GM.5.0.T-box.0014
# 10215   TBX5 GM.5.0.T-box.0026

gimme_annot.agr <- aggregate(Factor ~ Motif, gimme_annot, paste, collapse = ";")

gimme_annot.agr[gimme_annot.agr$Motif %in% Tbx5_motifs$Motif, ]
# Motif
# 1542 GM.5.0.T-box.0003
# 1544 GM.5.0.T-box.0005
# 1547 GM.5.0.T-box.0008
# 1553 GM.5.0.T-box.0014
# 1565 GM.5.0.T-box.0026
# Factor
# 1542                                                                                                                                                              Tbx4;Tbx5;TBX5
# 1544 Tbx5;Tbx20;Tbx21;Eomes;MGA;Mga;TBX10;TBX1;Tbx10;TBX15;Tbx15;Tbx22;Tbx18;TBX2;TBX4;TBX5;TBX5;TBX2;MGA;TBX15;TBX1;TBX4;TBX5;MGA;TBX1;TBX15;TBX1;TBX2;TBX2;TBX4;TBX4;TBX5;TBX5
# 1547                                                                                                                                                              Tbx19;T;T;TBX5
# 1553                                                                                                                                   MGA;Mga;TBX2;TBX4;TBX5;MGA;TBX2;TBX4;TBX5
# 1565                                                                                                                                                                        TBX5


gimme <- gimme[gimme$motif_id %in% Tbx5_motifs$Motif,] #
gimme$q.value_bin <- ""
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.01, 1, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.05 & gimme$q.value > 0.01, 2, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.1 & gimme$q.value > 0.05, 3, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.3 & gimme$q.value > 0.1, 4, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.5 & gimme$q.value > 0.3, 5, gimme$q.value_bin)

table(gimme$q.value_bin)
# 4   5
# 6 300

gimme$q.value_bin <- factor(x = gimme$q.value_bin, levels = c(1, 2, 3, 4, 5))

aggregate(q.value ~ q.value_bin, gimme, range)
# q.value_bin q.value.1 q.value.2
# 1           4     0.177     0.227
# 2           5     0.483     0.499


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

gimme_biochip <- cbind(gimme[x$queryHits,], mm10_summits[x$subjectHits, "mean", drop = F])# 148  17
# all unique motif instances
dim(unique(gimme_biochip[,c("CRE_chr", "motif_start", "motif_end", "strand")]))#  148   4

# I think I should keep unique motif instances based on the location
gimme_biochip$location_id <- with(gimme_biochip, paste(CRE_chr, motif_start, motif_end, strand, sep = "_"))
length(gimme_biochip$location_id[duplicated(gimme_biochip$location_id)])
# [1] 0

pdf("Tbx5_qval_versus_biochip_signal_violin.pdf")
ggplot(gimme_biochip, aes(x = q.value_bin, y = log10(mean), fill = q.value_bin)) +
  geom_violin(trim = FALSE, color="black") +
  geom_boxplot(width=0.1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_classic() +
  theme(legend.position="none",
        plot.title = element_text(size=11)) + scale_x_discrete(drop = FALSE)
dev.off()

table(gimme_biochip$q.value_bin)
# 1   2   3   4   5
# 0   0   0   4 144

n_motifs <- data.frame(quitile=c(1, 2, 3, 4, 5), n_motifs = c(0, 0, 0, 4, 144))
n_motifs$quitile <- factor(n_motifs$quitile ,levels = c(1, 2, 3, 4, 5))


pdf("Tbx5_qval_number_motifs.pdf")
ggplot(data=n_motifs, aes(x=quitile, y=n_motifs)) +
  geom_bar(stat="identity") + theme_classic() +
  scale_x_discrete(drop = FALSE)
dev.off()

gimme_biochip$location_id <- NULL
write.table(x = gimme_biochip, file = "Tbx5_gimme_biochip.txt"
            , sep ='\t', quote = F)
