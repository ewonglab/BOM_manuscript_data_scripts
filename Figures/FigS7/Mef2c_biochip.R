Mean biochip signal across q-value thresholds: Mef2c motifs

library("Hmisc", lib = "/g/data/zk16/software/Rpackages_paola/R_4.0.0")
library(GenomicRanges)
library(ggplot2)
library(viridis)

setwd("/g/data/zk16/cc3704/mouse_enh_grammar/data/mouse/e12.5")

## MOTIF: Mef2c
tf <- "Mef2c"
mm10_summits <- read.table(file = "chipseq_mouse_heart_fetal_mm10_0.bed"
                           , header = F, stringsAsFactors = F, sep = '\t')
mm10_summits$TF <- sub(".*_", "", mm10_summits$V4)
mm10_summits$id <- sub("_.*", "", mm10_summits$V4)
# read biochip signal for TF
mm10_summits <- mm10_summits[mm10_summits$TF == tf, ] # 30174     6

Mef2c_biochip <- read.table(file = "./biochip_mm9/Mef2c_combined_mean", header = T
                             , stringsAsFactors = F, sep ='\t') # 30178     4
# add mean signal to summits merging by mm19 id
mm10_summits <- merge(mm10_summits, Mef2c_biochip[,c("id", "mean")], by = "id")
table(is.na(mm10_summits$mean))
# FALSE
# 30174

# reading cardiomyocyte gimme
gimme <- read.table(file = "/g/data/zk16/cc3704/mouse_data/gottgens_scATAC/fimo/gimme/cardiom_Gimme_vertv5.0/fimo.tsv"
                    , header = T, stringsAsFactors = F, sep = '\t')
gimme <- gimme[gimme$q.value <= 0.5, ]

# read annotation and keep only Srf motifs
gimme_annot <- read.table(file = "/g/data/zk16/useful/gimmemotifs/gimme.vertebrate.v5.0.motif2factors.txt"
                          , header = T, stringsAsFactors = F, sep ='\t')
unique(grep("Mef2c", gimme_annot$Factor, ignore.case = T, value=T))
# "Mef2c" "MEF2C"

Mef2c_motifs <- unique(gimme_annot[gimme_annot$Factor %in% c("Mef2c", "MEF2C")
                                    , c("Factor", "Motif")])
(Mef2c_motifs)
#       Factor                Motif
# 1523   Mef2c GM.5.0.MADS_box.0002
# 1537   MEF2C GM.5.0.MADS_box.0002
# 2338   MEF2C     GM.5.0.bZIP.0013
# 2475   MEF2C     GM.5.0.Runt.0003
# 3195   Mef2c GM.5.0.MADS_box.0003
# 3202   MEF2C GM.5.0.MADS_box.0003
# 5515   MEF2C GM.5.0.MADS_box.0005
# 5629   MEF2C GM.5.0.MADS_box.0006
# 5774   MEF2C GM.5.0.MADS_box.0007
# 8013   Mef2c GM.5.0.MADS_box.0013
# 8181   Mef2c GM.5.0.MADS_box.0014
# 9337   MEF2C GM.5.0.MADS_box.0019
# 9356   Mef2c GM.5.0.MADS_box.0020
# 10081  Mef2c GM.5.0.MADS_box.0023


gimme_annot.agr <- aggregate(Factor ~ Motif, gimme_annot, paste, collapse = ";")

gimme_annot.agr[gimme_annot.agr$Motif %in% Mef2c_motifs$Motif, ]
# Motif
# 173      GM.5.0.bZIP.0013
# 1065 GM.5.0.MADS_box.0002
# 1066 GM.5.0.MADS_box.0003
# 1068 GM.5.0.MADS_box.0005
# 1069 GM.5.0.MADS_box.0006
# 1070 GM.5.0.MADS_box.0007
# 1076 GM.5.0.MADS_box.0013
# 1077 GM.5.0.MADS_box.0014
# 1082 GM.5.0.MADS_box.0019
# 1083 GM.5.0.MADS_box.0020
# 1086 GM.5.0.MADS_box.0023
# 1447     GM.5.0.Runt.0003
# Factor
# 173  FOS;BATF;CEBPB;EP300;FOSL1;Fra2;HMGN3;JDP2;JDP2;Jdp2;FOS;Fosb;Fos;FOSL1;JUN;JUNB;JUND;JUND;SMARCC1;Smarcc2;Smarcc1;FOS;Fosb;Fos;FOSL2;JUN;FOSL1;JUNB;ATF3;JDP2;NFE2;Nfe2;Jdp2;FOS;JUN;JDP2;NFE2;FOSL1;FOSL2;JUNB;FOSB;JUND;MEF2A;MEF2C;NFE2;NR3C1;PRDM1;SMARCB1;TAL1;TCF7L2;TRIM28;JUN;JUNB;JUND
# 1065                                                                                                           AC002126.6;MEF2A;Mef2d;Mef2c;Mef2a;AC002126.6;Mef2d;Mef2c;Mef2a;MEF2A;AC002126.6;MEF2A;Mef2d;Mef2c;Mef2a;MEF2B;MEF2D;MEF2C;MEF2A;MEF2B;MEF2D;MEF2B;MEF2B;MEF2D;MEF2A;MEF2A;MEF2A;MYEF2
# 1066                                                                                                                                                                Mef2a;Mef2c;Mef2d;AC002126.6;MEF2A;Mef2d;Mef2c;Mef2a;MEF2C;AC002126.6;MEF2A;Mef2d;Mef2c;Mef2a;MEF2C;MEF2A;MEF2C;MEF2B;MEF2C;MEF2D
# 1068                                                                                                                                                                                                                                                                     BCL11A;MEF2A;MEF2C;PAX5;RXRA
# 1069                                                                                                                                                                                                                                                                                            MEF2C
# 1070                                                                                                                                                                                                                                                                                      MEF2A;MEF2C
# 1076                                                                                                                                                                                                                                                         AC002126.6;MEF2A;Mef2d;Mef2c;Mef2a;MYEF2
# 1077                                                                                                                                                                                                                                                         AC002126.6;MEF2A;Mef2d;Mef2c;Mef2a;MYEF2
# 1082                                                                                                                                                                                                                                                                                            MEF2C
# 1083                                                                                                                                                                                                                                                               AC002126.6;MEF2A;Mef2d;Mef2c;Mef2a
# 1086                                                                                                                                                                                                                                                               AC002126.6;MEF2D;Mef2d;Mef2c;Mef2a
# 1447                                                                                                                                                                                                          RUNX1;ENSG00000250096;Runx1;Runx2;Runx3;MEF2A;MEF2C;RUNX1;Runx1;RUNX2;Runx2;RUNX3;Runx3
#

gimme <- gimme[gimme$motif_id %in% Mef2c_motifs$Motif,] #
gimme$q.value_bin <- ""
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.01, 1, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.05 & gimme$q.value > 0.01, 2, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.1 & gimme$q.value > 0.05, 3, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.3 & gimme$q.value > 0.1, 4, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.5 & gimme$q.value > 0.3, 5, gimme$q.value_bin)

table(gimme$q.value_bin)
# 1    2    3    4    5
# 216  755  496 1024  933

gimme$q.value_bin <- factor(x = gimme$q.value_bin, levels = c(1, 2, 3, 4, 5))

aggregate(q.value ~ q.value_bin, gimme, range)

#   q.value_bin q.value.1 q.value.2
# 1           1  0.000192  0.009940
# 2           2  0.010200  0.050000
# 3           3  0.050100  0.098500
# 4           4  0.101000  0.298000
# 5           5  0.301000  0.498000

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

gimme_biochip <- cbind(gimme[x$queryHits,], mm10_summits[x$subjectHits, "mean", drop = F])# 937  17
# all unique motif instances
dim(unique(gimme_biochip[,c("CRE_chr", "motif_start", "motif_end", "strand")]))# 857   4

# I think I should keep unique motif instances based on the location
gimme_biochip$location_id <- with(gimme_biochip, paste(CRE_chr, motif_start, motif_end, strand, sep = "_"))
length(gimme_biochip$location_id[duplicated(gimme_biochip$location_id)])
# [1] 80
gimme_biochip <- gimme_biochip[!duplicated(gimme_biochip$location_id),]# 857  18

pdf("Mef2c_qval_versus_biochip_signal_violin.pdf")
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
# 31 129 132 312 253

n_motifs <- data.frame(quitile=c(1, 2, 3, 4, 5), n_motifs = c(31, 129, 132, 312, 253))
n_motifs$quitile <- factor(n_motifs$quitile ,levels = c(1, 2, 3, 4, 5))


pdf("Mef2c_qval_number_motifs.pdf")
ggplot(data=n_motifs, aes(x=quitile, y=n_motifs)) +
  geom_bar(stat="identity") + theme_classic() +
  scale_x_discrete(drop = FALSE)
dev.off()

gimme_biochip$location_id <- NULL
write.table(x = gimme_biochip, file = "Mef2c_gimme_biochip.txt"
            , sep ='\t', quote = F)
