
library("Hmisc")
library(GenomicRanges)
library(ggplot2)
library(viridis)

## MOTIF: Mef2a
tf <- "Mef2a"
mm10_summits <- read.table(file = "chipseq_mouse_heart_fetal_mm10_0.bed"
                           , header = F, stringsAsFactors = F, sep = '\t')
mm10_summits$TF <- sub(".*_", "", mm10_summits$V4)
mm10_summits$id <- sub("_.*", "", mm10_summits$V4)
# read biochip signal for TF
mm10_summits <- mm10_summits[mm10_summits$TF == tf, ] # 15931     6

Mef2a_biochip <- read.table(file = "./biochip_mm9/Mef2a_combined_mean", header = T
                             , stringsAsFactors = F, sep ='\t') #
# add mean signal to summits merging by mm19 id
mm10_summits <- merge(mm10_summits, Mef2a_biochip[,c("id", "mean")], by = "id")
table(is.na(mm10_summits$mean))
# FALSE
# 15931

# reading cardiomyocyte gimme
gimme <- read.table(file = "/gimme/cardiom_Gimme_vertv5.0/fimo.tsv"
                    , header = T, stringsAsFactors = F, sep = '\t')
gimme <- gimme[gimme$q.value <= 0.5, ]

# read annotation and keep only Srf motifs
gimme_annot <- read.table(file = "/useful/gimmemotifs/gimme.vertebrate.v5.0.motif2factors.txt"
                          , header = T, stringsAsFactors = F, sep ='\t')
unique(grep("Mef2a", gimme_annot$Factor, ignore.case = T, value=T))
# [1] "MEF2A" "Mef2a"

Mef2a_motifs <- unique(gimme_annot[gimme_annot$Factor %in% c("MEF2A", "Mef2a")
                                    , c("Factor", "Motif")])
(Mef2a_motifs)
#       Factor                Motif
# 1521   MEF2A GM.5.0.MADS_box.0002
# 1524   Mef2a GM.5.0.MADS_box.0002
# 2337   MEF2A     GM.5.0.bZIP.0013
# 2474   MEF2A     GM.5.0.Runt.0003
# 3194   Mef2a GM.5.0.MADS_box.0003
# 3198   MEF2A GM.5.0.MADS_box.0003
# 5514   MEF2A GM.5.0.MADS_box.0005
# 5773   MEF2A GM.5.0.MADS_box.0007
# 8007   MEF2A GM.5.0.MADS_box.0012
# 8011   MEF2A GM.5.0.MADS_box.0013
# 8014   Mef2a GM.5.0.MADS_box.0013
# 8179   MEF2A GM.5.0.MADS_box.0014
# 8182   Mef2a GM.5.0.MADS_box.0014
# 9354   MEF2A GM.5.0.MADS_box.0020
# 9357   Mef2a GM.5.0.MADS_box.0020
# 10082  Mef2a GM.5.0.MADS_box.0023


gimme_annot.agr <- aggregate(Factor ~ Motif, gimme_annot, paste, collapse = ";")

gimme_annot.agr[gimme_annot.agr$Motif %in% Mef2a_motifs$Motif, ]
# Motif
# 173      GM.5.0.bZIP.0013
# 1065 GM.5.0.MADS_box.0002
# 1066 GM.5.0.MADS_box.0003
# 1068 GM.5.0.MADS_box.0005
# 1070 GM.5.0.MADS_box.0007
# 1075 GM.5.0.MADS_box.0012
# 1076 GM.5.0.MADS_box.0013
# 1077 GM.5.0.MADS_box.0014
# 1083 GM.5.0.MADS_box.0020
# 1086 GM.5.0.MADS_box.0023
# 1447     GM.5.0.Runt.0003
# Factor
# 173  FOS;BATF;CEBPB;EP300;FOSL1;Fra2;HMGN3;JDP2;JDP2;Jdp2;FOS;Fosb;Fos;FOSL1;JUN;JUNB;JUND;JUND;SMARCC1;Smarcc2;Smarcc1;FOS;Fosb;Fos;FOSL2;JUN;FOSL1;JUNB;ATF3;JDP2;NFE2;Nfe2;Jdp2;FOS;JUN;JDP2;NFE2;FOSL1;FOSL2;JUNB;FOSB;JUND;MEF2A;MEF2C;NFE2;NR3C1;PRDM1;SMARCB1;TAL1;TCF7L2;TRIM28;JUN;JUNB;JUND
# 1065                                                                                                           AC002126.6;MEF2A;Mef2d;Mef2c;Mef2a;AC002126.6;Mef2d;Mef2c;Mef2a;MEF2A;AC002126.6;MEF2A;Mef2d;Mef2c;Mef2a;MEF2B;MEF2D;MEF2C;MEF2A;MEF2B;MEF2D;MEF2B;MEF2B;MEF2D;MEF2A;MEF2A;MEF2A;MYEF2
# 1066                                                                                                                                                                Mef2a;Mef2c;Mef2d;AC002126.6;MEF2A;Mef2d;Mef2c;Mef2a;MEF2C;AC002126.6;MEF2A;Mef2d;Mef2c;Mef2a;MEF2C;MEF2A;MEF2C;MEF2B;MEF2C;MEF2D
# 1068                                                                                                                                                                                                                                                                     BCL11A;MEF2A;MEF2C;PAX5;RXRA
# 1070                                                                                                                                                                                                                                                                                      MEF2A;MEF2C
# 1075                                                                                                                                                                                                                                                                                      MEF2A;MYEF2
# 1076                                                                                                                                                                                                                                                         AC002126.6;MEF2A;Mef2d;Mef2c;Mef2a;MYEF2
# 1077                                                                                                                                                                                                                                                         AC002126.6;MEF2A;Mef2d;Mef2c;Mef2a;MYEF2
# 1083                                                                                                                                                                                                                                                               AC002126.6;MEF2A;Mef2d;Mef2c;Mef2a
# 1086                                                                                                                                                                                                                                                               AC002126.6;MEF2D;Mef2d;Mef2c;Mef2a
# 1447                                                                                                                                                                                                          RUNX1;ENSG00000250096;Runx1;Runx2;Runx3;MEF2A;MEF2C;RUNX1;Runx1;RUNX2;Runx2;RUNX3;Runx3

gimme <- gimme[gimme$motif_id %in% Mef2a_motifs$Motif,] #
gimme$q.value_bin <- ""
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.01, 1, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.05 & gimme$q.value > 0.01, 2, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.1 & gimme$q.value > 0.05, 3, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.3 & gimme$q.value > 0.1, 4, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.5 & gimme$q.value > 0.3, 5, gimme$q.value_bin)

table(gimme$q.value_bin)
# 1    2    3    4    5
# 216  755  511 1047  885

gimme$q.value_bin <- factor(x = gimme$q.value_bin, levels = c(1, 2, 3, 4, 5))

aggregate(q.value ~ q.value_bin, gimme, range)
#   q.value_bin q.value.1 q.value.2
# 1           1  0.000192  0.009940
# 2           2  0.010200  0.050000
# 3           3  0.050100  0.098500
# 4           4  0.101000  0.299000
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

gimme_biochip <- cbind(gimme[x$queryHits,], mm10_summits[x$subjectHits, "mean", drop = F])# 908  17
# all unique motif instances
dim(unique(gimme_biochip[,c("CRE_chr", "motif_start", "motif_end", "strand")]))# 807   4

# In this case I have motifs with different ids mapping to the same coordinate and strand:
# example:
# motif_id motif_alt_id          sequence_name start stop
# 21376 GM.5.0.MADS_box.0003           NA 14:111887671-111888192    97  106
# 32139 GM.5.0.MADS_box.0002           NA 14:111887671-111888192    97  106
# strand   score  p.value q.value matched_sequence q.value_bin CRE_chr
# 21376      - 15.7059 1.64e-06  0.0621       CTAAAAATAG           3      14
# 32139      - 17.1928 3.28e-06  0.0887       CTAAAAATAG           3      14
# CRE_start   CRE_end motif_start motif_end    mean
# 21376 111887671 111888192   111887768 111887777 5.61962
# 32139 111887671 111888192   111887768 111887777 5.61962

# I think I should keep unique motif instances based on the location
gimme_biochip$location_id <- with(gimme_biochip, paste(CRE_chr, motif_start, motif_end, strand, sep = "_"))
gimme_biochip <- gimme_biochip[!duplicated(gimme_biochip$location_id),]# 807  18

pdf("Mef2a_qval_versus_biochip_signal_violin.pdf")
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
# 12  96 128 352 219

n_motifs <- data.frame(quitile=c(1, 2, 3, 4, 5), n_motifs = c(12, 96, 128, 352, 219))
n_motifs$quitile <- factor(n_motifs$quitile ,levels = c(1, 2, 3, 4, 5))


pdf("Mef2a_qval_number_motifs.pdf")
ggplot(data=n_motifs, aes(x=quitile, y=n_motifs)) +
  geom_bar(stat="identity") + theme_classic() +
  scale_x_discrete(drop = FALSE)
dev.off()

gimme_biochip$location_id <- NULL
write.table(x = gimme_biochip, file = "Mef2a_gimme_biochip.txt"
            , sep ='\t', quote = F)
