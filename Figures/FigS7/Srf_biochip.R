library("Hmisc")
library(GenomicRanges)

setwd("data/mouse/e12.5")

## MOTIF: SRF
tf <- "Srf"
mm10_summits <- read.table(file = "chipseq_mouse_heart_fetal_mm10_0.bed"
                           , header = F, stringsAsFactors = F, sep = '\t')
mm10_summits$TF <- sub(".*_", "", mm10_summits$V4)
mm10_summits$id <- sub("_.*", "", mm10_summits$V4)
# read biochip signal for TF
mm10_summits <- mm10_summits[mm10_summits$TF == tf, ] # 18519     6

srf_biochip <- read.table(file = "./biochip_mm9/Srf_combined_mean", header = T
                          , stringsAsFactors = F, sep ='\t')
# add mean signal to summits merging by mm19 id
mm10_summits <- merge(mm10_summits, srf_biochip[,c("id", "mean")], by = "id")
table(is.na(mm10_summits$mean))
# FALSE
# 18519

# reading cardiomyocyte gimme
gimme <- read.table(file = "/gimme/cardiom_Gimme_vertv5.0/fimo.tsv"
                    , header = T, stringsAsFactors = F, sep = '\t')
gimme <- gimme[gimme$q.value <= 0.5, ]

# read annotation and keep only Srf motifs
gimme_annot <- read.table(file = "/useful/gimmemotifs/gimme.vertebrate.v5.0.motif2factors.txt"
                          , header = T, stringsAsFactors = F, sep ='\t')

srf_motifs <- unique(gimme_annot[gimme_annot$Factor %in% c("SRF", "Srf"), c("Factor", "Motif")])
# Factor                Motif
# 26       SRF    GM.5.0.Mixed.0001
# 528      SRF GM.5.0.MADS_box.0001
# 533      Srf GM.5.0.MADS_box.0001
# 5497     SRF GM.5.0.MADS_box.0004
# 6129     SRF GM.5.0.MADS_box.0008
# 6401     SRF GM.5.0.MADS_box.0009
# 8346     SRF GM.5.0.MADS_box.0015
# 8420     SRF GM.5.0.MADS_box.0016
# 8962     SRF GM.5.0.MADS_box.0018
# 10012    SRF GM.5.0.MADS_box.0022
# 10014    Srf GM.5.0.MADS_box.0022

gimme_annot.agr <- aggregate(Factor ~ Motif, gimme_annot, paste, collapse = ";")

gimme_annot.agr[gimme_annot.agr$Motif %in% srf_motifs$Motif, ]
#                     Motif                      Factor
# 1064 GM.5.0.MADS_box.0001     SRF;SRF;SRF;SRF;SRF;Srf
# 1067 GM.5.0.MADS_box.0004                         SRF
# 1071 GM.5.0.MADS_box.0008                         SRF
# 1072 GM.5.0.MADS_box.0009 SRF;SRF;SRF;SRF;SRF;SRF;SRF
# 1078 GM.5.0.MADS_box.0015                         SRF
# 1079 GM.5.0.MADS_box.0016                         SRF
# 1081 GM.5.0.MADS_box.0018                         SRF
# 1085 GM.5.0.MADS_box.0022             SRF;SRF;Srf;SRF
# 1090    GM.5.0.Mixed.0001                    EGR1;SRF

gimme <- gimme[gimme$motif_id %in% srf_motifs$Motif,] # [1] 1546   10
gimme$q.value_bin <- ""
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.01, 1, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.05 & gimme$q.value > 0.01, 2, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.1 & gimme$q.value > 0.05, 3, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.3 & gimme$q.value > 0.1, 4, gimme$q.value_bin)
gimme$q.value_bin <- ifelse(gimme$q.value <= 0.5 & gimme$q.value > 0.3, 5, gimme$q.value_bin)

table(gimme$q.value_bin)
#  1   2   3   4   5
# 90  73 119 734 530

gimme$q.value_bin <- as.factor(gimme$q.value_bin)

aggregate(q.value ~ q.value_bin, gimme, range)
#   q.value_bin q.value.1 q.value.2
# 1           1   0.00702   0.00972
# 2           2   0.01090   0.04400
# 3           3   0.05060   0.09700
# 4           4   0.10100   0.30000
# 5           5   0.30200   0.47300

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

x <- as.data.frame(findOverlaps(gimme_gr, mm10_summits_gr)) #[1] 698   2

gimme_biochip <- cbind(gimme[x$queryHits,], mm10_summits[x$subjectHits, "mean", drop = F])# 698  17
# all unique motif instances
dim(unique(gimme_biochip[,c("CRE_chr", "motif_start", "motif_end", "strand")]))# [1] 698   4
gimme_biochip$location_id <- with(gimme_biochip, paste(CRE_chr, motif_start, motif_end, strand, sep = "_"))

# No duplicated motif locations (exactly same location and strand but different id)

# length(gimme_biochip$location_id[duplicated(gimme_biochip$location_id)])
# [1] 0


library(ggplot2)
library(viridis)

pdf("Srf_qval_versus_biochip_signal_violin.pdf")
ggplot(gimme_biochip, aes(x = q.value_bin, y = log10(mean), fill = q.value_bin)) +
  geom_violin(trim = FALSE, color="black") +
  geom_boxplot(width=0.1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_classic() +
  theme(legend.position="none",
        plot.title = element_text(size=11))
# xlab("")
dev.off()

table(gimme_biochip$q.value_bin)
# 1   2   3   4   5
# 10  24  71 403 190

n_motifs <- data.frame(quitile=c(1, 2, 3, 4, 5), n_motifs = c(10, 24, 71, 403, 190))
n_motifs$quitile <- factor(n_motifs$quitile ,levels = c(1, 2, 3, 4, 5))

pdf("Srf_qval_number_motifs.pdf")
ggplot(data=n_motifs, aes(x=quitile, y=n_motifs)) +
  geom_bar(stat="identity") + theme_classic() +
  scale_x_discrete(drop = FALSE)
dev.off()

