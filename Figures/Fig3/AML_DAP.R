
# Differentially accessible peaks AML and normal condition
library(DiffBind)
library(tidyverse)

# SRR2920581 # SU501-pHSC
# SRR2920580 # SU501-Leuk-
# SRR2920576 # SU484-pHSC
# SRR2920575 # SU484-Leuk-
# SRR2920595 # SU654-pHSC-
# SRR2920593 # SU654-Leuk

setwd("/g/data/zk16/cc3704/replication_timing/human/AML/DAP")

np_path <- "human/AML/narrowPeak/"
bam_path <- "human/AML/align/"
samples <- data.frame(SampleID = c(paste0("AML_", 1:3), paste0("normal_", 1:3))
                      , Condition = c(rep("AML", 3), rep("normal", 3))
                      , Replicate = rep(1:3, 2)
                      , bamReads = paste0(bam_path
                                          , c("SRR2920575_1.fastq.bam" # # SRR2920575 # SU484-Leuk-
                                     , "SRR2920580_1.fastq.bam" # # SRR2920580 # SU501-Leuk-
                                     , "SRR2920593_1.fastq.bam" # SRR2920593 # SU654-Leuk
                                     , "SRR2920576_1.fastq.bam" # # SRR2920576 # SU484-pHSC
                                     , "SRR2920581_1.fastq.bam" # # SRR2920581 # SU501-pHSC
                                     , "SRR2920595_1.fastq.bam")) # # SRR2920595 # SU654-pHSC-
                      , Peaks = paste0(np_path, c("SU484_blast_peaks.narrowPeak" # blast = keukemia
                                                  , "SU501_blast_peaks.narrowPeak"
                                                  , "SU654_blast_peaks.narrowPeak"
                                                  , "SU484_pHSC_peaks.narrowPeak" # pre-leukemia/normal
                                                  , "SU501_pHSC_peaks.narrowPeak"
                                                  , "SU654_pHSC_peaks.narrowPeak"))
                      , PeakCaller = "narrow")

dba_object <- dba(sampleSheet = samples)

dba_object <- dba.count(dba_object, summits = 200, bUseSummarizeOverlaps = TRUE)
saveRDS(object = dba_object, file = "dba_object.Rds")
dba_object <- readRDS(file = "dba_object.Rds")
dba_object
# 6 Samples, 104146 sites in matrix:
#         ID Condition Replicate    Reads FRiP
# 1    AML_1       AML         1 11611699 0.42
# 2    AML_2       AML         2 20616394 0.46
# 3    AML_3       AML         3 26161826 0.47
# 4 normal_1    normal         1 16372202 0.54
# 5 normal_2    normal         2 12239175 0.54
# 6 normal_3    normal         3 22614388 0.54

# dba_object$
#  dba_object$peaks        dba_object$samples      dba_object$merged       dba_object$masks        dba_object$minCount
# dba_object$class        dba_object$called       dba_object$totalMerged  dba_object$SN           dba_object$summits
# dba_object$chrmap       dba_object$score        dba_object$attributes   dba_object$maxFilter    
# dba_object$config       dba_object$binding      dba_object$minOverlap   dba_object$filterFun

set.seed(123)
dba_object <- dba.contrast(dba_object, categories = DBA_CONDITION, minMembers = 2)
DBA_ALL_METHODS
# [1] "edgeRGLM" "DESeq2"
dba_object <- dba.analyze(dba_object, method = DBA_ALL_METHODS)
# Applying Blacklist/Greylists...
# Genome detected: Hsapiens.1000genomes.hs37d5
# Applying blacklist...
# No blacklist found for BSgenome.Hsapiens.1000genomes.hs37d5
# No intervals removed.
# Normalize edgeR with defaults...
# Normalize DESeq2 with defaults...
# Analyzing...
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates

# bContrasts: logical indicating whether peaksets or contrast attributes
# are to be retrieved. ‘TRUE’ retrieves a dataframe of contrast
# information instead of peakset attributes.  If no contrasts
# are set, returns possible contrasts. See ‘dba.contrast’.

dba.show(dba_object, bContrasts = TRUE)
#       Factor  Group Samples Group2 Samples2 DB.edgeR DB.DESeq2
# 1 Condition normal       3    AML        3    26051     45569

# default values stored in object
dba_object$config$th
# [1] 0.05
dba_object$config$bUsePval
# FALSE

res_deseq <- dba.report(dba_object, method = DBA_DESEQ2, contrast = 1)
length(res_deseq)
res_edger <- dba.report(dba_object, method=DBA_EDGER, contrast = 1)
length(res_edger)

saveRDS(object = res_deseq, file = "diffBind_DAP_deseq2.rds")
saveRDS(object = res_edger, file = "diffBind_DAP_edger.rds")

res_deseq.df <- as.data.frame(res_deseq)
res_edger.df <- as.data.frame(res_edger)

library(GenomicRanges)

res_deseq.df$id <- with(res_deseq.df, paste(seqnames, paste(start, end, sep = "-"), sep = ":"))
res_edger.df$id <- with(res_edger.df, paste(seqnames, paste(start, end, sep = "-"), sep = ":"))

x <- as.data.frame(findOverlaps(res_deseq, res_edger))

# filter by fold value
res_deseq.fdr <- res_deseq.df[abs(res_deseq.df$Fold) > 1, ]# 41690    12
res_edger.fdr <- res_edger.df[abs(res_edger.df$Fold) > 1, ]# 25800    12

keep <- intersect(res_deseq.fdr$id, res_edger.fdr$id)
deseq_conf <- res_deseq.fdr[res_deseq.fdr$id %in% keep, ]
edger_conf <- res_edger.fdr[res_edger.fdr$id %in% keep, ]

# making sure signifficant peaks have the same direction of enrichment
# CONSISTENT DIRECTION BETWEEN TWO METHODS
length(setdiff(deseq_conf[deseq_conf$Fold > 0, "id"], edger_conf[edger_conf$Fold > 0, "id"]))
# [1] 0
length(setdiff(edger_conf[edger_conf$Fold > 0, "id"], deseq_conf[deseq_conf$Fold > 0, "id"]))
# [1] 0
length(intersect(edger_conf[edger_conf$Fold > 0, "id"], deseq_conf[deseq_conf$Fold > 0, "id"]))
# [1] 11052

length(setdiff(deseq_conf[deseq_conf$Fold < 0, "id"], edger_conf[edger_conf$Fold < 0, "id"]))
# [1] 0
length(setdiff(edger_conf[edger_conf$Fold < 0, "id"], deseq_conf[deseq_conf$Fold < 0, "id"]))
# [1] 0
length(intersect(edger_conf[edger_conf$Fold < 0, "id"], deseq_conf[deseq_conf$Fold < 0, "id"]))
# [1] 14644

# keep only main chromosomes
table(deseq_conf$seqnames)
# 1         10         11         12         13         14         15
# 2261       1233       1200       1362        821        783        742
# 16         17         18         19          2         20         21
# 696        813        657        491       2249        708        345
# 22          3          4          5          6          7          8
# 318       1930       1356       1569       1613       1306       1393
# 9 GL000195.1 GL000202.1 GL000205.1 GL000209.1 GL000220.1          X
# 1111          1          1          3          1          2        730
# Y
# 1

deseq_conf <- deseq_conf[!deseq_conf$seqnames %in%
                           c("GL000195.1", "GL000202.1", "GL000205.1", "GL000209.1", "GL000220.1"),]

unique(deseq_conf$seqnames)
# [1] 4  10 20 6  9  5  12 21 8  1  16 13 11 2  X  15 7  18 3  14 17 22 19 Y
deseq_conf$seqnames <- paste0("chr", deseq_conf$seqnames)
dap_normal <- deseq_conf[deseq_conf$Fold > 0,c("seqnames", "start", "end", "id")]# 11046     4
dap_AML <- deseq_conf[deseq_conf$Fold < 0,c("seqnames", "start", "end", "id")]# 14642     4

write.table(x = dap_normal, file = "dap_normal.bed", col.names = F, row.names = F, sep = '\t', quote = F)
write.table(x = dap_AML, file = "dap_AML.bed", col.names = F, row.names = F, sep = '\t', quote = F)
