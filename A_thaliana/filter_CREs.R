library(data.table)
library(GenomicRanges)
library(GenomicFeatures)

# Function used to trim CREs
trim_region <- function(x, n){
  x$width <- with(x, end - start)
  x$centre <- with(x, start + round(width/2))
  x$start <- round(x$centre - n/2)
  x$end <- round(x$centre + n/2)
  print(unique(with(x, end-start)))
  return(x)
}

peaks <- read.csv(file = "Supplementary_Table_1_MarkersBy.csv", header =T, stringsAsFactors=F)

# based on Figure 1d
# cortex = CLUSTER 3
# c/e precursor = CLUSTER 7
# endodermis 3 = CLUSTER 8
# endodermis 2 = CLUSTER 10
# endodermis 1 = CLUSTER 4
# epidermis = CLUSTER 0
# stele 1 (xylem) = CLUSTER 1
# stele 2 (peri) = CLUSTER 2
# stele (phloem) = CLUSTER 11

cluster_annot <- data.frame(cluster = c(3, 7, 8, 10, 4, 0, 1, 2, 11)
                            , cellType = c("cortex", "c.e_precursor"
                                           , "endodermis_3", "endodermis_2"
                                           , "endodermis_1", "epidermis"
                                           , "stele_1.xylem", "stele_2.peri"
                                           , "stele.phloem")
                            , stringsAsFactors = F)

peaks <- merge(peaks, cluster_annot, by = "cluster", all = T)
unique(peaks[,c("cluster", "cellType")])
# cluster      cellType
# 1           0     epidermis
# 2066        1 stele_1.xylem
# 3854        2  stele_2.peri
# 4898        3        cortex
# 5854        4  endodermis_1
# 6807        5          <NA>
#   7025        6          <NA>
#   7120        7 c.e_precursor
# 8308        8  endodermis_3
# 9261        9          <NA>
#   9360       10  endodermis_2
# 10299      11  stele.phloem

peaks[peaks$cluster == 5, "cellType"] <- "NA_1"
peaks[peaks$cluster == 6, "cellType"] <- "NA_2"
peaks[peaks$cluster == 9, "cellType"] <- "NA_3"

peaks$chr <- sub(":.*", "", peaks$peak)
peaks$start <- sub("-.*", "", sub(".*:", "", peaks$peak))
peaks$end <- sub(".*-", "", sub(".*:", "", peaks$peak))

peaks$start <- as.integer(peaks$start)
peaks$end <- as.integer(peaks$end)


# Filter to keep only distal CREs
# Read TAIR10 annotation
genomePath <- "/g/data/zk16/useful/genomes/tair10/"
tair10_txdb <- makeTxDbFromGFF(file = paste0(genomePath, "Arabidopsis_thaliana.TAIR10.59.gff3.gz")
                               , format = "gff3", organism = "Arabidopsis thaliana")
tair10_txdb
# TxDb object:
# Db type: TxDb
# Supporting package: GenomicFeatures
# Data source: /g/data/zk16/useful/genomes/tair10/Arabidopsis_thaliana.TAIR10.59.gff3.gz
# Organism: Arabidopsis thaliana
# Taxonomy ID: 3702
# miRBase build ID: NA
# Genome: NA
# Nb of transcripts: 54013
# Db created by: GenomicFeatures package from Bioconductor
# Creation time: 2024-06-18 11:10:38 +1000 (Tue, 18 Jun 2024)
# GenomicFeatures version at creation time: 1.42.3
# RSQLite version at creation time: 2.2.9
# DBSCHEMAVERSION: 1.2

# +- 500bp as proximal definition
tair10_prom <- promoters(x = tair10_txdb, upstream = 500, downstream = 500)
tair10_exons <- GenomicFeatures::exons(x = tair10_txdb)

peaks.gr <- with(peaks, GRanges(chr, IRanges(start + 1, end)))
# do peaks overlap among themselves
x <- as.data.frame(findOverlaps(peaks.gr, peaks.gr))
x <- x[x$queryHits != x$subjectHits, ]
length(unique(x$queryHits))
# [1] 8323
length(unique(x$subjectHits))
# [1] 8323

# Get cell type-specific peaks
peaks.li <- lapply(unique(peaks$cellType), function(cellType) peaks[peaks$cellType == cellType, ])
names(peaks.li) <- unique(peaks$cellType)

cellT_spec <- list()
for(cellType in names(peaks.li)){
  print(cellType)
  Peaks <- peaks.li[[cellType]]
  Other <- peaks.li[names(peaks.li) != cellType]#
  Other <- do.call("rbind", Other)
  Peaks.gr <- with(Peaks, GRanges(chr, IRanges(start + 1, end)))
  Other.gr <- with(Other, GRanges(chr, IRanges(start + 1, end)))
  x <- as.data.frame(findOverlaps(Peaks.gr, Other.gr))
  Peaks <- Peaks[-unique(x$queryHits),]
  cellT_spec[[cellType]] <- Peaks
}

sapply(cellT_spec, nrow)
# epidermis stele_1.xylem  stele_2.peri        cortex  endodermis_1
# 384           389           155           134           211
# NA_1          NA_2 c.e_precursor  endodermis_3          NA_3
# 98            14           588           296            97
# endodermis_2  stele.phloem
# 381           485
sum(sapply(cellT_spec, nrow))
# [1] 3232

cellT_spec <- do.call("rbind", cellT_spec) 

# match chromosome notation
cellT_spec$chr <- sub("Chr", "", cellT_spec$chr)
cellT_spec.gr <- with(cellT_spec, GRanges(chr, IRanges(start + 1, end)))

# keep distal only
x <- as.data.frame(findOverlaps(cellT_spec.gr, tair10_prom))
distal <- cellT_spec[-unique(x$queryHits),]

distal.gr <- with(distal, GRanges(chr, IRanges(start + 1, end)))
x <- as.data.frame(findOverlaps(distal.gr, tair10_exons))
distalNonEx <- distal[-unique(x$queryHits), ] # not enough CREs

# as.data.frame(table(distalNonEx$cellType))
#             Var1 Freq
# 1  c.e_precursor   83
# 2         cortex   35
# 3   endodermis_1   30
# 4   endodermis_2   61
# 5   endodermis_3   65
# 6      epidermis   41
# 7           NA_1    8
# 8           NA_2    2
# 9           NA_3   16
# 10 stele_1.xylem   62
# 11  stele_2.peri   33
# 12  stele.phloem   46

distal_600bp <- trim_region(distal, 600)

# Save to bed
for(cellType in unique(distal_600bp$cellType)){
  ct_peaks <- distal_600bp[distal_600bp$cellType == cellType, c("chr", "start", "end")]
  out_bed <- paste0(cellType, "_distal_trim600bp.bed")
  write.table(x = unique(ct_peaks), file = out_bed, col.names = F, row.names = F
              , sep = '\t', quote = F)
}
