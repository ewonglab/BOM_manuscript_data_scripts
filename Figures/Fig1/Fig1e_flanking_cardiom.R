suppressMessages({
    library(stringr)
    library(GenomicRanges)
    library("GenomeInfoDb")
    library("GenomicFeatures")
    library(ggplot2)
    library(dplyr)
    library(data.table)
})

setwd("/data/mouse/e8.25")
mouse <- read.csv(file = "Pijuan_etal_table_S6.csv", header =T, stringsAsFactors=F)
mouse <- mouse[,c("peak_chr", "peak_start", "peak_end", "peakID", "celltype_specificity")]
mouse <- unique(mouse)

mouse_celltype <- mouse[!is.na(mouse$celltype_specificity),]#19453     5
mouse_celltype.spc <- mouse_celltype[grep(";", mouse_celltype$celltype_specificity, invert = TRUE),]# 15504     5

mm10_txdb <- makeTxDbFromGFF(file="/g/data/zk16/cc3704/mouse_data/Mus_musculus.GRCm38.92.gtf"
                             , format="gtf",organism="Mus musculus")
mm10_prox <- promoters(x = mm10_txdb, upstream = 1000, downstream = 1000)

mouse_celltype.spc_gr <- with(mouse_celltype.spc, GRanges(peak_chr, IRanges(peak_start+1, peak_end)))

# DISTAL CRES
x <- as.data.frame(findOverlaps(mouse_celltype.spc_gr, mm10_prox))
mouse_distal <- mouse_celltype.spc[-unique(x$queryHits),]#13014     5
dim(mouse_distal)
as.data.frame(table(mouse_distal$celltype_specificity))

# check overlap between peaks
peak_overlaps <- function(celltype){
  target <- mouse_distal[mouse_distal$celltype_specificity == celltype,]
  other <- mouse_distal[mouse_distal$celltype_specificity != celltype,]
  target_gr <- with(target, GRanges(peak_chr, IRanges(peak_start+1, peak_end)))
  other_gr <- with(other, GRanges(peak_chr, IRanges(peak_start+1, peak_end)))
  x <- as.data.frame(findOverlaps(target_gr, other_gr))
  print(paste(celltype, nrow(x), "overlapping pairs"))
}

celltypes_li <- c("Pharyngeal mesoderm", "Forebrain", "Gut"
                  , "Paraxial mesoderm", "Spinal cord", "Somitic mesoderm"
                  , "Endothelium", "Cardiomyocytes", "Mid/Hindbrain"
                  , "Neural crest", "Erythroid", "NMP", "Mesenchyme"
                  , "Surface ectoderm", "Mixed mesoderm", "Allantois"
                  , "ExE endoderm", "Notochord")

sapply(celltypes_li, peak_overlaps)#0 
mouse_distal_gr <- with(mouse_distal, GRanges(peak_chr, IRanges(peak_start+1, peak_end)))
mm10_exons <- GenomicFeatures::exons(x = mm10_txdb)
x <- as.data.frame(findOverlaps(mouse_distal_gr, mm10_exons))
mouse_distal_nc <- mouse_distal[-unique(x$queryHits),]#12161     5

mouse_distal_nc$flanking_start <- with(mouse_distal_nc, peak_start - 2000)
mouse_distal_nc$flanking_end <- with(mouse_distal_nc, peak_end + 2000)

# upstream regions
flanking_ups <- mouse_distal_nc[, c("peak_chr", "flanking_start", "peak_start", "peakID", "celltype_specificity")]
# downstream regions
flanking_dwn <- mouse_distal_nc[, c("peak_chr", "peak_end", "flanking_end", "peakID", "celltype_specificity")]

colnames(flanking_ups) <- c("peak_chr", "ups_flanking_start", "ups_flanking_end", "CRE_ID", "celltype_specificity")
colnames(flanking_dwn) <- c("peak_chr", "dwn_flanking_start", "dwn_flanking_end", "CRE_ID", "celltype_specificity")

summary(with(flanking_ups, ups_flanking_end - ups_flanking_start))
head(flanking_ups, 2)

chr_sizes <- read.table(file = "/g/data/zk16/cc3704/mouse_data/mouse_sizes_primary_genome_tab.txt"
                        , header = FALSE, stringsAsFactors = FALSE, sep = '\t')
# table(flanking_dwn$dwn_flanking_start < 0) 
tmp <- merge(flanking_dwn, chr_sizes, by.x = "peak_chr", by.y = "V1")

flanking_ups_gr <- with(flanking_ups, GRanges(peak_chr, IRanges(ups_flanking_start + 1, ups_flanking_end)))
flanking_dwn_gr <- with(flanking_dwn, GRanges(peak_chr, IRanges(dwn_flanking_start + 1, dwn_flanking_end)))

# same width as trimmed CREs; 50bp stride
flanking_ups_windows <- slidingWindows(flanking_ups_gr, width = 500, step = 50)
flanking_dwn_windows <- slidingWindows(flanking_dwn_gr, width = 500, step = 50)

flanking_ups_windows <- unlist(flanking_ups_windows)
flanking_dwn_windows <- unlist(flanking_dwn_windows)
flanking_ups_windows.df <- as.data.frame(flanking_ups_windows)
flanking_dwn_windows.df <- as.data.frame(flanking_dwn_windows)

flanking_ups_windows.df$flankID <- with(flanking_ups_windows.df, paste(seqnames, paste(start, end, sep = "_"), sep = ":"))
flanking_dwn_windows.df$flankID <- with(flanking_dwn_windows.df, paste(seqnames, paste(start, end, sep = "_"), sep = ":"))

ups_ov <- as.data.frame(findOverlaps(flanking_ups_gr, flanking_ups_windows))
dwn_ov <- as.data.frame(findOverlaps(flanking_dwn_gr, flanking_dwn_windows))

ups_to_windows <- cbind(flanking_ups[ups_ov$queryHits, ], flanking_ups_windows.df[ups_ov$subjectHits, ])
dnw_to_windows <- cbind(flanking_dwn[dwn_ov$queryHits, ], flanking_dwn_windows.df[dwn_ov$subjectHits, ])

setwd("/scratch/zk16/cc3704/E8.25_flanking")
getwd()
write.table(x = ups_to_windows, file = "dist_nc_upstream_windows_500bp_50str.txt", quote = FALSE, sep = "\t")
write.table(x = dnw_to_windows, file = "dist_nc_downstream_windows_500bp_50str.txt", quote = FALSE, sep = "\t")

mouse <- read.csv(file = "Pijuan_etal_table_S6.csv", header =T, stringsAsFactors=F)
mouse <- mouse[,c("peak_chr", "peak_start", "peak_end", "peakID", "celltype_specificity")]
mouse <- unique(mouse)
mouse_celltype <- mouse[!is.na(mouse$celltype_specificity),]#19453     5
mouse_celltype.spc <- mouse_celltype[grep(";", mouse_celltype$celltype_specificity, invert = TRUE),]# 15504     5

# gr object of enhancers
mouse_celltype.spc_gr <- with(mouse_celltype.spc, GRanges(peak_chr, IRanges(peak_start+1, peak_end)))
x <- as.data.frame(findOverlaps(mouse_celltype.spc_gr, mm10_prox))
mouse_distal <- mouse_celltype.spc[-unique(x$queryHits),]#13014     5
mouse_distal_gr <- with(mouse_distal, GRanges(peak_chr, IRanges(peak_start+1, peak_end)))
mm10_exons <- GenomicFeatures::exons(x = mm10_txdb)
x <- as.data.frame(findOverlaps(mouse_distal_gr, mm10_exons))#
mouse_distal_nc <- mouse_distal[-unique(x$queryHits),]#12161     5
dim(mouse_distal_nc)

mouse_distal_nc_gr <- with(mouse_distal_nc, GRanges(peak_chr, IRanges(peak_start + 1, peak_end)))
ups_to_windows_gr <- with(ups_to_windows, GRanges(peak_chr, IRanges(ups_flanking_start + 1, ups_flanking_end)))
dnw_to_windows_gr <- with(dnw_to_windows, GRanges(peak_chr, IRanges(dwn_flanking_start + 1, dwn_flanking_end)))

x <- as.data.frame(findOverlaps(mouse_distal_nc_gr, ups_to_windows_gr))
y <- as.data.frame(findOverlaps(mouse_distal_nc_gr, dnw_to_windows_gr))

CREs_ov_ups <- cbind(mouse_distal_nc[x$queryHits, ], ups_to_windows[x$subjectHits, ])
CREs_ov_dwn <- cbind(mouse_distal_nc[y$queryHits, ], dnw_to_windows[y$subjectHits, ])

cardiom_CREs <- mouse_distal_nc[mouse_distal_nc$celltype_specificity == "Cardiomyocytes", ]
cardiom_windows_ups <- ups_to_windows[ups_to_windows$celltype_specificity == "Cardiomyocytes", ]
cardiom_windows_dwn <- dnw_to_windows[dnw_to_windows$celltype_specificity == "Cardiomyocytes", ]

cardiom_CREs_gr <- with(cardiom_CREs, GRanges(peak_chr, IRanges(peak_start + 1, peak_end)))
cardiom_ups_gr <- with(cardiom_windows_ups, GRanges(peak_chr, IRanges(ups_flanking_start + 1, ups_flanking_end)))
cardiom_dnw_gr <- with(cardiom_windows_dwn, GRanges(peak_chr, IRanges(dwn_flanking_start + 1, dwn_flanking_end)))

x <- as.data.frame(findOverlaps(cardiom_CREs_gr, cardiom_ups_gr))
y <- as.data.frame(findOverlaps(cardiom_CREs_gr, cardiom_dnw_gr))

cardiomCREs_ov_ups <- (cbind(cardiom_CREs[x$queryHits, ], cardiom_windows_ups[x$subjectHits, ]))
cardiomCREs_ov_dwn <- (cbind(cardiom_CREs[y$queryHits, ], cardiom_windows_dwn[y$subjectHits, ]))
dim(cardiomCREs_ov_ups)
dim(cardiomCREs_ov_dwn)

cardiom_windows_ups.unique <- unique(cardiom_windows_ups[, c('seqnames', 'start', 'end', 'flankID')])
cardiom_windows_dwn.unique <- unique(cardiom_windows_dwn[, c('seqnames', 'start', 'end', 'flankID')])

cardiom_windows_ups.unique$start <- cardiom_windows_ups.unique$start - 1
cardiom_windows_dwn.unique$start <- cardiom_windows_dwn.unique$start - 1

cardiom_CREs_gr <- with(cardiom_CREs, GRanges(peak_chr, IRanges(peak_start + 1, peak_end)))
cardiom_windows_ups.unique_gr <- with(cardiom_windows_ups.unique, GRanges(seqnames, IRanges(start + 1, end)))
cardiom_windows_dwn.unique_gr <- with(cardiom_windows_dwn.unique, GRanges(seqnames, IRanges(start + 1, end)))

x <- as.data.frame(findOverlaps(cardiom_CREs_gr, cardiom_windows_ups.unique_gr))

cardiom_UPS <- cardiom_windows_ups.unique[-unique(x$subjectHits), ]
cardiom_DWN <- cardiom_windows_dwn.unique[-unique(x$subjectHits), ]

# update window IDs with 0-based coordinates
cardiom_UPS$flankID <- with(cardiom_UPS, paste(seqnames, paste(start, end, sep = "_"), sep = ":"))
cardiom_DWN$flankID <- with(cardiom_DWN, paste(seqnames, paste(start, end, sep = "_"), sep = ":"))

write.table(x = cardiom_UPS
            , file = "cardiom_windows_ups_50str.bed", col.names = FALSE, row.names = FALSE
            , sep = '\t', quote = FALSE)
write.table(x = cardiom_DWN
            , file = "cardiom_windows_dwn_50str.bed", col.names = FALSE, row.names = FALSE
            , sep = '\t', quote = FALSE)
