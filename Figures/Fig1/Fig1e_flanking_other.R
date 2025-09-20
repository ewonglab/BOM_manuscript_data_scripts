suppressMessages({
    library(stringr)
    library(GenomicRanges)
    library("GenomeInfoDb")
    library("GenomicFeatures")
    library(ggplot2)
    library(dplyr)
    library(data.table)
})

ups_to_windows <- read.table(file = "dist_nc_upstream_windows_500bp_50str.txt", header = TRUE
                             , stringsAsFactors = FALSE, sep = "\t")
dnw_to_windows <- read.table(file = "dist_nc_downstream_windows_500bp_50str.txt", header = TRUE
                             , stringsAsFactors = FALSE, sep = "\t")

mouse <- read.csv(file = "Pijuan_etal_table_S6.csv"
                  , header = TRUE, stringsAsFactors = FALSE)
mouse <- mouse[,c("peak_chr", "peak_start", "peak_end", "peakID", "celltype_specificity")]
mouse <- unique(mouse)
mouse_celltype <- mouse[!is.na(mouse$celltype_specificity),]
mouse_celltype.spc <- mouse_celltype[grep(";", mouse_celltype$celltype_specificity, invert = TRUE),]

mm10_txdb <- makeTxDbFromGFF(file="Mus_musculus.GRCm38.92.gtf"
                             , format="gtf",organism="Mus musculus")
mm10_prox <- promoters(x = mm10_txdb, upstream = 1000, downstream = 1000)

mouse_celltype.spc_gr <- with(mouse_celltype.spc, GRanges(peak_chr, IRanges(peak_start+1, peak_end)))
x <- as.data.frame(findOverlaps(mouse_celltype.spc_gr, mm10_prox))
mouse_distal <- mouse_celltype.spc[-unique(x$queryHits),]#13014     5
mouse_distal_gr <- with(mouse_distal, GRanges(peak_chr, IRanges(peak_start+1, peak_end)))
mm10_exons <- GenomicFeatures::exons(x = mm10_txdb)
x <- as.data.frame(findOverlaps(mouse_distal_gr, mm10_exons))#
mouse_distal_nc <- mouse_distal[-unique(x$queryHits),]

mouse_distal_nc_gr <- with(mouse_distal_nc, GRanges(peak_chr, IRanges(peak_start + 1, peak_end)))
ups_to_windows_gr <- with(ups_to_windows, GRanges(peak_chr, IRanges(ups_flanking_start + 1, ups_flanking_end)))
dnw_to_windows_gr <- with(dnw_to_windows, GRanges(peak_chr, IRanges(dwn_flanking_start + 1, dwn_flanking_end)))

x <- as.data.frame(findOverlaps(mouse_distal_nc_gr, ups_to_windows_gr))
y <- as.data.frame(findOverlaps(mouse_distal_nc_gr, dnw_to_windows_gr))

CREs_ov_ups <- cbind(mouse_distal_nc[x$queryHits, ], ups_to_windows[x$subjectHits, ])
CREs_ov_dwn <- cbind(mouse_distal_nc[y$queryHits, ], dnw_to_windows[y$subjectHits, ])


get_flanking <- function(celltype){
    # subset to cardiomyocyte CREs and flanking windows
    celltype_CREs <- mouse_distal_nc[mouse_distal_nc$celltype_specificity == celltype, ]
    celltype_windows_ups <- ups_to_windows[ups_to_windows$celltype_specificity == celltype, ]
    celltype_windows_dwn <- dnw_to_windows[dnw_to_windows$celltype_specificity == celltype, ]
    
    print("dim CREs")
    print(dim(celltype_CREs))
    print("dim upstream windows")
    print(dim(celltype_windows_ups))
    print("dim  downstream windows")
    print(dim(celltype_windows_dwn))
    
    # gr objects of CREs and flanking regions
    celltype_CREs_gr <- with(celltype_CREs, GRanges(peak_chr, IRanges(peak_start + 1, peak_end)))
    celltype_ups_gr <- with(celltype_windows_ups, GRanges(peak_chr, IRanges(ups_flanking_start + 1, ups_flanking_end)))
    celltype_dnw_gr <- with(celltype_windows_dwn, GRanges(peak_chr, IRanges(dwn_flanking_start + 1, dwn_flanking_end)))
    
    x <- as.data.frame(findOverlaps(celltype_CREs_gr, celltype_ups_gr))
    y <- as.data.frame(findOverlaps(celltype_CREs_gr, celltype_dnw_gr))
    
    # dim(x)
    # dim(y)

    celltypeCREs_ov_ups <- (cbind(celltype_CREs[x$queryHits, ], celltype_windows_ups[x$subjectHits, ]))
    celltypeCREs_ov_dwn <- (cbind(celltype_CREs[y$queryHits, ], celltype_windows_dwn[y$subjectHits, ]))
    
    # print("dim celltypeCREs_ov_ups")
    # print(dim(celltypeCREs_ov_ups))
    # print("dim celltypeCREs_ov_dwn")
    # print(dim(celltypeCREs_ov_dwn))
    # print(dim(x))
    # print(dim(y))
    
    # correct numbers
    print("Number of CREs overlapping UPS flanking")
    print(length(unique(celltypeCREs_ov_ups$peakID)))
    print("Number of UPS flanking overlapping CREs")
    print(nrow(unique(celltypeCREs_ov_ups[, c("peak_chr", "ups_flanking_start", "ups_flanking_end")])))

    print("Number of CREs overlapping DWN flanking")
    print(length(unique(celltypeCREs_ov_dwn$peakID)))
    print("Number of DWN flanking overlapping CREs")
    print(nrow(unique(celltypeCREs_ov_dwn[, c("peak_chr", "dwn_flanking_start", "dwn_flanking_end")])))
    
    # there are duplicated windows coordinates linked to different CREs; CREs do not overlap but flanking regions do
    celltype_windows_ups.unique <- unique(celltype_windows_ups[, c('seqnames', 'start', 'end', 'flankID')])
    celltype_windows_dwn.unique <- unique(celltype_windows_dwn[, c('seqnames', 'start', 'end', 'flankID')])
    
    # windows to 0-based
    celltype_windows_ups.unique$start <- celltype_windows_ups.unique$start - 1
    celltype_windows_dwn.unique$start <- celltype_windows_dwn.unique$start - 1

    print("windows width")
    print(unique(with(celltype_windows_ups.unique, end - start)))
    print(unique(with(celltype_windows_dwn.unique, end - start)))

    # remove windows overlapping CREs
    celltype_CREs_gr <- with(celltype_CREs, GRanges(peak_chr, IRanges(peak_start + 1, peak_end)))
    celltype_windows_ups.unique_gr <- with(celltype_windows_ups.unique, GRanges(seqnames, IRanges(start + 1, end)))
    celltype_windows_dwn.unique_gr <- with(celltype_windows_dwn.unique, GRanges(seqnames, IRanges(start + 1, end)))
    
    x <- as.data.frame(findOverlaps(celltype_CREs_gr, celltype_windows_ups.unique_gr))
    celltype_UPS <- celltype_windows_ups.unique[-unique(x$subjectHits), ]

    x <- as.data.frame(findOverlaps(celltype_CREs_gr, celltype_windows_dwn.unique_gr))
    celltype_DWN <- celltype_windows_dwn.unique[-unique(x$subjectHits), ]

    # update window IDs with 0-based coordinates
    celltype_UPS$flankID <- with(celltype_UPS, paste(seqnames, paste(start, end, sep = "_"), sep = ":"))
    celltype_DWN$flankID <- with(celltype_DWN, paste(seqnames, paste(start, end, sep = "_"), sep = ":"))
    
    print("dim celltype_UPS")
    print(dim(celltype_UPS))
    print("dim unique celltype_UPS")
    print(dim(unique(celltype_UPS)))
    
    print("dim celltype_DWN")
    print(dim(celltype_DWN))
    print("dim celltype_DWN")
    print(dim(unique(celltype_DWN)))
    
    return(list(upstream = celltype_UPS, downstream = celltype_DWN))

}

Endothelium <- get_flanking("Endothelium")
Forebrain <- get_flanking("Forebrain")
Gut <- get_flanking("Gut")
Erythroid <- get_flanking("Erythroid")
ExEendoderm <- get_flanking("ExE endoderm")

write.table(x = Endothelium$upstream
            , file = "Endothelium_windows_ups_50str.bed", col.names = FALSE, row.names = FALSE
            , sep = '\t', quote = FALSE)
write.table(x = Endothelium$downstream
            , file = "Endothelium_windows_dwn_50str.bed", col.names = FALSE, row.names = FALSE
            , sep = '\t', quote = FALSE)

write.table(x = Forebrain$upstream
            , file = "Forebrain_windows_ups_50str.bed", col.names = FALSE, row.names = FALSE
            , sep = '\t', quote = FALSE)
write.table(x = Forebrain$downstream
            , file = "Forebrain_windows_dwn_50str.bed", col.names = FALSE, row.names = FALSE
            , sep = '\t', quote = FALSE)

write.table(x = Gut$upstream
            , file = "Gut_windows_ups_50str.bed", col.names = FALSE, row.names = FALSE
            , sep = '\t', quote = FALSE)
write.table(x = Gut$downstream
            , file = "Gut_windows_dwn_50str.bed", col.names = FALSE, row.names = FALSE
            , sep = '\t', quote = FALSE)

write.table(x = Erythroid$upstream
            , file = "Erythroid_windows_ups_50str.bed", col.names = FALSE, row.names = FALSE
            , sep = '\t', quote = FALSE)
write.table(x = Erythroid$downstream
            , file = "Erythroid_windows_dwn_50str.bed", col.names = FALSE, row.names = FALSE
            , sep = '\t', quote = FALSE)
write.table(x = ExEendoderm$upstream
            , file = "ExEendoderm_windows_ups_50str.bed", col.names = FALSE, row.names = FALSE
            , sep = '\t', quote = FALSE)
write.table(x = ExEendoderm$downstream
            , file = "ExEendoderm_windows_dwn_50str.bed", col.names = FALSE, row.names = FALSE
            , sep = '\t', quote = FALSE)
