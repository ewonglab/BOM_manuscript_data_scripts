# human enhancers
#cardi <- read.table(file = "Cardiomyocytes_HEART_distal_nc_500bp.bed", header = F, stringsAsFactors = F, sep='\t')
#ery <- read.table(file = "Erythroblasts_MULTI_distal_nc_500bp.bed", header = F, stringsAsFactors = F, sep='\t')

#cardi$id <- paste0("id", 1:nrow(cardi))
#ery$id <- paste0("id", 1:nrow(ery))

#write.table(x = cardi, file = "Cardiomyocytes_HEART_distal_nc_500bp_0.bed", quote = F, sep="\t", col.names = F, row.names = F)
#write.table(x = ery, file = "Erythroblasts_MULTI_distal_nc_500bp_0.bed", quote = F, sep="\t", col.names = F, row.names = F)

# mouse enhancers
#cardi <- read.table(file = "cardiom_distal_nc_500bp.bed", header = F, stringsAsFactors = F, sep='\t')
#ery <- read.table(file = "erythroid_distal_nc_500bp.bed", header = F, stringsAsFactors = F, sep='\t')

#cardi$V1 <- sub("^", "chr", cardi$V1)
#ery$V1 <- sub("^", "chr", ery$V1)

#cardi <- cardi[,1:3]
#ery <- ery[,1:3]

#cardi$id <- paste0("id", 1:nrow(cardi))
#ery$id <- paste0("id", 1:nrow(ery))

#write.table(x = cardi, file = "cardiom_distal_nc_500bp_0.bed", quote = F, sep="\t", col.names = F, row.names = F)
#write.table(x = ery, file = "erythroid_distal_nc_500bp_0.bed", quote = F, sep="\t", col.names = F, row.names = F)

# map enhancers between hg19 and mm10 using 500bp CREs

# mapping human enhancers
minmatch=0.95
liftOver -minMatch=${minmatch} Cardiomyocytes_HEART_distal_nc_500bp_0.bed \
hg19ToMm10.over.chain.gz \
Cardiomyocytes_HEART_500bp_mm10_${minmatch}.bed Cardiomyocytes_HEART_500bp_unmap
liftOver -minMatch=${minmatch} Erythroblasts_MULTI_distal_nc_500bp_0.bed \
hg19ToMm10.over.chain.gz \
Erythroblasts_MULTI_500bp_mm10_${minmatch}.bed Erythroblasts_MULTI_500bp_unmap

minmatch=0.6
liftOver -minMatch=${minmatch} Cardiomyocytes_HEART_distal_nc_500bp_0.bed \
hg19ToMm10.over.chain.gz \
Cardiomyocytes_HEART_500bp_mm10_${minmatch}.bed Cardiomyocytes_HEART_500bp_unmap
liftOver -minMatch=${minmatch} Erythroblasts_MULTI_distal_nc_500bp_0.bed \
hg19ToMm10.over.chain.gz \
Erythroblasts_MULTI_500bp_mm10_${minmatch}.bed Erythroblasts_MULTI_500bp_unmap

wc -l Cardiomyocytes_HEART_500bp_mm10_*bed
# 405 Cardiomyocytes_HEART_500bp_mm10_0.6.bed
# 49 Cardiomyocytes_HEART_500bp_mm10_0.95.bed
wc -l Erythroblasts_MULTI_500bp_mm10_*bed
# 81 Erythroblasts_MULTI_500bp_mm10_0.6.bed
# 6 Erythroblasts_MULTI_500bp_mm10_0.95.bed


## mapping mouse enhancers
minmatch=0.95
liftOver -minMatch=${minmatch} cardiom_distal_nc_500bp_0.bed \
mm10ToHg19.over.chain.gz \
mouse_cardiom_500bp_hg19_${minmatch}.bed mouse_cardiom_500bp_unmap
liftOver -minMatch=${minmatch} erythroid_distal_nc_500bp_0.bed \
mm10ToHg19.over.chain.gz \
mouse_erythroid_500bp_hg19_${minmatch}.bed mouse_erythroid_500bp_unmap

minmatch=0.6
liftOver -minMatch=${minmatch} cardiom_distal_nc_500bp_0.bed \
mm10ToHg19.over.chain.gz \
mouse_cardiom_500bp_hg19_${minmatch}.bed mouse_cardiom_500bp_unmap
liftOver -minMatch=${minmatch} erythroid_distal_nc_500bp_0.bed \
mm10ToHg19.over.chain.gz \
mouse_erythroid_500bp_hg19_${minmatch}.bed mouse_erythroid_500bp_unmap

wc -l mouse_cardiom_500bp_hg19_*bed
# 624 mouse_cardiom_500bp_hg19_0.6.bed
# 308 mouse_cardiom_500bp_hg19_0.95.bed

wc -l mouse_erythroid_500bp_hg19_*bed
# 578 mouse_erythroid_500bp_hg19_0.6.bed
# 208 mouse_erythroid_500bp_hg19_0.95.bed
