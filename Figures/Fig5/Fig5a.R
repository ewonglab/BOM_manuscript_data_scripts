library(data.table)
library("Signac")
library("Seurat")
library("GenomeInfoDb")
library("GenomicFeatures")
library("AnnotationFilter")
library(GenomicRanges)
library(Biostrings)
library("BSgenome.Hsapiens.UCSC.hg19")
library("Gviz")

setwd("/data/mouse/e8.25/bigwig")

chr <- "chr17"
start_pos <- 26835000
end_pos <- 26871000

bigwig_files <- c("GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_Allantois.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_Cardiomyocytes.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_Endothelium.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_Erythroid.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_ExE_endoderm.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_Forebrain.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_Gut.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_Mesenchyme.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_MidHindbrain.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_Mixed_mesoderm.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_NMP.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_Neural_crest.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_Paraxial_mesoderm.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_Pharyngeal_mesoderm.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_Somitic_mesoderm.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_Spinal_cord.bw"
                  , "GSE133244_embryo_revision1_passQCAll_nuclearGenes_processed_norm_cluster_Surface_ectoderm.bw")
track_name <- sub(".*_cluster_", "", bigwig_files)
track_name <- sub(".bw", "", track_name)

gen <- "mm10"

data_tracks <- lapply(bigwig_files, function(file) {
  print(file)
  track_name <- sub(".*_cluster_", "", file)
  track_name <- sub(".bw", "", track_name)
  # track_name <- sub(".*_", "", file)
  bw_test <- import.bw(file, as="GRanges")
  accDT <- DataTrack(bw_test, chromosome=chr, name=track_name)
  return(accDT)
})

gtrack <- GenomeAxisTrack()


mm10_txdb <- makeTxDbFromGFF(file="/g/data/zk16/cc3704/mouse_data/Mus_musculus.GRCm38.92.gtf"
                             , format="gtf",organism="Mus musculus")
region_gr <- GRanges(17, IRanges(start_pos+1, end_pos))
mouse_genes <- genes(mm10_txdb)
x <- as.data.frame(findOverlaps(mouse_genes, region_gr))
mouse_genes <- mouse_genes[x$queryHits]


options(ucscChromosomeNames=FALSE)
genetrack <- GeneRegionTrack("/g/data/zk16/cc3704/mouse_data/Mus_musculus.GRCm38.92.gtf", chromosome="17")
genetrack <- genetrack[gene(genetrack) %in%  mouse_genes@elementMetadata$gene_id ]
seqlevels(genetrack@range) <- sub("17", "chr17", seqlevels(genetrack@range))
chromosome(genetrack) <- "chr17"
genome(genetrack) <- "mm10"
peakTrack <- AnnotationTrack(start = c(26850448, 26865164, 26868664)
                             , end = c(26850950, 26865668, 26869213)
                             , chromosome = c("chr17", "chr17",  "chr17")
                             , genome = "mm10", name = "peak")
unique(feature(genetrack))
# [1] "gene"            "transcript"      "exon"            "CDS"            
# [5] "start_codon"     "stop_codon"      "five_prime_utr"  "three_prime_utr"

unique(gene(genetrack))
# [1] "ENSMUSG00000015579" "ENSMUSG00000078838" "ENSMUSG00000092539"

pdf("Nkx2.5_tracks.5.pdf")
plotTracks(c(gtrack, data_tracks, genetrack[feature(genetrack) %in%
                                              c("exon", "five_prime_utr", "three_prime_utr", "CDS") ]
             , peakTrack)
           , from=start_pos+1, to=end_pos, type = "histogram", ylim=c(0,40)
           ,   background.panel = "white", cex.axis = 0.6)
dev.off()
