library("GenomeInfoDb")
library("AnnotationFilter")
library("SummarizedExperiment")
library(ggplot2)
library("Seurat")
library(Matrix)
library("tidyr")
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(graphics)

setwd("/xgb/e8.25")
# read files containing mean SHAP for the positive class of each model
positive_motifs.df <- read.table(file = "E8.25_top30positive_motifs.txt"
                                 , header = T, sep = '\t', stringsAsFactors = F)

rna_E8.25 <- readRDS("/data/mouse/EmbryoTimecourse2018/rna_E8.25_seurat_filt.rds")
RNAdata <- GetAssayData(object = rna_E8.25, slot = "data")
RNAdata <- as.data.frame(as.matrix(RNAdata))
RNAdata$gene_id <- rownames(RNAdata)
#remove ensemble id
RNAdata$gene_id <- sub(".*?-", "", RNAdata$gene_id)

### get TF names from the mouse CisBP motif ids
positive_motifs.df$gene_id <- sub("Homo_sapiens.*", "", positive_motifs.df$motif)
positive_motifs.df$gene_id <- sub(".*_2.00_", "", positive_motifs.df$gene_id)
positive_motifs.df$gene_id <- sub("^\\.", "", positive_motifs.df$gene_id)
positive_motifs.df$gene_id <- sub("\\._\\.$", "", positive_motifs.df$gene_id)
positive_motifs.df$gene_id <- sub("._.Arabidopsis_thaliana", "", positive_motifs.df$gene_id)
positive_motifs.df$gene_id <- sub("\\._\\.DBD_.*", "", positive_motifs.df$gene_id)
positive_motifs.df$gene_id <- sub("\\._\\..*", "", positive_motifs.df$gene_id)

# SELECTIVE TFs BASED ON MOTIF IMPORTANCE AND HIGH EXPRESSION IN TARGET CELL TYPES
# ( INLCUDING MOTIFS IMPORTANT IN CARDIOMYOCYTES, ENDOTHELIUM AND NEURAL CREST)

tf_li1 <- c("Sox12", "Zfp212", "Rara", "Tfap2b", "Tfap2c", "Sox11", "Tfap2a"
            , "Cdc5l", "Junb", "Etv6", "Sox7", "Hsf1", "Irf3", "Fli1", "Erg"
            , "Ets1", "Jund", "Gata4", "Gata6", "Gata5", "Nkx2-5", "Mef2c")

tf_li2 <- c("Bcl6b", "Sox12", "Zfp212", "Sox4", "Rara", "Tfap2b", "Tfap2c", "Sox11"
            , "Hbp1", "Tfap2a", "Cdc5l", "Junb", "Etv6", "Sox7", "Hsf1", "Irf3"
            , "Fli1", "Erg", "Ets1", "Jund", "Gata4", "Meis3", "Gata6", "Gata5"
            , "Nkx2-5", "Mef2c")

table(tf_li1 %in% positive_motifs.df$gene_id)
# FALSE  TRUE
# 1    21
table(tf_li2 %in% positive_motifs.df$gene_id)
# FALSE  TRUE
# 1    25

tf_li1[!tf_li1 %in% positive_motifs.df$gene_id]
#"Nkx2-5"
tf_li2[!tf_li2 %in% positive_motifs.df$gene_id]
#"Nkx2-5"

unique(grep("\\.", positive_motifs.df$gene_id, value=T))
#"Nkx2.5" "Nkx2.1" "Nkx2.9" "Nkx2.3" "Nkx2.2" "Nkx2.4"

# are top motif TFs present in gene expression data?
table(positive_motifs.df$gene_id %in% RNAdata$gene_id)

positive_motifs.df$gene_id[!positive_motifs.df$gene_id %in% RNAdata$gene_id]
# "Nkx2.5" "Nkx2.1" "Nkx2.9" "Nkx2.3" "Nkx2.2" "Nkx2.4"

positive_motifs.df$gene_id <- sub("\\.", "-", positive_motifs.df$gene_id)

### subset gene expression: keep only genes corresponting to top TFs
## subset expression for selected TFs
RNAdata <- RNAdata[RNAdata$gene_id %in% positive_motifs.df$gene_id,]#109 15936

RNAdata_tfs1 <- RNAdata[RNAdata$gene_id %in% tf_li1,]# 22 15936
RNAdata_tfs2 <- RNAdata[RNAdata$gene_id %in% tf_li2,]# 26 15936

positive_motifs.df <- positive_motifs.df[positive_motifs.df$target_celltype !="forebrain",]

RNAdata_tfs1 <- merge(RNAdata_tfs1, unique(positive_motifs.df[,c("gene_id", "target_celltype")])
                      , by="gene_id", all.x=T)# 23 15937
RNAdata_tfs2 <- merge(RNAdata_tfs2, unique(positive_motifs.df[,c("gene_id", "target_celltype")])
                      , by="gene_id", all.x=T)# 27 15937

RNAdata_tfs1 <- RNAdata_tfs1[-which(RNAdata_tfs1$gene_id=="Gata5" &
                                 RNAdata_tfs1$target_celltype =="endothelium"),]#22 15937

RNAdata_tfs2 <- RNAdata_tfs2[-which(RNAdata_tfs2$gene_id=="Gata5" &
                                 RNAdata_tfs2$target_celltype =="endothelium"),]#26 15937


## melt and add cell identities
RNAdata_tfs1_long <- melt(RNAdata_tfs1, c("gene_id", "target_celltype"))#350570      4
RNAdata_tfs2_long <- melt(RNAdata_tfs2, c("gene_id", "target_celltype"))#414310      4

## read metadata
# add cell type identity to every cell
meta <- read.table(file = "/data/mouse/EmbryoTimecourse2018/atlas/meta.tab"
                   , header =T, stringsAsFactors = F, sep ='\t')#139331     28

RNAdata_tfs1_long$variable <- sub("-", "_", RNAdata_tfs1_long$variable)
RNAdata_tfs2_long$variable <- sub("-", "_", RNAdata_tfs2_long$variable)

RNAdata_tfs1_long <- merge(RNAdata_tfs1_long, unique(meta[,c("cell", "celltype")])
                           , by.x="variable", by.y="cell")#350570      5
RNAdata_tfs2_long <- merge(RNAdata_tfs2_long, unique(meta[,c("cell", "celltype")])
                           , by.x="variable", by.y="cell")#414310      5


## CALCULATE MEAN PER TF, PER TARGET CELL TYPE, PER CELL TYPE
tfs1_means <- aggregate(value~., RNAdata_tfs1_long[,c("gene_id", "target_celltype"
                                                    , "value", "celltype")], mean)#726   4
tfs2_means <- aggregate(value~., RNAdata_tfs2_long[,c("gene_id", "target_celltype"
                                                    , "value", "celltype")], mean)#858   4

# make matrix
tfs1_means <- tidyr::spread(tfs1_means, celltype, value)#22 35
tfs2_means <- tidyr::spread(tfs2_means, celltype, value)#26 35

# order by cell type
tfs1_means <- tfs1_means[with(tfs1_means, order(target_celltype)),]
tfs2_means <- tfs2_means[with(tfs2_means, order(target_celltype)),]

## remove cell types not present in BOM models

remove_ct <- c("Anterior Primitive Streak", "Blood progenitors 1"
               , "Blood progenitors 2", "Caudal Mesoderm", "Caudal epiblast"
               , "Caudal neurectoderm", "Def. endoderm", "ExE mesoderm"
               , "Forebrain/Midbrain/Hindbrain", "Haematoendothelial progenitors"
               , "Intermediate mesoderm", "Notochord", "PGC", "Parietal endoderm"
               , "Rostral neurectoderm", "Visceral endoderm", "ExE ectoderm")

tfs1_means <- tfs1_means[,colnames(tfs1_means)[!colnames(tfs1_means) %in% remove_ct]]#22 18
tfs2_means <- tfs2_means[,colnames(tfs2_means)[!colnames(tfs2_means) %in% remove_ct]]#26 18

rownames(tfs1_means) <- with(tfs1_means, paste(gene_id, target_celltype, sep=":"))
rownames(tfs2_means) <- with(tfs2_means, paste(gene_id, target_celltype, sep=":"))

tfs1_means$gene_id <- NULL
tfs1_means$target_celltype <- NULL

tfs2_means$gene_id <- NULL
tfs2_means$target_celltype <- NULL

expr_breaks <- function(x){
  min_expr <- min(unlist(apply(x, 2, min)))
  max_expr <- max(unlist(apply(x, 2, max)))
 
  expr_breaks <- seq(min_expr, max_expr, length.out = 101)
  return(expr_breaks)
}

order_motifs <- function(x){
  tmp <- positive_motifs.df
  tmp$label <- with(tmp, paste(gene_id, target_celltype, sep=":"))
  tmp <- tmp[tmp$label %in% rownames(x),]
  x <- x[unique(tmp$label),]
  return(x)
}

tfs1_means <- order_motifs(tfs1_means)
tfs2_means <- order_motifs(tfs2_means)

make_ct_annot <- function(x){
  annot <- data.frame(celltype = sub(".*:", "", rownames(x)))
  annot$celltype <- factor(annot$celltype
                           , levels = c("cardiom", "endothelium", "neuralcrest"))
  rownames(annot) <- rownames(x)
  return(annot)
}

tfs1_row_annot <- make_ct_annot(tfs1_means)
tfs2_row_annot <- make_ct_annot(tfs2_means)

tfs1_breaks <- expr_breaks(tfs1_means[,c("Cardiomyocytes", "Endothelium", "Neural crest")])
my_palette2 <- colorRampPalette(c("Light Blue", "red"))

pdf("selectedTFs1_expr.2.pdf")
pheatmap((tfs1_means[,c("Cardiomyocytes", "Endothelium", "Neural crest")])
         , annotation_row = tfs1_row_annot
         , cluster_cols=FALSE
         , cluster_rows=FALSE, show_rownames =T
         , breaks = tfs1_breaks, color=my_palette2(length(tfs1_breaks))
         , fontsize = 8)
dev.off()
