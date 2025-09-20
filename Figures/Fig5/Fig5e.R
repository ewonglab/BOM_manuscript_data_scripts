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
positive_motifs.df <- read.table(file = "E8.25_top30positive_motifs.txt"
                                 , header = T, sep = '\t', stringsAsFactors = F)
positive_motifs.df <- positive_motifs.df[positive_motifs.df$target_celltype != "forebrain",]

### get TF names from the mouse CisBP motif ids
positive_motifs.df$gene_id <- sub("Homo_sapiens.*", "", positive_motifs.df$motif)
positive_motifs.df$gene_id <- sub(".*_2.00_", "", positive_motifs.df$gene_id)
positive_motifs.df$gene_id <- sub("^\\.", "", positive_motifs.df$gene_id)
positive_motifs.df$gene_id <- sub("\\._\\.$", "", positive_motifs.df$gene_id)
positive_motifs.df$gene_id <- sub("._.Arabidopsis_thaliana", "", positive_motifs.df$gene_id)
positive_motifs.df$gene_id <- sub("\\._\\.DBD_.*", "", positive_motifs.df$gene_id)
positive_motifs.df$gene_id <- sub("\\._\\..*", "", positive_motifs.df$gene_id)


#### read shap values, only CARDIOMYOCYTES, ENDOTHELIUM AND NEURAL CREST
shap_f <- c("cardiom_nc_vs_other_mouseCisBP_qval0.5_219T_shap.txt"
            , "endothelium_nc_vs_other_mouseCisBP_qval0.5_532T_shap.txt"
            , "neuralcrest_nc_vs_other_mouseCisBP_qval0.5_23T_shap.txt")
train_f <- c("cardiom_nc_vs_other_mouseCisBP_qval0.5_219T_trainSet.csv"
                 , "endothelium_nc_vs_other_mouseCisBP_qval0.5_532T_trainSet.csv"
                 , "neuralcrest_nc_vs_other_mouseCisBP_qval0.5_23T_trainSet.csv")

shap_tabs <- lapply(shap_f, read.table, header = F, stringsAsFactors = F, sep = '\t')
training <- lapply(train_f, read.csv, header = T, stringsAsFactors = F, row.names = 1)

shap_celltypes <- sub("_.*", "", shap_f)# same order in shap and training files
# "cardiom"     "endothelium" "neuralcrest"
train_celltypes <- sub("_.*", "", train_f)
# "cardiom"     "endothelium" "neuralcrest"

names(shap_tabs) <- shap_celltypes
names(training) <- train_celltypes

lapply(shap_tabs, nrow)
# $cardiom
# [1] 1039
#
# $endothelium
# [1] 2116
#
# $neuralcrest
# [1] 764

# add colnames and row names to SHAP values, also subset to positive set
for(i in 1:length(shap_tabs)){
  rownames(shap_tabs[[i]]) <- rownames(training[[i]])
  colnames(shap_tabs[[i]]) <- colnames(training[[i]])[-ncol(training[[i]])]
  shap_tabs[[i]] <- shap_tabs[[i]][which(training[[i]]$binary_celltype==1),]
}

lapply(shap_tabs, nrow)
# $cardiom
# [1] 516
#
# $endothelium
# [1] 1059
#
# $neuralcrest
# [1] 372

# GET MOTIF IDS FOR THE SELECTED TFs
tf_li1 <- c("Sox12", "Zfp212", "Rara", "Tfap2b", "Tfap2c", "Sox11", "Tfap2a"
            , "Cdc5l", "Junb", "Etv6", "Sox7", "Hsf1", "Irf3", "Fli1", "Erg"
            , "Ets1", "Jund", "Gata4", "Gata6", "Gata5", "Nkx2-5", "Mef2c")

positive_motifs.df$gene_id <- sub("\\.", "-", positive_motifs.df$gene_id)

motifs_1 <- positive_motifs.df[positive_motifs.df$gene_id %in% tf_li1,]#23  6
motifs_1$gene_id[duplicated(motifs_1$gene_id)]#"Gata5"

motifs_1[motifs_1$gene_id == "Gata5", c("motif", "gene_id", "target_celltype")]
# motif gene_id target_celltype
# cardiom.M00167_2.00_Gata5     M00167_2.00_Gata5   Gata5         cardiom
# endothelium.M00167_2.00_Gata5 M00167_2.00_Gata5   Gata5     endothelium

# Gata5 is only important AND highly expressed in cardiomyocytes, remove from endothelium
motifs_1 <- motifs_1[-which(motifs_1$gene_id=="Gata5" &
                              motifs_1$target_celltype =="endothelium"),]#22  6

shap_tfs1 <- list()
for(i in 1:length(shap_tabs)){
  shap_tfs1[[i]] <- shap_tabs[[i]][,motifs_1$motif]
}

for(i in 1:length(shap_tfs1)){
  shap_tfs1[[i]]$celltype <- names(shap_tabs)[i]
}

shap_tfs1 <- do.call("rbind", shap_tfs1)#1947   23
shap_tfs1 <- melt(shap_tfs1, "celltype")

shap_tfs1.mean <- aggregate(value~., shap_tfs1, mean)

shap_tfs1.mean <- tidyr::spread(shap_tfs1.mean, celltype, value)#22  4

shap_tfs1.mean <- merge(shap_tfs1.mean, unique(motifs_1[,c("motif", "target_celltype")])
                        , by.y="motif", by.x="variable", all.x=T)
shap_tfs1.mean <- shap_tfs1.mean[with(shap_tfs1.mean, order(target_celltype)),]

# add rownames (motif id + target cell type)
rownames(shap_tfs1.mean) <- with(shap_tfs1.mean, paste(variable, target_celltype, sep=":"))

shap_tfs1.mean$variable <- NULL
shap_tfs2.mean$variable <- NULL

order_motifs <- function(x){
  tmp <- positive_motifs.df
  tmp$label <- with(tmp, paste(motif, target_celltype, sep=":"))
  tmp <- tmp[tmp$label %in% rownames(x),]
  x <- x[unique(tmp$label),]
  return(x)
}

shap_tfs1.mean <- order_motifs(shap_tfs1.mean)

expr_breaks <- function(x){
  min_shap <- min(unlist(apply(x, 2, min)))
  max_shap <- max(unlist(apply(x, 2, max)))
 
  shap_breaks <- seq(min_shap, max_shap, length.out = 101)
  return(shap_breaks)
}


make_ct_annot <- function(x){
  annot <- data.frame(celltype = sub(".*:", "", rownames(x)))
  annot$celltype <- factor(annot$celltype
                           , levels = c("cardiom", "endothelium", "neuralcrest"))
  rownames(annot) <- rownames(x)
  return(annot)
}

shap_tfs1_row_annot <- make_ct_annot(shap_tfs1.mean)
shap_tfs1.mean$target_celltype <- NULL

# Light Blue: #8dd3c7
#   Sky Blue: #6cb4ce
#   Steel Blue: #4d87a0
#   Navy Blue: #1f4e79

my_palette2 <- colorRampPalette(c("Light Blue", "red"))

pdf("selectedTFs1_SHAP.2.pdf")
pheatmap(shap_tfs1.mean
         , annotation_row = shap_tfs1_row_annot
         , cluster_cols=FALSE
         , cluster_rows=FALSE, show_rownames =T
         , breaks = expr_breaks(shap_tfs1.mean)
         , fontsize = 8
         , color = my_palette2(length(expr_breaks(shap_tfs1.mean))))
dev.off()
