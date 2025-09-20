library(pheatmap)
library(data.table)
library(RColorBrewer)
library(viridis)

motif_mean <- function(x){
  x$celltype <- NULL
  colnames(x) <- sub("_NA$", "", colnames(x))
  counts_1 <- x[x$binary_celltype==1,]
  counts_0 <- x[x$binary_celltype==0,]
  counts_1$binary_celltype <- NULL
  counts_0$binary_celltype <- NULL
  #column sums
  counts_1_colSums <- as.data.frame(apply(counts_1, 2, mean))
  counts_0_colSums <- as.data.frame(apply(counts_0, 2, mean))
  counts_colSums <- merge(counts_1_colSums, counts_0_colSums, by=0)
  colnames(counts_colSums) <- c("motif", "sum_1", "sum_0")
  counts_colSums$direction <- ifelse(counts_colSums$sum_1 > counts_colSums$sum_0
                                     , "enriched_in_class_1", "enriched_in_class_0")
  rownames(counts_colSums) <- counts_colSums$motif
  return(counts_colSums)  
}


add_missing_motifs <- function(motifs, counts){
  missing <- setdiff(motifs, colnames(counts))
  if(length(missing) > 0){
    for(i in missing){
      counts[,i] <- 0
    }
  }
  return(counts)
}

std1 <- function(x){
  return ((x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

# function to normalize motif counts
normalize_counts <- function(counts){
  motifs.std <- counts
  motifs.std$peakid <- rownames(motifs.std)
  #reshape counts
  motifs.std <- reshape2::melt(motifs.std)#432620      3
  #Using peakid as id variables
 
  colnames(motifs.std) <-c("peakid","motif", "count")
  motifs.std <- as.data.table(motifs.std)
  motifs.std[, count_std := std1(count), by = "motif"]
 
  motifs.std$count <- NULL
  motifs.std <- tidyr::spread(motifs.std, motif, count_std)
  motifs.std <- as.data.frame(motifs.std)#3528   38
  rownames(motifs.std) <- motifs.std$peakid
  motifs.std$peakid <- NULL
  motifs.std <- motifs.std[,colnames(counts)]
  motifs.std <- motifs.std[rownames(counts),]
  # when a motif has zero counts across all CREs the result is NaN
  # changing NaN to zero
  for(i in 1:ncol(motifs.std)){
    if(unique(is.na(motifs.std[,i]))){
      print(paste("Motif", colnames(motifs.std)[i], "has zero counts. Transforming NaNs to zeros"))
      motifs.std[,i] <- 0
    }
  }
  return(motifs.std)
}

combine_counts <- function(model_counts, test_counts, shuffled_counts
                           , top_motifs, pred, test_align, shuffled_pred
                           ,shuffled_align, model_sp, test_sp, out_name){ #,
  top_motifs <- sub("_NA$", "", top_motifs)
  colnames(model_counts) <- sub("_NA$", "", colnames(model_counts))
  colnames(test_counts) <- sub("_NA$", "", colnames(test_counts))
 
  # only possitive instances in model and test counts and test instances
  model_counts <- model_counts[model_counts$binary_celltype == 1, ]
  test_counts <- test_counts[test_counts$binary_celltype == 1, ]
  pred <- pred[pred$actual == 1, ]
  # from the test species counts keep only test instances
  # test_counts <- test_counts[rownames(pred),]
  # add perdition and alignment information to test set
  test_counts <- merge(test_counts, pred[,"predicted",drop=F], by =0)
  # return(test_counts)
  test_counts$align <- ifelse(test_counts$Row.names %in% test_align, 1, 0)
 
  # order by alignment and by prediction
  test_counts <- test_counts[with(test_counts, order(align, predicted)),]
  rownames(test_counts) <- test_counts$Row.names
  test_counts$Row.names <- NULL
 
  shuffled_counts <- merge(shuffled_counts, shuffled_pred[,"predicted",drop=F], by =0)
  # return(test_counts)
  shuffled_counts$align <- ifelse(shuffled_counts$Row.names %in% shuffled_align, 1, 0)
  shuffled_counts <- shuffled_counts[with(shuffled_counts, order(align, predicted)),]
  rownames(shuffled_counts) <- shuffled_counts$Row.names
  shuffled_counts$Row.names <- NULL
 
  # add missing motifs from test and shuffled sets
  test_counts <- add_missing_motifs(top_motifs, test_counts)
  shuffled_counts <- add_missing_motifs(top_motifs, shuffled_counts)
 
  ## sample a similar number of CREs from each dataset
  n <- min(c(nrow(model_counts), nrow(test_counts), nrow(shuffled_counts)))
  print(paste("Selecting", n, "CREs from each dataset"))
  # sample from each dataset
  set.seed(1)
  model_sample <- sample(x = rownames(model_counts), size = n, replace = F)
  set.seed(1)
  test_sample <- sample(x = rownames(test_counts), size = n, replace = F)
  set.seed(1)
  shuffled_sample <- sample(x = rownames(shuffled_counts), size = n, replace = F)
 
  model_counts <- model_counts[rownames(model_counts) %in% model_sample,]
  test_counts <- test_counts[rownames(test_counts) %in% test_sample,]
  shuffled_counts <- shuffled_counts[rownames(shuffled_counts) %in% shuffled_sample,]
 
  model_annot <- data.frame(dataset = rep(model_sp, nrow(model_counts))
                            , prediction = rep(model_sp, nrow(model_counts))
                            , align = rep(model_sp, nrow(model_counts)))
  test_annot <- data.frame(dataset = test_sp
                           , prediction = test_counts$predicted
                           , align = test_counts$align)
  shuffled_annot <- data.frame(dataset = paste0("shuffled_", test_sp)
                               , prediction = shuffled_counts$predicted
                               , align = shuffled_counts$align)
  rownames(model_annot) <- rownames(model_counts)
  rownames(test_annot) <- rownames(test_counts)
  rownames(shuffled_annot) <- rownames(shuffled_counts)
 
  # keep only top motifs
  model_counts <- model_counts[,top_motifs]
  test_counts <- test_counts[,top_motifs]
  shuffled_counts <- shuffled_counts[,top_motifs]
 
  annot_row <- rbind(model_annot, test_annot, shuffled_annot)
  annot_row$dataset <- factor(annot_row$dataset)
  annot_row$prediction <- factor(annot_row$prediction)
  annot_row$align <- factor(annot_row$align)
  all_counts <- rbind(model_counts, test_counts, shuffled_counts)
  # noramlize motif counts
  all_counts <- normalize_counts(all_counts)
 
  my_palette <- colorRampPalette(c("white", "red"))
 
  pheatmap(all_counts, annotation_row = annot_row
           , cluster_rows = F, cluster_cols = T, filename = out_name
           , show_rownames =F, color = my_palette(50))
}

setwd("/xgb/cross_sp")

human_xgb_path <- "/xgb/fetal/"
mouse_xgb_path <- "/xgb/e8.25/bin_500bp/"
human_data_path <- "/data/human/fetal/"
mouse_data_path <- "/data/mouse/e8.25/"

cm_mouse <- read.table(file = paste0(mouse_xgb_path, "cardiom_500bp_vert5.0_qval0.5_578T_shap_rank.txt")
                       , header = T, stringsAsFactors = F, sep ='\t')
er_mouse <- read.table(file = paste0(mouse_xgb_path, "erythroid_500bp_vert5.0_qval0.5_20T_shap_rank.txt")
                       , header = T, stringsAsFactors = F, sep ='\t')

cm_human <- read.table(file = paste0(human_xgb_path, "Cardiomyocytes_HEART_v_other_distal_nc_500bp_vert5.0_304T_shap_rank.txt")
                       , header = T, stringsAsFactors = F, sep ='\t')
er_human <- read.table(file = paste0(human_xgb_path, "Erythroblasts_MULTI_v_other_distal_nc_500bp_vert5.0_10T_shap_rank.txt")
                       , header = T, stringsAsFactors = F, sep ='\t')

cm_mouse_counts <- read.table(file = paste0(mouse_data_path, "cardiom_500bp_vert5.0_qval0.5.txt")
                              , header = T, stringsAsFactors = F, sep ='\t')
er_mouse_counts <- read.table(file = paste0(mouse_data_path, "erythroid_500bp_vert5.0_qval0.5.txt")
                              , header = T, stringsAsFactors = F, sep ='\t')

cm_human_counts <- read.table(file = paste0(human_data_path, "Cardiomyocytes_HEART_v_other_distal_nc_500bp_vert5.0.txt")
                              , header = T, stringsAsFactors = F, sep ='\t')
er_human_counts <- read.table(file = paste0(human_data_path, "Erythroblasts_MULTI_v_other_distal_nc_500bp_vert5.0.txt")
                              , header = T, stringsAsFactors = F, sep ='\t')

# random sequences motif counts
shuff_path <- "/data/random/"
human_cm_shuff <- read.table(file = paste0(shuff_path, "human_cm_counts.txt")
                             , sep ="\t", stringsAsFactors = F, header = T)
human_er_shuff <- read.table(file = paste0(shuff_path, "human_er_counts.txt")
                             , sep ="\t", stringsAsFactors = F, header = T)
mouse_cm_shuff <- read.table(file = paste0(shuff_path, "mouse_cm_counts.txt")
                             , sep ="\t", stringsAsFactors = F, header = T)
mouse_er_shuff <- read.table(file = paste0(shuff_path, "mouse_er_counts.txt")
                             , sep ="\t", stringsAsFactors = F, header = T)

# read cross - species predictions
mouse_cm_pred <- read.table(file = "mouse_cm_vs_other_byHumanFetal_vert5.0_q0.5_500bp.txt"
                            , header = T, stringsAsFactors = F)
mouse_er_pred <- read.table(file = "mouse_eryth_vs_other_byFetalErythroblasts_vert5.0_q0.5_500bp.txt"
                            , header = T, stringsAsFactors = F)
human_cm_pred <- read.table(file = "human_cm_vs_other_byMouseE8.25_vert5.0_q0.5_500bp.txt"
                            , header = T, stringsAsFactors = F)
human_er_pred <- read.table(file = "human_Erythroblasts_vs_other_byMouseE8.25Eryth_vert5.0_q0.5_500bp.txt"
                            , header = T, stringsAsFactors = F)

human_CRE_path <- "/data/human/fetal/"
mouse_CRE_path <- "/mouse_data/"

human_cm_CRE <- read.table(file = paste0(human_CRE_path
                                         , "Cardiomyocytes_HEART_distal_nc_500bp_0.bed")
                           , header = F, stringsAsFactors = F, sep ='\t')
human_er_CRE <- read.table(file = paste0(human_CRE_path
                                         , "Erythroblasts_MULTI_distal_nc_500bp_0.bed")
                           , header = F, stringsAsFactors = F, sep ='\t')
mouse_cm_CRE <- read.table(file = paste0(mouse_CRE_path
                                         , "cardiom_distal_nc_500bp_0.bed")
                           , header = F, stringsAsFactors = F, sep ='\t')
mouse_er_CRE <- read.table(file = paste0(mouse_CRE_path
                                         , "erythroid_distal_nc_500bp_0.bed")
                           , header = F, stringsAsFactors = F, sep ='\t')

human_cm_align <- read.table(file = paste0(human_CRE_path
                                           , "Cardiomyocytes_HEART_500bp_mm10_0.6.bed")
                             , header = F, stringsAsFactors = F, sep ='\t')
human_er_align <- read.table(file = paste0(human_CRE_path
                                           , "Erythroblasts_MULTI_500bp_mm10_0.6.bed")
                             , header = F, stringsAsFactors = F, sep ='\t')
mouse_cm_align <- read.table(file = paste0(mouse_CRE_path
                                           , "mouse_cardiom_500bp_hg19_0.6.bed")
                             , header = F, stringsAsFactors = F, sep ='\t')
mouse_er_align <- read.table(file = paste0(mouse_CRE_path
                                           , "mouse_erythroid_500bp_hg19_0.6.bed")
                             , header = F, stringsAsFactors = F, sep ='\t')
# keep only CREs that align
human_cm_CRE <- human_cm_CRE[human_cm_CRE$V4 %in% human_cm_align$V4, ]
human_er_CRE <- human_er_CRE[human_er_CRE$V4 %in% human_er_align$V4, ]
mouse_cm_CRE <- mouse_cm_CRE[mouse_cm_CRE$V4 %in% mouse_cm_align$V4, ]
mouse_er_CRE <- mouse_er_CRE[mouse_er_CRE$V4 %in% mouse_er_align$V4, ]

mouse_cm_CRE$V1 <- sub("chr", "", mouse_cm_CRE$V1)
mouse_er_CRE$V1 <- sub("chr", "", mouse_er_CRE$V1)

# reading predictions on random sequences
mouse_cm_shuff_pred <- read.table(file = "mouse_cm_shuff_byHumanFetal_vert5.0_q0.5_500bp.txt"
                                  , header = T, stringsAsFactors = F)
mouse_er_shuff_pred <- read.table(file = "mouse_eryth_shuff_byFetalErythroblasts_vert5.0_q0.5_500bp.txt"
                                  , header = T, stringsAsFactors = F)
human_cm_shuff_pred <- read.table(file = "human_cm_shuff_byMouseE8.25_vert5.0_q0.5_500bp.txt"
                                  , header = T, stringsAsFactors = F)
human_er_shuff_pred <- read.table(file = "human_Erythroblasts_shuff_byMouseE8.25Eryth_vert5.0_q0.5_500bp.txt"
                                  , header = T, stringsAsFactors = F)


human_cm_shuff_align <- read.table(file = paste0(human_CRE_path
                                                 , "Cardiomyocytes_HEART_500bp_shuffled_mm10_0.6.bed")
                                   , header = F, stringsAsFactors = F, sep ='\t')
human_er_shuff_align <- read.table(file = paste0(human_CRE_path
                                                 , "Erythroblasts_MULTI_500bp_shuffled_mm10_0.6.bed")
                                   , header = F, stringsAsFactors = F, sep ='\t')
mouse_cm_shuff_align <- read.table(file = paste0(mouse_CRE_path
                                                 , "mouse_cardiom_500bp_shuffled_hg19_0.6.bed")
                                   , header = F, stringsAsFactors = F, sep ='\t')
mouse_er_shuff_align <- read.table(file = paste0(mouse_CRE_path
                                                 , "mouse_erythroid_500bp_shuffled_hg19_0.6.bed")
                                   , header = F, stringsAsFactors = F, sep ='\t')

cm_human_motifMean <- motif_mean(cm_human_counts)
er_human_motifMean <- motif_mean(er_human_counts)
cm_mouse_motifMean <- motif_mean(cm_mouse_counts)
er_mouse_motifMean <- motif_mean(er_mouse_counts)

# select top motifs
select_top <- function(motif_enrich, motif_rank){
  rownames(motif_rank) <- sub("_NA", "", rownames(motif_rank))
  #rank and select top 5 enriched in target set and 5 enriched in background
  motif_enrich <- motif_enrich[rownames(motif_rank), ]
  print(head(rownames(motif_rank), 20))
  print(head(rownames(motif_enrich), 20))
  pos <- head(rownames(motif_enrich[motif_enrich$direction ==
                                      "enriched_in_class_1", ]),5)
  neg <- head(rownames(motif_enrich[motif_enrich$direction ==
                                      "enriched_in_class_0", ]),5)
  print(pos)
  print(neg)
  return(c(pos, neg))
}

# CARDIOMYOCYTE MODELS
combine_counts(model_counts = cm_human_counts
               , test_counts = cm_mouse_counts
               , shuffled_counts = mouse_cm_shuff
               , top_motifs = select_top(cm_human_motifMean, cm_human_shap)
               , pred = mouse_cm_pred
               , test_align = with(mouse_cm_CRE, paste(V1, paste(V2, V3, sep = "-"), sep = ":"))
               , shuffled_pred = mouse_cm_shuff_pred
               , shuffled_align = mouse_cm_shuff_align$V4
               , model_sp = "human"
               , test_sp = "mouse"
               , out_name = "fetal_cm_cross_sp_heatmap_top5each.pdf")

combine_counts(model_counts = cm_mouse_counts
               , test_counts = cm_human_counts
               , shuffled_counts = human_cm_shuff
               , top_motifs = select_top(cm_mouse_motifMean, cm_mouse_shap)
               , pred = human_cm_pred
               , test_align = with(human_cm_CRE, paste(V1, paste(V2, V3, sep = "-"), sep = ":"))
               , shuffled_pred = human_cm_shuff_pred
               , shuffled_align = human_cm_shuff_align$V4
               , model_sp = "mouse"
               , test_sp = "human"
               , out_name = "mouse_cm_cross_sp_heatmap_top5each.pdf")

# ERYTHROID MODELS
combine_counts(model_counts = er_human_counts
               , test_counts = er_mouse_counts
               , shuffled_counts = mouse_er_shuff
               , top_motifs = select_top(er_human_motifMean, er_human_shap)
               , pred = mouse_er_pred
               , test_align = with(mouse_er_CRE, paste(V1, paste(V2, V3, sep = "-")
                                                       , sep = ":"))
               , shuffled_pred = mouse_er_shuff_pred
               , shuffled_align = mouse_er_shuff_align$V4
               , model_sp = "human"
               , test_sp = "mouse"
               , out_name = "fetal_er_cross_sp_heatmap_top5each.pdf")

combine_counts(model_counts = er_mouse_counts
               , test_counts = er_human_counts
               , shuffled_counts = human_er_shuff
               , top_motifs = select_top(er_mouse_motifMean, er_mouse_shap)
               , pred = human_er_pred
               , test_align = with(human_er_CRE, paste(V1, paste(V2, V3, sep = "-")
                                                       , sep = ":"))
               , shuffled_pred = human_er_shuff_pred
               , shuffled_align = human_er_shuff_align$V4
               , model_sp = "mouse"
               , test_sp = "human"
               , out_name = "mouse_er_cross_sp_heatmap_top5each.pdf")

# get TF names of top motifs
cm_human_top <- select_top(cm_human_motifMean, cm_human_shap)
cm_mouse_top <- select_top(cm_mouse_motifMean, cm_mouse_shap)
er_human_top <- select_top(er_human_motifMean, er_human_shap)
er_mouse_top <- select_top(er_mouse_motifMean, er_mouse_shap)

gimme_annot<- read.table(file = "/useful/gimmemotifs/gimme.vertebrate.v5.0.motif2factors_annot"
                         , header = T, stringsAsFactors = F, sep ='\t')
gimme_annot_unique <- unique(gimme_annot[,c("Motif", "Factor")])
gimme_agr <- aggregate(Factor ~ Motif, gimme_annot_unique, paste, sep = ",")

rownames(gimme_agr) <- gimme_agr$Motif

gimme_agr[cm_human_top, "Factor", drop = F]
# Factor
# GM.5.0.MADS_box.0002                AC002126.6, MEF2A, Mef2d, Mef2c, Mef2a, MEF2B, MEF2D, MEF2C, MYEF2
# GM.5.0.T.box.0004                                                                                Tbx20
# GM.5.0.MADS_box.0003                       Mef2a, Mef2c, Mef2d, AC002126.6, MEF2A, MEF2C, MEF2B, MEF2D
# GM.5.0.MADS_box.0020                                            AC002126.6, MEF2A, Mef2d, Mef2c, Mef2a
# GM.5.0.Homeodomain.0130                                                         NKX2-3, Nkx3-2, NKX3-2
# NA                                                                                                NULL
# GM.5.0.Myb_SANT.0017                                                                           SMARCA5
# NA.1                                                                                              NULL
# GM.5.0.RFX.0006         RFX2, RFX7, RFX6, RFX8, Rfx5, Rfx6, Rfx1, Rfx8, MXI1, RFX1, Rfx2, Rfx3, SREBF1
# GM.5.0.STAT.0025                                                                          STAT1, STAT3

gimme_agr[cm_mouse_top, "Factor", drop = F]
# Factor
# GM.5.0.T.box.0004                                                                                                                                                                                                                                                                                                           Tbx20
# GM.5.0.MADS_box.0003                                                                                                                                                                                                                                                  Mef2a, Mef2c, Mef2d, AC002126.6, MEF2A, MEF2C, MEF2B, MEF2D
# GM.5.0.Homeodomain.0054                                                                                                                                                                                                    Nkx2.5, Nkx2.2, NKX2-3, NKX2-4, NKX2-6, NKX2-8, Nkx2-1, Nkx2-2, Nkx2-4, NKX2-5, Nkx2-5, Nkx2-6, Nkx2.1
# GM.5.0.Homeodomain.0130                                                                                                                                                                                                                                                                                    NKX2-3, Nkx3-2, NKX3-2
# GM.5.0.Homeodomain.0174                                                                                                                                                                                                                                                        Meis1, Meis3, Pknox1, PKNOX2, Pknox2, Tgif1, Tgif2
# GM.5.0.C2H2_ZF.0296                                                                                                                                                                                                                                                                                                         PRDM6
# GM.5.0.Forkhead.0008    FOXA1, Foxa2, EP300, FOXA2, Foxa1, FOXA3, Foxa3, FOXM1, Foxm1, GATA3, FOXK2, FOXD4, FOXE1, FOXS1, FOXD4L1, FOXI2, FOXE3, FOXD4L3, FOXB2, FOXD4L5, FOXD4L6, KIAA0415, Foxf2, Foxq1, Foxk2, Foxf1a, Foxe3, Foxi1, Foxi2, Foxl2, Foxd4, Foxd2, Foxi3, Foxb2, Foxb1, Foxd3, Foxe1, Foxs1, RXRA, TCF12, TCF7L2
# GM.5.0.Forkhead.0013                                                            FOXC1, Foxc1, FOXD1, Foxd3, FOXK2, FOXD4, FOXE1, FOXS1, FOXD4L1, FOXI2, FOXE3, FOXD4L3, FOXB2, FOXD4L5, FOXD4L6, KIAA0415, Foxa1, Foxf2, Foxq1, Foxk2, Foxa3, Foxf1a, Foxe3, Foxi1, Foxi2, Foxl2, Foxd4, Foxd2, Foxi3, Foxb2, Foxe1, Foxs1, Foxd1
# GM.5.0.Ets.0002                                                                                                                                                                                                                                                                                  ELF5, EHF, Ehf, ELF3, Elf3, Elf1
# GM.5.0.Forkhead.0012                                                                                                                                                                                                                                                                                                        FOXO3

gimme_agr[er_human_top, "Factor", drop = F]
# Factor
# GM.5.0.GATA.0001                      CCNT2, GATA1, GATA2, Gata1, Gata, Gata2, TAL1, Tal1
# GM.5.0.GATA.0003     Gata3, GATA1, Gata5, Gata6, HDAC2, GATA3, Gata2, Gata4, GATA2, GATA6
# NA                                                                                   NULL
# GM.5.0.Forkhead.0024                                                         FOXA1, FOXA2
# GM.5.0.GATA.0005                                        GATA6, GATA1, HMGN3, Gata2, GATA2
# GM.5.0.IRF.0002                                                     IRF1, IRF, Irf1, IRF2
# GM.5.0.C2H2_ZF.0295                                                                  OSR2
# GM.5.0.C2H2_ZF.0192                                                                 ZNF41
# GM.5.0.bZIP.0020                                                             GABP, NFE2L2
# GM.5.0.bHLH.0097                                                                MAX, MXI1

gimme_agr[er_mouse_top, "Factor", drop = F]
# Factor
# GM.5.0.GATA.0001                  CCNT2, GATA1, GATA2, Gata1, Gata, Gata2, TAL1, Tal1
# NA                                                                               NULL
# GM.5.0.Runt.0001             RUNX2, ENSG00000250096, Runx2, Runx3, RUNX3, RUNX1, RUNX
# GM.5.0.Homeodomain.0047                                 SIX3, SIX6, Gata2, SIX1, SIX2
# GM.5.0.C2H2_ZF.0314                                      KLF7, Klf6, Klf5, Klf3, Klf2
# GM.5.0.Sox.0015                                                    LEF1, TCF7L1, Tcf7
# GM.5.0.TEA.0007                                                         TEAD2, BCL11A
# GM.5.0.Sox.0018                                              SOX10, Sox10, SOX9, Sox9
# GM.5.0.Sox.0019                     Sox2, SOX4, SOX11, SOX1, SOX13, Sox11, Sox4, SOX9
# GM.5.0.Nuclear_receptor.0035                                Hnf4a, Hnf4g, NR2F2, Rxra


cm_human_top
# [1] "GM.5.0.MADS_box.0002"    "GM.5.0.T.box.0004"      
# [3] "GM.5.0.MADS_box.0003"    "GM.5.0.MADS_box.0020"  
# [5] "GM.5.0.Homeodomain.0130" "GM.5.0.Unknown.0010"    
# [7] "GM.5.0.Myb_SANT.0017"    "GM.5.0.Unknown.0165"    
# [9] "GM.5.0.RFX.0006"         "GM.5.0.STAT.0025"      
cm_mouse_top
# [1] "GM.5.0.T.box.0004"       "GM.5.0.MADS_box.0003"  
# [3] "GM.5.0.Homeodomain.0054" "GM.5.0.Homeodomain.0130"
# [5] "GM.5.0.Homeodomain.0174" "GM.5.0.C2H2_ZF.0296"    
# [7] "GM.5.0.Forkhead.0008"    "GM.5.0.Forkhead.0013"  
# [9] "GM.5.0.Ets.0002"         "GM.5.0.Forkhead.0012"  
er_human_top
# [1] "GM.5.0.GATA.0001"     "GM.5.0.GATA.0003"     "GM.5.0.Unknown.0192"
# [4] "GM.5.0.Forkhead.0024" "GM.5.0.GATA.0005"     "GM.5.0.IRF.0002"    
# [7] "GM.5.0.C2H2_ZF.0295"  "GM.5.0.C2H2_ZF.0192"  "GM.5.0.bZIP.0020"    
# [10] "GM.5.0.bHLH.0097"    
er_mouse_top
# [1] "GM.5.0.GATA.0001"             "GM.5.0.Unknown.0192"        
# [3] "GM.5.0.Runt.0001"             "GM.5.0.Homeodomain.0047"    
# [5] "GM.5.0.C2H2_ZF.0314"          "GM.5.0.Sox.0015"            
# [7] "GM.5.0.TEA.0007"              "GM.5.0.Sox.0018"            
# [9] "GM.5.0.Sox.0019"              "GM.5.0.Nuclear_receptor.0035"
