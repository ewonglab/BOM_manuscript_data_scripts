# Example: Rscript get_nonOverlapping_motifs.R cardiom cardiom_vert5.0_nonOv.txt

library(data.table)
library(tidyr)
library(dplyr)

get_non_overlap_fimo <- function(one_enh){
  one_enh$motif_coor <- with(one_enh, paste(start, stop, sep = "_"))
  one_enh <- one_enh[with(one_enh, order(sequence_name, start
                                         , stop, -score)),]
  one_enh <- one_enh[!duplicated(one_enh[,c("motif_coor", "sequence_name")]), ]
  one_enh$keep <- ""
  for(i in 1:nrow(one_enh)){
    # if(i%%50 ==0){print(i)}
    if(one_enh[i,"keep"]!="no"){
      if(one_enh[i,"keep"]=="yes"){
        #the motif has been evaluated before and I just need to check if it overlaps
        # with motifs not evaluated yet
        #I also have to avoid taking previous selected motifs
        #I can just take i and do r bind with motifs not evaluated yet
        tmp <- rbind(one_enh[i,], one_enh[one_enh$keep == "",])
        motifs_sub <- tmp[tmp$start <= one_enh[i,"stop"], "motif_coor"]}
      if(one_enh[i,"keep"]==""){
        #the motif has not been tested before, so it does not overlap with previous motifs
        motifs_sub <- one_enh[one_enh$start <= one_enh[i,"stop"] &
                                one_enh$keep == ""
                              , "motif_coor"]}
      #find the motif with the best score
      best_score <- (max(one_enh[one_enh$motif_coor %in% motifs_sub, "score"]))
      best_motif <- one_enh[one_enh$motif_coor %in% motifs_sub &
                              one_enh$score == best_score, "motif_coor"]
      if(length(best_motif) > 1){
        set.seed(1)
        best_motif <- sample(x = best_motif, size = 1, replace = F)
      }
      #keep motif(s) with the best score
      one_enh[one_enh$motif_coor == best_motif, "keep"] <- "yes"
      one_enh[one_enh$motif_coor %in% motifs_sub &
                one_enh$motif_coor != best_motif, "keep"] <- "no"}
  }
  return(one_enh[one_enh$keep=="yes",])}

setwd("/gimme")
args <- commandArgs(trailingOnly=TRUE)
#species <- args[1]
target_ct <- args[1]
out_filename <- args[2]

data_path <- paste0(target_ct, "_Gimme_vertv5.0", "/fimo.tsv")
counts_df <- fread(data_path)
counts_df <- as.data.frame(counts_df)
counts_df <- counts_df[counts_df$sequence_name!="",]
counts_df$enh_start <- gsub("-.*", "", gsub(".*:", "", counts_df$sequence_name))
counts_df$enh_start <- as.integer(counts_df$enh_start)
counts_df$abs_start <- with(counts_df, enh_start + start)
counts_df$abs_stop <- with(counts_df, enh_start + stop)

uniq_enhancers <- unique(counts_df$sequence_name)
# x <- motifs_enh_df[motifs_enh_df$enh_id=="chr1_106199556_106200149",]#326  12
nonov_enh <- lapply(uniq_enhancers, function(x)
  get_non_overlap_fimo(counts_df[counts_df$sequence_name == x,]))

nonov_enh.df <- do.call("rbind", nonov_enh)
write.table(x = nonov_enh.df, file = out_filename, quote = F, sep = '\t')
