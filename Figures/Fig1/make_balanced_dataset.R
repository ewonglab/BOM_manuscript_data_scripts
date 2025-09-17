library(data.table)
library(tidyr)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
target_ct <- args[1]
data_path <- args[2]
qval_thresh <- as.numeric(args[3])
out_filename <- args[4]

directories <- list.dirs(path = data_path, full.names = T, recursive=F)
celltypes <- basename(directories)
counts_files <- paste0(directories, "/fimo.tsv")

counts <- lapply(counts_files, fread)
names(counts) <- celltypes

counts <- (lapply(counts, as.data.frame))
#subset by q-value (column 9)
counts <- (lapply(counts, function(x) x[x[,9] <= qval_thresh,]))

n_enh_by_celltype <- unlist(lapply(counts, function(x) length(unique(x$sequence_name))))
names(n_enh_by_celltype) <- celltypes

n <- n_enh_by_celltype[target_ct]
other_n <- round(n/(length(celltypes) -1))
                                   

if(other_n > min(n_enh_by_celltype)){
  print("Warning: reducing the number of enhancers selected from each cell type")
  other_n <- min(n_enh_by_celltype)
  n <- other_n * (length(celltypes) -1)
  # selecting positive set at random
  #set.seed(123)
  tmp <- counts[[target_ct]]
  set.seed(123)
  pos_sample <- sample(x = unique(tmp$sequence_name), size = n, replace = F)
  counts[[target_ct]] <- tmp[tmp$sequence_name %in% pos_sample,]
  }


negative_set <- list()
for(celltype in celltypes[celltypes!=target_ct]){
  tmp <- counts[[celltype]]
  set.seed(123)
  ct_sample <- sample(x = unique(tmp$sequence_name), size = other_n, replace = F)
  negative_set[[length(negative_set)+1]] <- tmp[tmp$sequence_name %in% ct_sample,]}

names(negative_set) <- celltypes[celltypes!=target_ct]

for(i in 1:length(negative_set)){
  negative_set[[i]]$celltype <- names(negative_set[i])}
  
negative_set.df <- do.call("rbind", negative_set)
negative_set.df <- negative_set.df[,c("motif_id", "sequence_name", "celltype")]
#count motif instances in every enhancer
tmp <- as.data.frame(table(negative_set.df[,1:2]))


tmp <- tidyr::spread(tmp, motif_id, Freq)
tmp <- merge(tmp, unique(negative_set.df[,2:3]), by="sequence_name")
rownames(tmp) <- tmp$sequence_name
tmp$sequence_name <- NULL


positive <- counts[[target_ct]]
positive$celltype <- target_ct
positive <- positive[,c("motif_id", "sequence_name", "celltype")]
tmp2 <- as.data.frame(table(positive[,1:2]))
tmp2 <- tidyr::spread(tmp2, motif_id, Freq)
tmp2 <- merge(tmp2, unique(positive[,2:3]), by="sequence_name")
#enhancer IDs as rownames
rownames(tmp2) <- tmp2$sequence_name
tmp2$sequence_name <- NULL

final_set <- dplyr::bind_rows(tmp, tmp2)
final_set[is.na(final_set)] <- 0
final_set <- final_set[,c(colnames(final_set)[colnames(final_set)!="celltype"],"celltype")]
final_set$binary_celltype <- ifelse(final_set$celltype==target_ct, 1, 0)
print("Content of output table: ")
print(table(final_set$binary_celltype))

write.table(x = final_set, file = out_filename, quote = F, sep = '\t')
