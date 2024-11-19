suppressMessages({
  library(data.table)
  library(tidyr)
  library(dplyr)
})

data_path <- "trimmed600bp/atha_CisBP2.0/"
qval_thresh <- 0.5
out_filename <- "atha_CisBP2.0_4CellT_600bp_q0.5.txt"

celltypes <- c("c.e_precursor", "endodermis_2", "endodermis_3", "stele_1.xylem")
counts_files <- paste0(data_path, celltypes, "/fimo.tsv")

counts <- lapply(counts_files, fread)
counts <- (lapply(counts, as.data.frame))

# subset by q-value (column 9)
counts <- (lapply(counts, function(x) x[x[,9] <= qval_thresh,]))

# add cell type
for(i in 1:length(counts)){
  counts[[i]]$celltype <- celltypes[i]
}

# list to data frame
counts <- do.call("rbind", counts)
counts$motif.id <- paste(counts$motif_id, counts$motif_alt_id, sep='_')
final_dataset <- as.data.frame(table(counts[,c("motif.id", "sequence_name")]))
final_dataset <- tidyr::spread(final_dataset, motif.id, Freq)
final_dataset <- merge(final_dataset, unique(counts[,c("sequence_name", "celltype")])
                       , by = "sequence_name", all.x = T)

### cell type as numeric
unique_ct <- unique(final_dataset$celltype)
unique_ct <- unique_ct[order(unique_ct)]
ct_annot <- data.frame(celltype = unique_ct
                       , celltype_numeric= 0:(length(unique_ct)-1))
ct_annot
#        celltype celltype_numeric
# 1 c.e_precursor                0
# 2  endodermis_2                1
# 3  endodermis_3                2
# 4 stele_1.xylem                3

final_dataset <- merge(final_dataset, ct_annot, by = "celltype", all.x=T)

rownames(final_dataset) <- final_dataset$sequence_name
final_dataset$sequence_name <- NULL

final_dataset <- final_dataset[,c(colnames(final_dataset)[!colnames(final_dataset) %in% c("celltype", "celltype_numeric")]
                                  ,"celltype", "celltype_numeric")] 

# save dataset
write.table(x = final_dataset, file = out_filename, quote = F, sep = '\t')

print("Content of output table: ")
print(as.data.frame(table(final_dataset$celltype)))
#            Var1 Freq
# 1 c.e_precursor  216
# 2  endodermis_2  127
# 3  endodermis_3  127
# 4 stele_1.xylem  139
