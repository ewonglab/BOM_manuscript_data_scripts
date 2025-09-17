suppressMessages(library("Biostrings"))
args <- commandArgs(trailingOnly=TRUE)

seq2kmers <- function(sequence){
  sequence <- toupper(sequence)
  tmp <- sapply(seq(from=1, to=nchar(sequence), by=1), function(x) substr(sequence, x, x+(k-1)))
  tmp <- tmp[which(nchar(tmp)==k)]
  tmp <- paste(tmp, collapse = " ")
  return(tmp)}

celltype <- args[1]
k <- as.integer(args[2])
input_path <- args[3]
out_path <- args[4]

pos_train <- paste0(celltype, "_positive_500bp_train.fa")
neg_train <- paste0(celltype, "_negative_500bp_train.fa")
pos_val <- paste0(celltype, "_positive_500bp_val.fa")
neg_val <- paste0(celltype, "_negative_500bp_val.fa")
pos_test <- paste0(celltype, "_positive_500bp_test.fa")
neg_test <- paste0(celltype, "_negative_500bp_test.fa")

setwd(input_path)

pos_train.bs <- readDNAStringSet(pos_train)
neg_train.bs <- readDNAStringSet(neg_train)
pos_val.bs <- readDNAStringSet(pos_val)
neg_val.bs <- readDNAStringSet(neg_val)
pos_test.bs <- readDNAStringSet(pos_test)
neg_test.bs <- readDNAStringSet(neg_test)
### get kmers: iterate over the sequences of every set and get "k-mers", remove those kmers 
# shorter than k 

pos_train.kmers <- lapply(1:length(pos_train.bs)
                          , function(id) seq2kmers(as.data.frame(pos_train.bs[id])$x))

neg_train.kmers <- lapply(1:length(neg_train.bs)
                          , function(id) seq2kmers(as.data.frame(neg_train.bs[id])$x))

pos_val.kmers <- lapply(1:length(pos_val.bs)
                          , function(id) seq2kmers(as.data.frame(pos_val.bs[id])$x))

neg_val.kmers <- lapply(1:length(neg_val.bs)
                          , function(id) seq2kmers(as.data.frame(neg_val.bs[id])$x))

pos_test.kmers <- lapply(1:length(pos_test.bs)
                         , function(id) seq2kmers(as.data.frame(pos_test.bs[id])$x))

neg_test.kmers <- lapply(1:length(neg_test.bs)
                         , function(id) seq2kmers(as.data.frame(neg_test.bs[id])$x))


## make a dataframe of kmers and label
# label 1 for positive set
# label 0 for negative set

pos_train.df <- data.frame(sequence=unlist(pos_train.kmers), label=1)
neg_train.df <- data.frame(sequence=unlist(neg_train.kmers), label=0)
train_set.df <- rbind(pos_train.df, neg_train.df)

pos_val.df <- data.frame(sequence=unlist(pos_val.kmers), label=1)
neg_val.df <- data.frame(sequence=unlist(neg_val.kmers), label=0)
val_set.df <- rbind(pos_val.df, neg_val.df)

pos_test.df <- data.frame(sequence=unlist(pos_test.kmers), label=1)
neg_test.df <- data.frame(sequence=unlist(neg_test.kmers), label=0)
test_set.df <- rbind(pos_test.df, neg_test.df)

## create a new directory where training and test sets will be saved with proper names
# for model fine tunning

system(paste("mkdir", out_path))
setwd(out_path)
print(paste("kmers saved to", out_path))
write.table(x = train_set.df, file = "train.tsv", row.names = F, sep = '\t', quote = F)
write.table(x = val_set.df, file = "val.tsv", row.names = F, sep = '\t', quote = F)
write.table(x = test_set.df, file = "dev.tsv", row.names = F, sep = '\t', quote = F)
