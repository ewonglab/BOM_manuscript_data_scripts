# Load packages
required_packages <- c("data.table","tidyr","dplyr")

load_required_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is not installed. Please install it before proceeding."))
    } else {
      suppressMessages(library(pkg, character.only = TRUE))
    }
  }
}

load_required_packages(required_packages)

# Read parameters

params <- list(
  target_ct = NULL,
  data_path = NULL,
  qval_thresh = NULL,
  out_filename = NULL)


display_help <- FALSE

args <- commandArgs(trailingOnly = TRUE)

for (arg in args) {
  param <- strsplit(arg, "=", fixed = TRUE)[[1]]
  if (param[1] == "--help") {
    display_help <- TRUE
    break
  }
  param[1] <- sub("--", "", param[1])
  if (param[1] %in% names(params)) {
    params[[param[1]]] <- param[2]
  }
}

if (display_help) {
  
  cat("Usage: Rscript matrix_for_binary_model.R [options]\n")
  cat("\n")
  cat("Options:\n")
  cat("  --target_ct=<target_cell_type>    Name of the target cell type\n")
  cat("  --data_path=<data_directory>     Path to the data directory\n")
  cat("  --qval_thresh=<qvalue_threshold> Q-value threshold for filtering\n")
  cat("  --out_filename=<output_filename> Name for the output file\n")
  cat("  --help                           Display this help message\n")
  
  q("no", status = 0)
} else {
  target_ct <- params$target_ct
  data_path <- params$data_path
  qval_thresh <- as.numeric(params$qval_thresh)
  out_filename <- params$out_filename
}


# list directories containing FIMO output
directories <- list.dirs(path = data_path, full.names = T, recursive=F)
#print(directories)
celltypes <- basename(directories)

counts_files <- (paste0(directories, "/fimo.tsv"))

# num_lines <- system("wc -l < file_path", intern = TRUE)
# num_lines_to_read <- num_lines - num_lines_to_skip
# data <- fread(file_path, skip = num_lines_to_read)

# Read results of motif search
cat(paste("Reading data...\n"))
suppressWarnings({
  counts <- lapply(counts_files, fread)
})

#counts <- suppressMessages(lapply(counts_files, fread))
names(counts) <- celltypes

counts <- (lapply(counts, as.data.frame))
counts <- (lapply(counts, function(x) x[x[,9] <= qval_thresh,]))

## count the number of unique CREs per condition
n_enh_by_celltype <- unlist(lapply(counts, function(x) length(unique(x$sequence_name))))
names(n_enh_by_celltype) <- celltypes

## define the number of CREs from target and background conditions
if(! target_ct %in% names(n_enh_by_celltype)){
  stop(paste(target_ct, "is not among the contexts provided. Please check the spelling and case."))
}

n <- n_enh_by_celltype[target_ct]
other_n <- round(n/(length(celltypes) -1))

#### ***************************************************************************** ###
#       check that there is enough enhancers of every cell type to select          ###
#                         the necessary contribution                               ###
#### ***************************************************************************** ###

if(other_n > min(n_enh_by_celltype)){
  print("Warning: reducing the number of enhancers selected from each cell type")
  other_n <- min(n_enh_by_celltype)
  n <- other_n * (length(celltypes) -1)
  tmp <- counts[[target_ct]]
  set.seed(123)
  pos_sample <- sample(x = unique(tmp$sequence_name), size = n, replace = F)
  counts[[target_ct]] <- tmp[tmp$sequence_name %in% pos_sample,]
}

## *****************************************###
###      PREPARING NEGATIVE DATASET         ###
## *****************************************###

negative_set <- lapply(celltypes[celltypes != target_ct], function(ct) {
  tmp <- counts[[ct]]
  set.seed(123)
  ct_sample <- sample(unique(tmp$sequence_name), size = other_n, replace = FALSE)
  tmp[tmp$sequence_name %in% ct_sample, ]
})

# collapse segative and positive datasets as a single data frame for each
# add cell types to the negative sets
names(negative_set) <- celltypes[celltypes!=target_ct]

for (i in seq_along(negative_set)) {
  negative_set[[i]]$celltype <- celltypes[celltypes != target_ct][i]
}

negative_set.df <- do.call("rbind", negative_set)

negative_set.df <- negative_set.df[,c("motif_id", "sequence_name", "celltype")]
tmp <- as.data.frame(table(negative_set.df[,1:2]))

tmp <- tidyr::spread(tmp, motif_id, Freq)
tmp <- merge(tmp, unique(negative_set.df[,2:3]), by="sequence_name")
rownames(tmp) <- tmp$sequence_name
tmp$sequence_name <- NULL

## *****************************************###
###      PREPARING POSITIVE DATASET         ###
## *****************************************###

positive <- counts[[target_ct]]
positive$celltype <- target_ct
positive <- positive[,c("motif_id", "sequence_name", "celltype")]
tmp2 <- as.data.frame(table(positive[,1:2]))
tmp2 <- tidyr::spread(tmp2, motif_id, Freq)
tmp2 <- merge(tmp2, unique(positive[,2:3]), by="sequence_name")
#enhancer IDs as rownames
rownames(tmp2) <- tmp2$sequence_name
tmp2$sequence_name <- NULL

## ****************************************************************###
## make a single dataframe with the positive and negative datasets  ##
## ****************************************************************###

final_set <- dplyr::bind_rows(tmp, tmp2)
final_set[is.na(final_set)] <- 0
final_set <- final_set[,c(colnames(final_set)[colnames(final_set)!="celltype"],"celltype")]
final_set$binary_celltype <- ifelse(final_set$celltype==target_ct, 1, 0)

## ****************************************************************###

# save dataset
cat(paste("Saving matrix of motif counts for binary classification...\n"))
write.table(x = final_set, file = out_filename, quote = F, sep = '\t')

cat("Content of output table:\n")
out_content <- as.data.frame(table(final_set$binary_celltype))
colnames(out_content) <- c("Context", "Number of elements")
print(out_content)
cat('\n')

## ****************************************************************###

