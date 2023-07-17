# Read parameters

params <- list(
  input_bed = NULL,
  annot = NULL,
  chr_sizes = NULL,
  u = NULL,
  d = NULL,
  nbp = NULL,
  keep_proximal = FALSE,
  remove_proximal = FALSE,
  non_exonic = FALSE,
  out_bed = NULL
)

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
  cat("Usage: filter_CREs.R [parameters]\n")
  cat("\n")
  cat("Parameters:\n")
  cat("--input_bed=<file> input BED file\n")
  cat("--annot=<file> assembly annotation file (.gtf) \n")
  cat("--chr_sizes=<file> File with chromosome length values \n")
  cat("--u=<integer>   Number of base pairs upstream of TSS for the definition of proximal regions\n")
  cat("--d=<integer>   Number of base pairs downstream of TSS for the definition of proximal regions\n")
  cat("--nbp=<integer> Number of central base pairs in adjusted CREs \n")
  cat("--out_bed=<file> output BED file\n")
  cat("--keep_proximal=<logical> whether only proximal regions to TSS should be retained (default: FALSE)\n")
  cat("--remove_proximal=<logical> whether proximal regions to TSS should be removed (default: FALSE)\n")
  cat("--non_exonic=<logical> whether regions overlapping exons should be removed (default: FALSE)\n")
  cat("--help                           Display this help message\n")
  cat("\n")
  q("no", status = 0)
} else {
  input_bed <- params$input_bed
  annot <- params$annot
  chr_sizes <- params$chr_sizes

  if (!is.null(params$u)) {
    u <- as.integer(params$u)
  }

  if (!is.null(params$d)) {
    d <- as.integer(params$d)
  }

  if (!is.null(params$nbp)) {
    nbp <- as.integer(params$nbp)
  }

  keep_proximal <- as.logical(params$keep_proximal)
  remove_proximal <- as.logical(params$remove_proximal)
  non_exonic <- as.logical(params$non_exonic)
  out_bed <- params$out_bed
}

# Load packages

cat("Loading packages...\n")
required_packages <- c(
#  "rlang",
  "GenomicFeatures",
  "GenomicRanges"
)

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

# Function to trim (or extend) CREs 
adjust_CREs <- function(x, N){
  
  x$width <- with(x, V3-V2)
  #x$centre <- with(x, V2 + (width/2))
  x$centre <- with(x, V2 + round(width/2))
  #x$start <- ceiling(x$centre - N/2)
  x$start <- round(x$centre - N/2)
  #x$end <- ceiling(x$centre + N/2)
  x$end <- round(x$centre + N/2)

  # add chromosome size
  x$ord <- 1:nrow(x)
  x <- merge(x, chrom_sizes, by.x="V1", by.y="chr")

  # Look for CREs at the extremes of chromosomes
  if(any(x$start < 0) | any(x$end > x$chr_size)){
    warning(paste("Warning: some CREs are shorter than ", paste0(N, "bp")
                  , "because they are at the extremes of chromosomes"))
    if(any(x$start < 0)){
      # special case 1: negative start position 
      x$start <- ifelse(x$start < 0, 0, x$start)
    }
    if(any(x$end > x$chr_size)){
      # special case 2: CRE end position exceeds chromosome length
      x$end <- with(x, ifelse(end > chr_size, chr_size, end))
    }
  }
  
  #return to original order
  x <- x[with(x, order(ord)), ]
  x <- x[,colnames(x)[!colnames(x) %in% c("chr_size", "centre", "width", "V2", "V3", "ord")]]
  x <- x[,c("V1", "start", "end", setdiff(colnames(x), c("V1", "start", "end")))]
  return(x)
  
}


cat("Reading CREs...\n")
cres <- read.table(file = input_bed, header = F, stringsAsFactors = F, sep = '\t')
cres_gr <- with(cres, GRanges(V1, IRanges(V2+1, V3)))

if(keep_proximal | remove_proximal | non_exonic){

  cat("Reading genome annotation...\n")
  txdb_obj <- makeTxDbFromGFF(file = annot, format = "gtf")

  if(non_exonic){
    cat("Removing exonic regions...\n")
    exons <- exons(txdb_obj)
    x <- as.data.frame(findOverlaps(cres_gr, exons))
    if(nrow(x) > 0){
      cres <- cres[-unique(x$queryHits),]
    } 
  }
  
  if(keep_proximal & remove_proximal){
    stop(paste("'keep_proximal' and 'remove_proximal' are mutually exclussive"))
  }
  if(keep_proximal){
    cat("Keeping proximal regions to TSSs...\n")
    proximal <- promoters(x = txdb_obj, upstream = u, downstream = d)
    cres_gr <- with(cres, GRanges(V1, IRanges(V2+1, V3)))
    x <- as.data.frame(findOverlaps(cres_gr, proximal))
    if(nrow(x) > 0){
      cres <- cres[unique(x$queryHits),]
    } else {
      stop(paste("No remaining CREs"))
    }
  }
  if(remove_proximal){
    cat("Removing proximal regions to TSSs...\n")
    proximal <- promoters(x = txdb_obj, upstream = u, downstream = d)
    cres_gr <- with(cres, GRanges(V1, IRanges(V2+1, V3)))
    x <- as.data.frame(findOverlaps(cres_gr, proximal))
    if(nrow(x) > 0){
      cres <- cres[-unique(x$queryHits),]
    } 
    
  }
}

#if((u == NULL & d != NULL) | (d == NULL & u != NULL)){
#  stop(paste("u AND d should be provided to adjust CREs"))
if (exists("nbp") && !is.null(nbp)) {
  cat("Reading chromosome sizes...\n")
  chrom_sizes <- read.table(file = chr_sizes, header = F, stringsAsFactors = F, sep ='\t')
  colnames(chrom_sizes) <- c("chr", "chr_size")
  
  # Adjust CREs
  cat("Adjusting CRE length...\n")
  cres <- adjust_CREs(cres, nbp)
}


# Save filtered regions
cat("Saving CREs...\n")

if(nrow(cres) == 0){
  stop(paste("No remaining CREs"))
}

write.table(x = cres, file = out_bed, quote = F, col.names = F
            , row.names = F, sep ='\t')


cat(paste("Done\n"))


