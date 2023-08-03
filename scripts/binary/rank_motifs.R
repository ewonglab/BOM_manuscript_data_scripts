params <- list(
  shap_file = NULL
  , out_file = NULL
  , rank_type = "sum"
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
  cat("Usage: Rscript rank_motifs.R [parameters]\n")
  cat("\n")
  cat("Parameters:\n")
  cat("  --shap_file=<file>       Path to the SHAP values file\n")
  cat("  --out_file=<file>      Path to save the motifs ranked by SHAP\n")
  cat("  --rank_type=<file>      Rank type. Either 'sum' or 'mean' \n")
  cat("\n")
  q("no", status = 1)
} 


shap_file <- params$shap_file
out_file <- params$out_file
rank_type <- params$rank_type

rank_shap <- function(shap.df, type){
  if(!type %in% c("sum", "mean")){
    stop(paste("Invalid ranking type:", type))
  }
  if(type == "sum"){
    shap_per_motif <- apply(shap.df, 2, function(x) sum(abs(x)))
  }
  if(type == "mean"){
    shap_per_motif <- apply(shap.df, 2, function(x) mean(abs(x)))
  }
  
  shap_per_motif <- shap_per_motif[order(-shap_per_motif)]
  return(shap_per_motif)
  }


# READING SHAP VALUES
cat("Reading SHAP values...\n")

shap <- read.table(file = shap_file, header = T, stringsAsFactors = F, sep ='\t', row.names = 1)

cat("Ranking motifs by SHAP...\n")

ranked_motifs <- rank_shap(shap, rank_type)
ranked_motifs <- as.data.frame(ranked_motifs)

if(rank_type == "sum"){
  colnames(ranked_motifs) <- "sum_abs_SHAP"
}
if(rank_type == "mean"){
  colnames(ranked_motifs) <- "mean_abs_SHAP"
}


cat("Saving ranked motifs...\n")

write.table(x = ranked_motifs, file = out_file, quote = F)

cat("Done\n")


