# Load required packages

required_packages <- c("rsample","xgboost")

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
  input_data = NULL
  , xgb_model = NULL
  , pred = "predictions.txt"
  , training_set = NULL
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
  cat("Usage: Rscript xgboost_predictions.R [parameters]\n")
  cat("\n")
  cat("Parameters:\n")
  cat("--input_data=<file>       Path to the input data file\n")
  cat("--xgb_model=<file>        Path to the xgboost model file\n")
  cat("--predictions=<file>      Path to save the predicted values\n")
  cat("--training_set=<file>     Path to save the training set (optional)\n")
  cat("--help                    Display this help message\n")
  cat("\n")
  q("no", status = 0)
  } else {
  input_data <- params$input_data
  xgb_model <- params$xgb_model
  predictions <- params$predictions
  training_set <- params$training_set
}

add.missing.vars_xgb <- function(xgb.model , data){
  missing.vars <- setdiff(xgb.model$feature_names, colnames(data))
  data[, missing.vars] <- lapply(missing.vars, function(var) as.integer(0))
  return(data)}


# READING XGBOOST MODEL

cat("Reading model...\n")

#xgb <- xgb.load(xgb_model)
xgb <- readRDS(xgb_model)
xgb.save(xgb, gsub(".rds", ".bin", xgb_model))

cat(paste("Best tree:", xgb$best_iteration, "\n"))

# Reading table of motif counts

cat("Reading data...\n")
counts.tab <- read.table(file = input_data, header =T, stringsAsFactors = F, sep = '\t')
counts.tab$celltype <- NULL

counts.tab.NAs <- sapply(counts.tab, function(x) sum(is.na(x)))

if(any(counts.tab.NAs) > 0){
  cat("There are NAs in input data...\n")  
}

# Binary label as numeric
counts.tab$binary_celltype <-as.numeric(counts.tab$binary_celltype)

# Split dataset into training, validation and test sets
set.seed(123)
motifs_split <- initial_split(counts.tab, prop = .6)
motifs_train <- training(motifs_split)
motifs_test <- testing(motifs_split)

set.seed(123)
motifs_split2 <- initial_split(motifs_test, prop = .5)
motifs_val <- training(motifs_split2)
motifs_test <- testing(motifs_split2)

# Removing non-variable motifs from training set
motifs_train.sd <- apply(motifs_train, 2, sd)
motifs_train <- motifs_train[, names(which(motifs_train.sd != 0))]

# Save training set if a file name is provided
if (!is.null(training_set)) {
  cat("Saving training set...\n")
  write.table(x = motifs_train, file = training_set, quote = FALSE, sep ='\t')
}

motifs_test <- add.missing.vars_xgb(xgb, motifs_test)
test_labels <- motifs_test$binary_celltype
motifs_test <- motifs_test[,xgb$feature_names]  

set.seed(123)
y_pred <- predict(xgb, data.matrix(motifs_test), type="response")
predicted.class <- y_pred > 0.5 
predicted.class <- gsub("TRUE", 1, predicted.class)
predicted.class <- gsub("FALSE", 0, predicted.class)

actual.vs.predicted <- data.frame(true_class = test_labels
                                  , predicted_class = predicted.class, prob = y_pred
                                  , stringsAsFactors = F)

rownames(actual.vs.predicted) <- rownames(motifs_test)

cat(paste("Saving predicted values...\n"))

write.table(x = actual.vs.predicted, file = params$pred, quote = F)

cat("Done\n")

