# Load packages

required_packages <- c("rsample","xgboost")

load_required_packages <- function(packages) { for (pkg in packages) { if (!requireNamespace(pkg, quietly = TRUE)) { 
      stop(paste("Package", pkg, "is not installed. Please install it before proceeding."))
    } else {
      suppressMessages(library(pkg, character.only = TRUE))
    }
  }
}


load_required_packages(required_packages)

# Read parameters

params <- list(input_data = NULL
               , data = NULL
               , nrounds = 10000
               , eta = 0.01
               , max_depth = 6
               , subsample = 0.5
               , colsample_bytree = 0.5
               , objective = "binary:logistic"
               , watchlist = NULL
               , early_stopping_rounds = NULL
               , nthread = 1
               , eval_metric = "error"
               , maximize = F
               , params = list()
               , feval = NULL
               , verbose = 1
               , print_every_n = 1L
               , save_period = NULL
               , save_name = "xgboost.model"
               , xgb_model = NULL
               , callbacks = list())


display_help <- FALSE

# Parse command-line arguments
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

### correct data type for parameters
params$nrounds <- as.integer(params$nrounds)

params$max_depth <- as.integer(params$max_depth)

if(!is.null(params$early_stopping_rounds)){
  params$early_stopping_rounds <- as.integer(params$early_stopping_rounds)
}

params$nthread <- as.integer(params$nthread)

if(!is.null(params$save_period)){
  params$save_period <- as.integer(params$save_period)
}

if(!is.null(params$nthread)){
  params$nthread <- as.integer(params$nthread)
}

params$eta <- as.numeric(params$eta)
params$subsample <- as.numeric(params$subsample)
params$colsample_bytree <- as.numeric(params$colsample_bytree)
params$maximize <- as.logical(params$maximize)


if (display_help) {
  
  cat("Usage: Rscript training_binary.R [options]\n")
  
  cat("Options:\n")
  cat("--input_data=<file>\t\tPath to the input matrix of motif counts (required)\n")
  cat("--data=<data>\t\t\tThe training data (default: dtrain)\n")
  cat("--nrounds=<n>\t\t\tNumber of boosting rounds (default: 10000)\n")
  cat("--eta=<value>\t\t\tLearning rate (default: 0.01)\n")
  cat("--max_depth=<n>\t\tMaximum tree depth (default: 6)\n")
  cat("--subsample=<value>\t\tSubsample ratio of the training instances (default: 0.5)\n")
  cat("--colsample_bytree=<value>\tSubsample ratio of columns when constructing each tree (default: 0.5)\n")
  cat("--objective=<name>\t\tObjective function (default: binary:logistic)\n")
  cat("--early_stopping_rounds=<n>\tPerform early stopping if no improvement for this many rounds (default: NULL)\n")
  cat("--nthread=<n>\t\t\tNumber of parallel threads (default: 1)\n")
  cat("--eval_metric=<name>\t\tEvaluation metric (default: error)\n")
  cat("--maximize=<bool>\t\tWhether to maximize the evaluation metric (default: FALSE)\n")
  cat("--save_period=<n>\t\tSave the model for every given number of rounds (default: NULL)\n")
  cat("--save_name=<file>\t\tName of the saved model file (default: xgboost.model)\n")
  cat("--feval=<file>\t\tCustomized evaluation metric (default: NULL)\n")
  cat("--Verbose=<file>\t\tHow much details on the progress to print (default: 1)\n")
  cat("--print_every_n=<file>\t\tPrint evaluation messages each n-th iterations (default: 1)\n")
  cat("--save_name=<file>\t\tName of the saved model file (default: xgboost.model)\n")
  cat("--help                    Display this help message\n")
  cat("\n")
  
  q("no", status = 0)
} 

input_data <- params$input_data
params <- params[-match("input_data", names(params))]

# Reading table of motif counts
cat("Reading data...\n")
counts.tab <- read.table(file = input_data, header =T
                         , stringsAsFactors = F, sep = '\t')
counts.tab$celltype <- NULL

counts.tab.NAs <- sapply(counts.tab, function(x) sum(is.na(x)))
if(any(counts.tab.NAs) > 0){
 print("Warning: NAs in input data...")  
}

# Binary label as numeric
counts.tab$binary_celltype <- as.numeric(counts.tab$binary_celltype)

# Split dataset into training, validation and test sets
cat("Splitting data into training, validation and test sets...\n")
set.seed(123)
motifs_split <- initial_split(counts.tab, prop = .6)
motifs_train <- training(motifs_split)
motifs_test <- testing(motifs_split)

set.seed(123)
motifs_split2 <- initial_split(motifs_test, prop = .5)
motifs_val <- training(motifs_split2)
motifs_test <- testing(motifs_split2)

# REMOVING NON-VARIABLE MOTIFS FROM TRAINING SET
motifs_train.sd <- apply(motifs_train, 2, sd)
motifs_train <- motifs_train[, names(which(motifs_train.sd != 0)), drop = F]
motifs_val <- motifs_val[,colnames(motifs_train), drop = F]

# Prepare training and validation DMatrix objects
dtrain <- xgb.DMatrix(label = as.numeric(motifs_train$binary_celltype)
                      , data = as.matrix(motifs_train[, colnames(motifs_train)[colnames(motifs_train)!="binary_celltype"]]))
dvalid <- xgb.DMatrix(label = as.numeric(motifs_val$binary_celltype)
                      , data = as.matrix(motifs_val[, colnames(motifs_val)[colnames(motifs_val)!="binary_celltype"]]))

params$data <- dtrain
params$watchlist <- list(train = dtrain, validation = dvalid) 


cat("Training model...\n")
set.seed(123) 
model <- xgboost::xgb.train(
  data = params$data,
  nrounds = params$nrounds,
  watchlist = params$watchlist,
  objective = params$objective,
  eta = params$eta,
  max_depth = params$max_depth,
  subsample = params$subsample,
  colsample_bytree = params$colsample_bytree,
  nthread = params$nthread,
  eval_metric = params$eval_metric,
  params = params$params,
  feval = params$feval,
  verbose = params$verbose,
  print_every_n = params$print_every_n,
  early_stopping_rounds = params$early_stopping_rounds,
  maximize = params$maximize,
  save_period = params$save_period,
  save_name = params$save_name,
  xgb_model = params$xgb_model,
  callbacks = params$callbacks
)

cat("Saving model...\n")

#xgb.save(model, params$save_name)
saveRDS(model, params$save_name)
cat("Done\n")
