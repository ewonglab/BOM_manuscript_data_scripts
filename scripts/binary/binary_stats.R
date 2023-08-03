params <- list(input_file = NULL)

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
  cat("Usage: binary_stats.R [parameters]\n")
  cat("\n")
  cat("Options:\n")
  cat("--input_file:<file> file with predicted probabilities produced using a binary model\n")
  cat("\n")
  q("no", status = 0)
}

cat("Loading packages...\n")
required_packages <- c("yardstick", "cvAUC", "pROC")

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

#print(params[["input_file"]])


pred_tab <- read.table(file = params[["input_file"]], header =T, stringsAsFactors = F)

TP <- nrow(pred_tab[pred_tab$true_class == 1 & pred_tab$predicted_class == 1, ])

TN <- nrow(pred_tab[pred_tab$true_class == 0 & pred_tab$predicted_class == 0, ])

FP <- nrow(pred_tab[pred_tab$true_class == 0 & pred_tab$predicted_class == 1, ])

FN <- nrow(pred_tab[pred_tab$true_class == 1 & pred_tab$predicted_class == 0, ])

recall <- TP/(TP + FN)

precision <- TP/(TP + FP)

f1 <- 2*((precision*recall)/(precision+recall))

pred_tab$true_class <- factor(pred_tab$true_class, levels=c(0,1))

pr_auc_val <- pr_auc(pred_tab, truth=true_class, prob, event_level="second", estimator="binary")


## Print binary prediction statistics

cat(paste("auROC:", round(AUC(pred_tab$prob, pred_tab$true_class),4), '\n'))

cat(paste("auPR:", round(pr_auc_val$.estimate,4), '\n'))

cat(paste("Accuracy:", round(nrow(pred_tab[pred_tab$predicted_class==pred_tab$true_class,])/nrow(pred_tab),4), '\n'))

cat(paste("F1 score:", round(f1,4), '\n'))

cat(paste("Recall:", round(recall, 4), '\n'))

cat(paste("Precision:", round(precision,4), '\n'))
