# ROC curves of binary methods for Fig.1 - DNABERT, Enformer, gkm, BOM - 500bp

suppressMessages({
  library("yardstick") 
  library("cvAUC")
  library(pROC)
  library(ggplot2)
})

dnabert_path <- "/dnabert/e8.25_500bp/"
bom_path <- "/xgb/e8.25/bin_500bp/"
gkmsvm_path <- "/gkmsvm/e8.25"
enformer_path <- "/enformer/"

args <- commandArgs(trailingOnly=TRUE)

bom_file <- args[1]
bert_out <- args[2]
celltype <- sub("_500bp.*", "", bom_file)

test_1 <- read.table(file = paste0(bom_path, paste0(celltype, "_positive_500bp_test")), header = F, stringsAsFactors = F)

bom <- read.table(file = paste0(bom_path, bom_file), header = T, stringsAsFactors = F)

dnabert <- read.table(file = paste0(dnabert_path, bert_out), header = F, stringsAsFactors = F)
gkmsvm <- read.table(file = paste0(gkmsvm_path, paste0(celltype, "_500bp_t4_l11_k7_d3_pred.txt"))
                     , header = F, stringsAsFactors = F)

enformer <- read.csv(file = paste0(enformer_path, paste0(celltype, "/AUROC.csv"))
                     , header = T, stringsAsFactors = F)
colnames(enformer) <- c("raw", "actual")

gkmsvm$actual <- ifelse(gkmsvm$V1 %in% test_1$V1, 1, 0)

dnabert_labels <- read.table(file = paste0(dnabert_path, paste0(celltype, "_k6/dev.tsv"))
                             , header = T, stringsAsFactors = F, sep='\t')

dnabert$actual <- dnabert_labels$label

colnames(dnabert)[1] <- "raw"
colnames(gkmsvm)[2] <- "raw"

pred_li <- list(bom = bom, dnabert = dnabert, gkmsvm = gkmsvm, Enformer = enformer)

roc_li <- lapply(pred_li, function(x) roc(x$actual, x$raw, direction="<"))

method_colors <- c("#B89BC9", "#8B0000", "#AAD9BC", "#FF4500")

my_theme <-  theme(panel.border = element_blank(), panel.grid.major = element_blank()
                   , panel.grid.minor = element_blank()#legend.position="none"
                   , axis.line = element_line(colour = "black", size = 1)
                   , legend.title = element_blank()
                   , legend.key=element_blank()
                   , legend.position = c(0.75, 0.25)
                   , legend.key.width = unit(1.5,"cm")
                   , panel.background = element_blank()
                   , text = element_text(size=23)
                   , axis.text.x=element_text(colour="black")
                   , axis.text.y=element_text(colour="black")
                   , legend.text=element_text(size=10)
                   , axis.ticks = element_line(colour = "black"))

my_roc <- function(x){
  p <- ggroc(x, size = 1) + my_theme +
    guides(linetype = guide_legend(override.aes = list(size = 3))) +
    geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.8, color = "grey") +
    coord_equal() +
    scale_colour_manual(values=method_colors, aesthetics = c("colour", "fill")
                        , labels = c(names(x))) +
    labs(y= "Sensitivity", x = "Specificity")
  return(p)
}
                 
pdf(paste0(celltype, "_bin_methods_500bp_rocs_revised.pdf"))
my_roc(roc_li)
dev.off()

