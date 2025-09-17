library(ggplot2)
library(ggpubr)

# Initial dataset (not overlapping MPRA test set)
enh_test_dev <- read.table(file = "Dev_log2_enrichment_pred.1.txt"
                           , header = T, stringsAsFactors = F)#1258    2
enh_test_hkp <- read.table(file = "Hk_log2_enrichment_pred.1.txt"
                           , header = T, stringsAsFactors = F)#1258    2

p_dev2 <- ggplot(enh_test_dev, aes(actual, predicted)) +
  geom_density2d(color="#B27200") + geom_point(color="#0072B2")+
  theme_classic() + ggtitle("developmental")

p_hk2 <- ggplot(enh_test_hkp, aes(actual, predicted)) +
  geom_density2d(color="#B27200") + geom_point(color="#0072B2")+
  theme_classic() + ggtitle("housekeeping")

pdf("Fly_enh_activities_notOv_MPRAtest.pdf")
ggarrange(plotlist = list(p_dev2, p_hk2))
dev.off()
