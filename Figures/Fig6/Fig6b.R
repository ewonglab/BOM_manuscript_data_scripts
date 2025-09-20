library(ggplot2)
# hepg2 replicate 1
hep_rep1 <- data.frame(SRE = c("HSL1", "HSL2", "HSL3", "HSL5", "HSL6", "GSL1", "GSL2", "GSL3", "GSL4", "GSL6")
                       , FC = c(106.73, 115.02, 181.65, 27.62, 335.51, 8.71, 62.13, 8.92, 14.64, 15.82), celline = rep("HepG2", 10), replicate = rep(1, 10))
#Gm replicate 1
gm_rep1 <- data.frame(SRE = c("HSL1", "HSL2", "HSL3", "HSL5", "HSL6", "GSL1", "GSL2", "GSL3", "GSL4", "GSL6")
                      , FC = c(0.96, 0.90, 1.43, 0.70, 1.07, 6.26, 7.46, 1.98, 2.30, 3.24), celline = rep("Gm12878", 10), replicate=rep(1, 10))
# hepG2 replicate 2
hep_rep2 <- data.frame(SRE = c("HSL1", "HSL2", "HSL3", "HSL5", "HSL6", "GSL1", "GSL2", "GSL3", "GSL4", "GSL6")
                       , FC = c(18.78963480916450, 11.66315226323090, 7.10328064676170, 8.71315415852130
                                , 24.54302743732440, 0.82529248693936, 10.40343501178290, 1.41268534393619, 1.16477431131699, 1.06032828408605)
                       , celline = rep("HepG2", 10), replicate=rep(2, 10))
hep_rep2$FC <- round(hep_rep2$FC, 2)
# Gm replicate 2
gm_rep2 <- data.frame(SRE = c("HSL1", "HSL2", "HSL3", "HSL5", "HSL6", "GSL1", "GSL2", "GSL3", "GSL4", "GSL6")
                      , FC = c(0.46, 0.26, 0.69, 0.22, 0.58, 3.12, 2.81, 1.33, 1.85, 0.99), celline = rep("Gm12878", 10), replicate=rep(2, 10))

FC <- rbind(hep_rep1, gm_rep1, hep_rep2, gm_rep2)
FC$target <- gsub('[0-9]+', '', FC$SRE)

FC$log2FC <- log2(FC$FC)

p <- ggplot(FC, aes(x=target, y=log2FC, fill=target)) +
  geom_boxplot() + facet_wrap(.~ replicate+celline) + theme_classic()

pdf("luciferase_boxplots0.pdf")
p
dev.off()

p <- ggplot(FC, aes(x=celline, y=log2FC, fill=celline)) +
  geom_boxplot() + facet_wrap(.~ replicate+target) + theme_classic()

pdf("luciferase_boxplots.pdf")
p
dev.off()

p <- ggplot(FC, aes(x=celline, y=log2FC, fill=celline)) +
  geom_boxplot() + facet_wrap(.~ target) + theme_classic()

pdf("luciferase_boxplots1.pdf")
p
dev.off()



colnames(sd_values)[3] <- "FC_sd"
colnames(mean_values)[3] <- "FC_mean"

sd_values$log2FC_sd <- log2(sd_values$FC_sd)
mean_values$log2FC_mean <- log2(mean_values$FC_mean)

df <- merge(sd_values, mean_values, by = c("SRE", "celline"))
df$celline <- as.factor(df$celline)

p <- ggplot(df, aes(x = SRE, y = log2FC_mean, fill = celline)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(position=position_dodge(width=0.9), aes(ymin = log2FC_mean - log2FC_sd
                                                        , ymax = log2FC_mean + log2FC_sd), width = 0.25) +
  xlab("SRE") +
  ylab("FC") +
  # ggtitle("Barplot with Error Bars") +
  theme_classic()

pdf("luciferase_replicates1.pdf")
p
dev.off()
