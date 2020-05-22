#'##############################################################################
#' Generate eQTM catalogue figures and tables
#' This script contains the code for
#'  - Volcano plot
#'  - QC checks
#'  - Get number significants pairs and make subgroups
#'  - Make distance plots
#'  - Compare Illumina annotation vs eQTMs
#'  - Compare Non-cell adjusted  (reference) vs cell adjusted models
#'  - Compare with other eQTMs studies
#'##############################################################################

## Load libraries ####
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(S4Vectors)
library(MultiDataSet)
library(ggforce)
#library(topGO)
library(FlowSorted.Blood.450k)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(limma)

## Load datasets ####
load("results/preprocessFiles/allOverlaps.Rdata")
load("results/preprocessFiles/methyAnnotation.Rdata")
load("results/preprocessFiles/gexpAnnotation.Rdata")

## Change name when loading results
load("results/MethComBatExpResidualsNoCellAdj/allres_simP_cpgs.Rdata")
modU <- df
featsU <- featStatsDF

load("results/MethComBatExpResidualsCellAdj/allres_simP_cpgs.Rdata")
modC <- df
featsC <- featStatsDF

## Create useful vars
### data.frame only with significant pairs
sigDf <- modU %>%
  as_tibble() %>%
  filter(sigPair)

codingTCs <- subset(expAnnot, Coding == "coding")$transcript_cluster_id
sigTCs <- unique(sigDf$TC)

### data.frame with results and overlaps info
modU_comp <- as_tibble(inner_join(modU, overDF, by = c("CpG", "TC")))

### Modify methylation annotation
methyAnnot <- methyAnnot %>%
  as_tibble() %>%
  mutate(CpG = Name)

# General Overview ####
## Statistics ####
### Significant pairs
nrow(sigDf)
# [1] 63831

### Significant TC
length(unique(sigDf$TC))
# [1] 11071

### Significant coding TC
sum(unique(sigDf$TC) %in% codingTCs)
# [1] 7874

### Significant CpGs
length(unique(sigDf$CpG))
# [1] 35228

### Proportion negative pairs
sum(sigDf$FC < 0)
# [1] 38310

mean(sigDf$FC < 0)*100
# [1] 60.01786

## Distribution CpGs/TC
CpG_plot <- sigDf %>%
  group_by(CpG) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = n)) + geom_histogram(binwidth = 1) +
  scale_x_continuous("", limits = c(0, 50)) +
  scale_y_continuous("")  +
  ggtitle("TCs associated with each CpG")

TC_plot <- sigDf %>%
  group_by(TC) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = n)) + geom_histogram(binwidth = 1) +
  scale_x_continuous("", limits = c(0, 50)) +
  scale_y_continuous("")  +
  ggtitle("CpGs associated with each TC")

png("paper/eQTMs_CpGs_TC_distr.png", width = 2000, height = 1500, res = 300)
plot_grid(CpG_plot, TC_plot, labels = "AUTO", nrow = 2)
dev.off()

sigDf %>%
  group_by(TC) %>%
  summarize(n = n()) %>%
  `$`("n") %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.000   2.000   5.766   6.000 176.000

sigDf %>%
  group_by(CpG) %>%
  summarize(n = n()) %>%
  `$`("n") %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.000   1.000   1.812   2.000  26.000



## Volcano plot ####
vol_all <- ggplot(modU, aes(x = FC/10, y = -log10(p.value))) +
  geom_point(alpha = 0.01) + 
  scale_x_continuous(name = "-log2 FC per 10% methylation change") +
  ylab(expression(-log[10](p.value))) +
  theme_bw()

vol_filt <- modU %>%
  filter(abs(FC) < 5) %>%
  ggplot(aes(x = FC/10, y = -log10(p.value))) +
  geom_point(alpha = 0.01) + 
  scale_x_continuous(name = "-log2 FC per 10% methylation change") +
  ylab(expression(-log[10](p.value))) +
  theme_bw()

png("paper/volcano_all.png", width = 2000, height = 2000, res = 300)
volcano_plot(modU$p.value, modU$FC/100, paste(modU$CpG, modU$TC), 
             tPV = -log10(1e-8), tFC = 0.01, show.labels = FALSE) +
  geom_point(alpha = 0.1)
dev.off()


## QCs ####
### CpGs variability
featsU_var <- methyAnnot %>%
  as_tibble() %>%
  mutate(feat = Row.names) %>%
  dplyr::select(feat, meth_range, variability) %>%
  right_join(featsU, by = "feat") %>%
  mutate(sig = ifelse(p.val.adj < 0.05, "significant", "random"))
table(featsU_var$variability, featsU_var$sig)
#           random significant
# invariant 169446        4647
# variant   181744       30581
chisq.test(table(featsU_var$variability, featsU_var$sig))
# data:  table(featsU_var$variability, featsU_var$sig)
# X-squared = 15894, df = 1, p-value < 2.2e-16

t <- table(featsU_var$variability, featsU_var$sig)
t[1]/t[2]/t[3]*t[4]

png("paper/CpGVar_eQTMs.png", width = 2000, height = 1000, res = 300)
featsU_var %>%
  mutate(sigVar = ifelse(sig == "random", "CpG not in eQTMs", "CpG in eQTMs")) %>%
  ggplot(aes(x = sigVar, y = meth_range, fill = sigVar)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(name = "") +
  theme(legend.position = "none") +
  scale_y_continuous(name = "CpG variability") +
  scale_fill_manual(values = c("#999999", "#FFFFFF"))
dev.off()



### TC call rate
int_TC <- expAnnot %>%
  dplyr::select(probeset_id, Expressed, CallRate) %>%
  mutate(sig = ifelse(probeset_id %in% sigTCs, "TCs in eQTMs", "TCs without eQTMs"))
table(int_TC$CallRate < 90, int_TC$sig)
# TCs in eQTMs TCs without eQTM
# FALSE       8584        19317
# TRUE        2487        30304

chisq.test(table(int_TC$CallRate > 90, int_TC$sig))
# X-squared = 1149, df = 1, p-value < 2.2e-16
t <- table(int_TC$CallRate < 90, int_TC$sig)
t[1]/t[2]/t[3]*t[4]

png("paper/CallRate_eQTMs.png", width = 2000, height = 1000, res = 300)
int_TC %>%
  ggplot(aes(x = sig, y = CallRate, fill = sig)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(name = "") +
  theme(legend.position = "none") +
  scale_y_continuous(name = "TC call rate") +
  scale_fill_manual(values = c("#999999", "#FFFFFF"))
dev.off()

summary(glm(CallRate/100 ~ sig, int_TC, family = "binomial"))

### Classify CpG in groups ####
CpGsSum <- modU %>%
  group_by(CpG) %>%
  summarise(Type = ifelse(sum(sigPair) == 0, "Non-significant",
                          ifelse(sum(sigPair) == 1, "Mono", "Multi")),
            Direction = ifelse(sum(sigPair) == 0, "Non-significant",
                               ifelse(all(FC[sigPair] > 0), "Positive", 
                                      ifelse(all(FC[sigPair] < 0), "Inverse", "Both")))) %>%
  mutate(Combined = ifelse(Type == "Non-significant", 
                           "Non-significant", 
                           paste(Type, Direction, sep = "_")))
t <- table(CpGsSum$Type, CpGsSum$Direction)
t <- t[c("Mono", "Multi"), c("Inverse", "Positive", "Both")]
t
#       Inverse Positive  Both
# Mono    12637     8961     0
# Multi    6135     3530  3965
chisq.test(t[, -3])
# X-squared = 68.441, df = 1, p-value < 2.2e-16
t[2]/t[1]/t[4]*t[3]
# [1] 0.8114236
addmargins(prop.table(t))
# Inverse  Positive      Both       Sum
# Mono  0.3587203 0.2543715 0.0000000 0.6130919
# Multi 0.1741512 0.1002044 0.1125525 0.3869081
# Sum   0.5328716 0.3545759 0.1125525 1.0000000




sink("paper/CpGs_type.txt")
addmargins(t)
sink()

# Distance plots ####
# Distance distribution
## Signif vs no-signif
dist_all <- ggplot(modU_comp, aes(x = Distance, color = sigPair)) + geom_density() + 
  theme_bw() + 
  scale_color_discrete(name = "", labels = c("non-eQTM", "eQTM")) +
  scale_x_continuous(breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "") + 
  ggtitle("CpG-TC Distance (all pairs)") + 
  theme(plot.title = element_text(hjust = 0.5))

ks.test(filter(modU_comp, sigPair) %>% pull(., Distance),
        filter(modU_comp, !sigPair) %>% pull(., Distance))
# Two-sample Kolmogorov-Smirnov test
# 
# data:  filter(modU_comp, sigPair) %>% pull(., Distance) and filter(modU_comp, !sigPair) %>% pull(., Distance)
# D = 0.15302, p-value < 2.2e-16
# alternative hypothesis: two-sided

wilcox.test(filter(modU_comp, sigPair) %>% pull(., Distance))
# Wilcoxon rank sum test with continuity correction
# 
# data:  filter(modU_comp, sigPair) %>% pull(., Distance) and filter(modU_comp, !sigPair) %>% pull(., Distance)
# V = 932277653, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0


modU_comp %>% 
  group_by(sigPair) %>%
  summarize(m = median(Distance), 
            l = quantile(Distance, 0.25), 
            h = quantile(Distance, 0.75))

## Distance per direction
dist_dir <- modU_comp %>%
  filter(sigPair) %>%
  mutate(Direction = ifelse(FC > 0, "Positive", "Inverse")) %>%
  ggplot(aes(x = Distance, color = Direction)) + geom_density() + 
  theme_bw() + 
  scale_color_manual(name = "", values = c("#000000", "#009E73")) +
  scale_x_continuous(breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "") + 
  ggtitle("CpG-TC Distance (significant pairs)") + 
  theme(plot.title = element_text(hjust = 0.5))

modU_comp %>% 
  filter(sigPair) %>%
  mutate(Direction = ifelse(FC > 0, "Positive", "Inverse")) %>%
  group_by(Direction) %>%
  summarize(m = median(Distance), 
            l = quantile(Distance, 0.25), 
            h = quantile(Distance, 0.75))


png("paper/distance_distr.png", width = 2000, height = 2000, res = 300)
plot_grid(dist_all, dist_dir, nrow = 2, labels = c("A", "B"))
dev.off()


## Distance per direction + eQTM 
png("paper/distance_distr.png", width = 2000, height = 1300, res = 300)
modU_comp %>%
  mutate(Direction = ifelse(!sigPair, "non-eQTM", ifelse(FC > 0, "Positive", "Inverse")),
         Direction = factor(Direction, levels = c("Inverse", "Positive", "non-eQTM"))) %>%
  ggplot(aes(x = Distance, color = Direction)) + geom_density() + 
  theme_bw() + 
  scale_color_manual(name = "", values = c("#E69F00", "#009E73", "#000000")) +
  scale_x_continuous(breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "") + 
  ggtitle("CpG-TC Distance") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()



## Signif per groups
dist_groups <- modU_comp %>%
  filter(sigPair) %>%
  left_join(CpGsSum, by = "CpG") %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both"))) %>%
  ggplot(aes(y = Distance, x = Direction, fill = Direction)) + 
  geom_boxplot() + 
  theme_bw() + 
  scale_fill_manual(name = "", 
                     breaks = c("Inverse", "Positive", "Both"),
                     values = c("turquoise", "coral", "grey70")) +
  facet_grid(~ Type, scales = "free_x") +
  scale_y_continuous(breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  ggtitle("CpG-TC Distance (Significant pairs)") + 
  theme(plot.title = element_text(hjust = 0.5))

modU_comp %>%
  filter(sigPair) %>%
  left_join(CpGsSum, by = "CpG") %>%
  group_by(Direction) %>%
  summarize(m = median(Distance), 
            l = quantile(Distance, 0.25), 
            h = quantile(Distance, 0.75))
# Direction     m        l       h
# <chr>     <int>    <dbl>   <dbl>
#   1 Both      -1679 -138102  112303
# 2 Inverse    -686  -98790.  74250.
# 3 Positive  -7771 -131119   82650.

modU_comp %>%
  filter(sigPair) %>%
  left_join(CpGsSum, by = "CpG") %>%
  group_by(Type) %>%
  summarize(m = median(Distance), 
            l = quantile(Distance, 0.25), 
            h = quantile(Distance, 0.75))
# Type       m        l      h
# <chr>  <dbl>    <dbl>  <dbl>
#   1 Mono  -2142. -121440. 75250.
# 2 Multi  -967  -115541  89436

modU_comp %>%
  filter(sigPair) %>%
  left_join(CpGsSum, by = "CpG") %>%
  group_by(Combined) %>%
  summarize(m = median(Distance), 
            l = quantile(Distance, 0.25), 
            h = quantile(Distance, 0.75)) %>%
  mutate(d = h - l)
# Combined            m        l      h       d
# 1 Mono_Inverse    -1094  -95897   62932 158829
# 2 Mono_Positive  -10839 -150837   86773 237610
# 3 Multi_Both      -1679 -138102  112303 250405
# 4 Multi_Inverse    -342 -101410.  80146 181556.
# 5 Multi_Positive  -4760 -110583   78594 189177

## Test for Multi Inverse - data with median closest to 0.
modU_comp %>%
  filter(sigPair) %>%
  left_join(CpGsSum, by = "CpG") %>%
  filter(Combined == "Multi_Inverse") %>%
  pull(., Distance) %>%
  wilcox.test()
# Wilcoxon signed rank test with continuity correction
# 
# data:  .
# V = 82114591, p-value = 2.781e-10
# alternative hypothesis: true location is not equal to 0

modU_comp %>%
  filter(sigPair) %>%
  left_join(CpGsSum, by = "CpG") %>%
  lapply(c("Mono_Inverse", "Mono_Positive", "Multi_Both", "Multi_Inverse", "Multi_Positive"),
         function(x, y) wilcox.test(subset(y, Combined == x)$Distance), y = .)
# 
# [[1]]
# 
# Wilcoxon signed rank test with continuity correction
# 
# data:  subset(y, Combined == x)$Distance
# V = 35386210, p-value < 2.2e-16
# alternative hypothesis: true location is not equal to 0
# 
# 
# [[2]]
# 
# Wilcoxon signed rank test with continuity correction
# 
# data:  subset(y, Combined == x)$Distance
# V = 17397764, p-value < 2.2e-16
# alternative hypothesis: true location is not equal to 0
# 
# 
# [[3]]
# 
# Wilcoxon signed rank test with continuity correction
# 
# data:  subset(y, Combined == x)$Distance
# V = 44268445, p-value = 4.402e-12
# alternative hypothesis: true location is not equal to 0
# 
# 
# [[4]]
# 
# Wilcoxon signed rank test with continuity correction
# 
# data:  subset(y, Combined == x)$Distance
# V = 82114591, p-value = 2.781e-10
# alternative hypothesis: true location is not equal to 0
# 
# 
# [[5]]
# 
# Wilcoxon signed rank test with continuity correction
# 
# data:  subset(y, Combined == x)$Distance
# V = 21816864, p-value = 7.207e-16
# alternative hypothesis: true location is not equal to 0
# 
png("paper/distance_distribution.png", width = 2000, height = 2000, res = 300)
plot_grid(dist_all, dist_groups, nrow = 2)
dev.off()

# Effect size ####
summary(abs(subset(modU_comp, sigPair)$FC/100))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000237 0.006348 0.011786 0.021641 0.022983 1.643445

effectAll <- modU_comp %>%
  filter(sigPair) %>%
  mutate(Dir = ifelse(FC > 0, "Positive", "Inverse")) %>%
  ggplot(aes(x = Dir, y = abs(FC/10), fill = Dir)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 0.5), name = "abs(log2 FC)/10% Methylation") +
  theme_bw() +
  scale_x_discrete(name = "eQTM Effect Direction") +
  theme(legend.position = "none")

ks.test(abs(subset(sigDf, FC > 0)$FC)/10,
        abs(subset(sigDf, FC < 0)$FC)/10)

modU_comp %>%
  filter(sigPair) %>%
  mutate(Dir = ifelse(FC > 0, "Positive", "Inverse"),
         FC = FC/10) %>%
  group_by(Dir) %>%
  summarize(m = median(FC),
            l = quantile(FC, 0.25),
            h = quantile(FC, 0.75))

# Dir           m       l       h
# <chr>     <dbl>   <dbl>   <dbl>
#   1 Inverse  -0.122 -0.240  -0.0655
# 2 Positive  0.112  0.0607  0.214
# 

effect_groups <- modU_comp %>%
  filter(sigPair) %>%
  left_join(CpGsSum, by = "CpG") %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both"))) %>%
  ggplot(aes(x = Direction, y = abs(FC/10), fill = Direction)) + 
  geom_boxplot() +
  theme_bw() + 
  scale_fill_manual(name = "", 
                    breaks = c("Inverse", "Positive", "Both"),
                    values = c("turquoise", "coral", "grey70")) +
  facet_grid(~ Type, scale = "free_x") +
  scale_y_continuous(name = "abs(log2 FC)/10% Methylation", 
                     limits = c(0, 0.5)) +
  scale_x_discrete(name = "CpG Type") + 
  ggtitle("Effect size in significant pairs") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

modU_comp %>%
  filter(sigPair) %>%
  left_join(CpGsSum, by = "CpG") %>%
  mutate(FC = abs(FC)/10) %>%
  group_by(Combined) %>%
  summarize(m = median(FC),
            l = quantile(FC, 0.25),
            h = quantile(FC, 0.75))
# Combined             m       l      h
# <chr>            <dbl>   <dbl>  <dbl>
#   1 Mono_Inverse   0.113  0.0628 0.219
# 2 Mono_Positive  0.111  0.0605 0.209
# 3 Multi_Both     0.0992 0.0543 0.187
# 4 Multi_Inverse  0.136  0.0727 0.274
# 5 Multi_Positive 0.127  0.0678 0.251


png("paper/effect_distribution.png", width = 2000, height = 1500, res = 300)
plot_grid(effectAll, effect_groups, nrow = 2)
dev.off()





# Distance vs Effect size
png("paper/distance_effect.png", width = 2000, height = 1300, res = 300)
modU_comp %>%
  filter(sigPair) %>%
  mutate(Direction = ifelse(FC > 0, "Positive", "Inverse")) %>%
  ggplot(aes(x = Distance, y = FC/10, color = Direction)) + 
  geom_point(alpha = 0.1) + 
  theme_bw() + 
  scale_color_manual(name = "", 
                     breaks = c("Inverse", "Positive"),
                     values = c("#E69F00", "#009E73")) +
  scale_x_continuous(breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "log2 FC/0.1 Methylation") + 
  ggtitle("CpG-TC distance vs effect size") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


modU_comp %>%
  filter(sigPair) %>%
  mutate(Direction = ifelse(FC > 0, "Positive", "Inverse")) %>%
  cor.test(x = pull(., Distance), 
           y = pull(., FC),
           cot = .)

modU_comp %>%
  filter(sigPair) %>%
  lm(abs(FC) ~ abs(Distance), .) %>%
  summary()
# Call:
#   lm(formula = abs(FC) ~ abs(Distance), data = .)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max
# -2.139  -1.529  -0.985   0.134 162.182
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)    2.165e+00  2.114e-02 102.407   <2e-16 ***
#   abs(Distance) -7.580e-09  9.968e-08  -0.076    0.939
# 


## Illumina Annotation ####
modU_Annot <- modU %>%
  as_tibble() %>%
  dplyr::select(CpG, TC, FC, p.value, sigPair) %>%
  left_join(CpGsSum, by = "CpG") %>%
  left_join(dplyr::select(methyAnnot, CpG, UCSC_RefGene_Group, UCSC_RefGene_Name), by = "CpG") %>%
  left_join(dplyr::select(expAnnot, TC, GeneSymbol_Affy), by = "TC") %>%
  mutate(GeneAffy = strsplit(GeneSymbol_Affy, ";"))

## Create list with genes in common in Affy and Illumina Annotation
methGenes <- unique(unlist(modU_Annot$UCSC_RefGene_Name))
expGenes <- unique(unlist(modU_Annot$GeneAffy))
comGenes <- intersect(methGenes, expGenes)

## Create list with genes in common in Affy and Illumina Annotation that are
## part of an eQTM
modU_Annot_f <- modU_Annot %>%
  filter(sigPair)

methGenesf <- unique(unlist(modU_Annot_f$UCSC_RefGene_Name))
expGenesf <- unique(unlist(modU_Annot_f$GeneAffy))
comGenesf <- intersect(methGenesf, expGenesf)

modU_Annot_f %>%
  dplyr::select(CpG, UCSC_RefGene_Name) %>%
  distinct() %>%
  summarize(n = n(),
            n_g = sum(UCSC_RefGene_Name != ""),
            n_GE = sum(sapply(UCSC_RefGene_Name, function(x) any(x %in% comGenes))))


modU_Annot_f2 <- modU_Annot_f %>%
  filter(sapply(UCSC_RefGene_Name, function(x) any(x %in% comGenes))) %>%
  mutate(match = sapply(seq_len(n()), function(x) any(GeneAffy[[x]] %in% UCSC_RefGene_Name[[x]])))


table(modU_Annot_f2$match)
modU_Annot_f2 %>%
  mutate(Dir = ifelse(FC > 0, "Positive", "Inverse")) %>%
  group_by(Dir) %>%
  summarize(N = sum(match),
            P = mean(match))



annotated <- modU_Annot_f2 %>%
  group_by(CpG) %>%
  summarize(found = any(match)) %>%
  left_join(CpGsSum, by = "CpG")
mean(annotated$found)
table(annotated$found)


annotated %>%
  group_by(Combined) %>%
  summarize(N = sum(found),
            P = mean(found))


# Enrichment by gene position based on annotation ####
modU_Annot <- modU_Annot %>%
  mutate(gene_match = sapply(seq_len(n()), function(x) any(GeneAffy[[x]] %in% UCSC_RefGene_Name[[x]])))

t <- table(modU_Annot$sigPair, modU_Annot$gene_match)
t[1]/t[2]/t[3]*t[4]
chisq.test(t)

## Expand CpGs and TCs to make all pairs between CpG Gene and TC
## Restrict pairs where CpG gene and TC gene match
modU_gene_pairs <- modU_Annot %>%
  filter(gene_match) %>%
  mutate(Gene_Group = strsplit(UCSC_RefGene_Group, ";")) %>%
  unnest(Gene_Group, UCSC_RefGene_Name) %>%
  mutate(GeneAffy = strsplit(GeneSymbol_Affy, ";")) %>%
  unnest(GeneAffy) %>%
  dplyr::select(CpG, TC, FC, sigPair, p.value, UCSC_RefGene_Name, Gene_Group, GeneAffy) %>%
  filter(UCSC_RefGene_Name == GeneAffy) %>%
  distinct()
  

gpos <- c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR")
names(gpos) <- gpos
#gpos <- c(as.list(gpos), list(Comb = c("TSS1500", "5'UTR", "Body")))

ORs <- lapply(gpos, function(x) {
  t <- table(modU_gene_pairs$Gene_Group %in% x, modU_gene_pairs$sigPair)
  list(p.value = chisq.test(t)$p.value, OR = t[1]/t[2]/t[3]*t[4])
  }) %>%
  Reduce(f = cbind)
colnames(ORs) <- gpos  


ORs2 <- lapply(gpos, function(x) {
  cot <- filter(modU_gene_pairs, sigPair)
  t <- table(cot$Gene_Group %in% x, sign(cot$FC))
  list(p.value = chisq.test(t)$p.value, OR = t[1]/t[2]/t[3]*t[4])
}) %>%
  Reduce(f = cbind)
colnames(ORs2) <- gpos  


# Compare models ####
sigDf2 <- modC %>%
  as_tibble() %>%
  filter(sigPair)
## Statistics ####
### Significant pairs
nrow(sigDf2)
# [1] 39749

table(sign(sigDf2$FC))

### Significant TC
length(unique(sigDf2$TC))
# [1] 8886

### Significant coding TC
sum(unique(sigDf2$TC) %in% codingTCs)
# [1] 6288

### Significant CpGs
length(unique(sigDf2$CpG))
# [1] 21966

sigDf2 %>%
  group_by(TC) %>%
  summarize(n = n()) %>%
  `$`("n") %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.000   2.000   4.473   5.000 128.000
sigDf2 %>%
  group_by(CpG) %>%
  summarize(n = n()) %>%
  `$`("n") %>%
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.00    1.00    1.00    1.81    2.00   24.00


### Classify CpG in groups ####
CpGsSum2 <- modC %>%
  group_by(CpG) %>%
  summarise(Type = ifelse(sum(sigPair) == 0, "Non-significant",
                          ifelse(sum(sigPair) == 1, "Mono", "Multi")),
            Direction = ifelse(sum(sigPair) == 0, "Non-significant",
                               ifelse(all(FC[sigPair] > 0), "Positive", 
                                      ifelse(all(FC[sigPair] < 0), "Inverse", "Both")))) %>%
  mutate(Combined = ifelse(Type == "Non-significant", 
                           "Non-significant", 
                           paste(Type, Direction, sep = "_")))
t <- table(CpGsSum2$Type, CpGsSum2$Direction)
t <- t[c("Mono", "Multi"), c("Inverse", "Positive", "Both")]
t


#       Inverse Positive Both
# Mono     8084     5681    0
# Multi    3738     2400 2063
chisq.test(t[, -3])
# X-squared = 8.2039, df = 1, p-value = 0.00418

sink("paper/CpGs_type2.txt")
addmargins(t)
sink()



## Data.frame with groups comparison
### General
adjList <- list(pairs = paste(sigDf$CpG, sigDf$TC, sep = "-"),
                 TCs = unique(sigDf$TC),
                 CpG = unique(sigDf$CpG))
adjList$codTC <- adjList$TCs[adjList$TCs %in% codingTCs]
cellList <- list(pairs = paste(sigDf2$CpG, sigDf2$TC, sep = "-"),
                 TCs = unique(sigDf2$TC),
                 CpG = unique(sigDf2$CpG))
cellList$codTC <- cellList$TCs[cellList$TCs %in% codingTCs]


modCompDF <- data.frame(Adj = lengths(adjList),
                    Cell = lengths(cellList),
                    Shared = unlist(Map(function(x, y) length(intersect(x, y)), adjList, cellList)),
                    AdjSp = unlist(Map(function(x, y) length(setdiff(x, y)), adjList, cellList)),
                    CellSp = unlist(Map(function(x, y) length(setdiff(x, y)), cellList, adjList))) %>%
  mutate(AdjSpProp = round(AdjSp/Adj*100, 1),
         CellSpProp = round(CellSp/Cell*100, 1))
write.table(modCompDF[c(1, 2, 4, 3), c(1:4, 6, 5, 7)], sep = "\t",
            file = "paper/ModelsStatisticComparison.txt",
            row.names = FALSE, quote = FALSE)


### Per CpG type
adjCpgs <- CpGsSum %>% 
  group_by(Combined) %>%
  group_split()
names(adjCpgs) <- sapply(adjCpgs, function(x) x$Combined[1])
adjCpgs <- sapply(adjCpgs, pull, CpG)

cellCpgs <- CpGsSum2 %>% 
  group_by(Combined) %>%
  group_split()
names(cellCpgs) <- sapply(cellCpgs, function(x) x$Combined[1])
cellCpgs <- sapply(cellCpgs, pull, CpG)

modCompDF2 <- data.frame(Adj = lengths(adjCpgs),
                        Cell = lengths(cellCpgs),
                        Shared = unlist(Map(function(x, y) length(intersect(x, y)), adjCpgs, cellCpgs)),
                        AdjSp = unlist(Map(function(x, y) length(setdiff(x, y)), adjCpgs, cellCpgs)),
                        CellSp = unlist(Map(function(x, y) length(setdiff(x, y)), cellCpgs, adjCpgs))) %>%
  mutate(AdjSpProp = round(AdjSp/Adj*100, 1),
         CellSpProp = round(CellSp/Cell*100, 1))

write.table(modCompDF2[1:5, c(1:4, 6, 5, 7)], sep = "\t",
            file = "paper/ModelsStatisticComparisonCpG.txt",
            row.names = FALSE, quote = FALSE)

parPlot <- CpGsSum %>%
  left_join(select(CpGsSum2, CpG, Combined), by = "CpG") %>%
  group_by(Combined.x, Combined.y) %>%
  summarize(n = n()) %>%
  filter(!(Combined.x == "Non-significant" & Combined.y == "Non-significant")) %>%
  ungroup() %>%
  mutate(Combined.x = factor(Combined.x, 
                             levels = c("Mono_Inverse", "Mono_Positive", "Multi_Inverse", 
                                        "Multi_Positive", "Multi_Both", "Non-significant")),
         Combined.y = factor(Combined.y, 
                             levels = c("Mono_Inverse", "Mono_Positive", "Multi_Inverse", 
                                        "Multi_Positive", "Multi_Both", "Non-significant"))) %>%
  gather_set_data(1:2) %>%
  ggplot(aes(x, id = id, split = y, value = n)) +
    geom_parallel_sets(aes(fill = Combined.y), alpha = 0.6, axis.width = 0.4) +
    geom_parallel_sets_axes(axis.width = 0.4, fill = "gray95",
                            color = "gray80") +
    geom_parallel_sets_labels(colour = 'gray35', size = 4.5, angle = 0, fontface = "bold") +
  scale_x_discrete(labels = c("Adjusted", "Cell Adjusted")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.title.x  = element_blank()
  )

png("paper/CompModelsCpGs.png", width = 2500, height = 2500, res = 300)
parPlot
dev.off()


### Merge dataset ####
mergeTB <- modU %>%
  left_join(modC, by = c("CpG", "TC")) %>%
  as_tibble() %>%
  filter(sigPair.x == TRUE | sigPair.y == TRUE) %>%
  mutate(sigType = ifelse(sigPair.x == TRUE, ifelse(sigPair.y == TRUE, "Shared", "Main"), "Cell"))


t <- mergeTB %>%
  group_by(TC) %>%
  summarize(main = sum(sigPair.x),
            cell = sum(sigPair.y)) %>%
  gather(model, nCpG, 2:3) %>%
  group_by(model) %>%
  summarize(mono = sum(nCpG == 1),
            multi = sum(nCpG > 1))
prop.table(data.matrix(t[,-1]), margin = 1)


## Distribution CpGs/TC
CpG_plot2 <- mergeTB %>%
  group_by(CpG) %>%
  summarize(main = sum(sigPair.x),
            cell = sum(sigPair.y)) %>%
  gather(model, nTC, 2:3) %>%
  filter(nTC > 0) %>%
  mutate(model = factor(model, levels = c("main", "cell")),
         TCs = ifelse(nTC > 10, "11+", 
                      ifelse(nTC > 5, "6-10", nTC)),
         TCs = factor(TCs, levels = c(1:5, "6-10", "11+"))) %>%
  group_by(model, TCs) %>%
  summarize(n = n()) %>%
  group_by(model) %>%
  mutate(p = n/sum(n)*100) %>%
  ggplot(aes(x = TCs, y = p, fill = model)) +  
  geom_bar(position = "dodge", stat = "identity") +
  scale_y_continuous("Percentage CpGs")  +
  scale_fill_manual(name = "Model", labels = c("Main", "Cell"), 
                    values = c("#F5793A", "#0F2080")) +
  ggtitle("TCs associated with each CpG") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")


TC_plot2 <- mergeTB %>%
  group_by(TC) %>%
  summarize(Main = sum(sigPair.x),
            Cell = sum(sigPair.y)) %>%
  gather(model, nCpG, 2:3) %>%
  filter(nCpG > 0) %>%
  mutate(model = factor(model, levels = c("Main", "Cell")),
         CpGs = ifelse(nCpG > 20, "21+", 
                       ifelse(nCpG > 10, "11-20", 
                              ifelse(nCpG > 5, "6-10", nCpG))),
         CpGs = factor(CpGs, levels = c(1:5, "6-10", "11-20", "21+"))) %>%
  group_by(model, CpGs) %>%
  summarize(n = n()) %>%
  group_by(model) %>%
  mutate(p = n/sum(n)*100) %>%
  ggplot(aes(x = CpGs, y = p, fill = model)) +  
  geom_bar(position = "dodge", stat = "identity") +
  scale_y_continuous("Percentage TCs")  +
  scale_fill_manual(name = "Model",
                    values = c("#F5793A", "#0F2080")) +
  ggtitle("CpGs associated with each TC") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
        
png("paper/eQTMs_CpGs_TC_distr_cell.png", width = 2500, height = 1500, res = 300)
plot_grid(CpG_plot2, TC_plot2, labels = "AUTO", ncol = 2)
dev.off()

pvals_comp <- mergeTB %>%
  mutate(sigType = factor(sigType, levels = c("Shared", "Main", "Cell"))) %>%
  ggplot(aes(x = -log10(p.value.y), y = -log10(p.value.x), col = sigType)) +
  geom_point() +
  scale_x_continuous(name = "Cell adjusted") + 
  scale_y_continuous("Main model") + 
  ggtitle("-log10 p-values comparative") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  geom_abline(slope = 1, linetype = "dashed") +
  facet_wrap(. ~ sigType) +
  scale_color_manual(name = "", values = c("#85C0F9", "#F5793A", "#0F2080"))
  
png("paper/CompModelsP_values.png", width = 1500, height = 1200, res = 300)
pvals_comp
dev.off()

## Distance distribution
png("paper/dist_distr_main_eQTMs.png", width = 1700, height = 800, res = 300)
mergeTB %>%
  filter(sigPair.x) %>%
  inner_join(overDF, by = c("TC", "CpG")) %>%
  mutate(sigType = factor(sigType, levels = c("Shared", "Main"))) %>%
  ggplot(aes(x = Distance, color = sigType)) +
  geom_density() +
  theme_bw() +
  scale_color_manual(name = "", values = c("#85C0F9", "#F5793A")) +
  scale_x_continuous(breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "") + 
  ggtitle("CpG-TC Distance") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("paper/dist_distr_cell.png", width = 2500, height = 1500, res = 300)
rbind(mutate(sigDf, mod = "Main"), 
                    mutate(sigDf2, mod = "Cell")) %>%
  inner_join(overDF, by = c("TC", "CpG")) %>%
  ggplot(aes(x = Distance, color = mod)) +
  geom_density() +
  theme_bw() +
  scale_color_manual(name = "Model", labels = c("Main", "Cell"), 
                    values = c("#E69F00", "#0072B2")) +
  scale_x_continuous(breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "") + 
  ggtitle("CpG-TC Distance") + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

mergeTB %>%
  filter(sigPair.x) %>%
  inner_join(overDF, by = c("TC", "CpG")) %>%
  mutate(mod = ifelse(sigType == "Both", "Common", sigType)) %>%
  group_by(mod) %>%
  summarize(m = median(Distance), 
            l = quantile(Distance, 0.25), 
            h = quantile(Distance, 0.75))

# mergeTB <- mergeTB %>%
#   mutate(Diff.pval = pnorm(abs(FC.x - FC.y)/sqrt(SD.x**2 + SD.y**2), lower.tail = FALSE), 
#          isDiff = ifelse(Diff.pval < 0.05, "Different", "Equal"))
## Compare estimates
top <- filter(mergeTB, sigType == "Common") %>%
  ggplot(aes(x = FC.y/10, y = FC.x/10)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_point(color = "#0072B2") +
  scale_x_continuous(name = "Cell adjusted") + 
  scale_y_continuous("Main model") + 
  ggtitle("Common") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(name = "") +
  geom_smooth(method = "lm")

bottom <- filter(mergeTB, sigType != "Common") %>%
  ggplot(aes(x = FC.y/10, y = FC.x/10, col = sigType)) +
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed") +
  scale_x_continuous(name = "Cell adjusted") + 
  scale_y_continuous("Main model") + 
 # ggtitle("Estimates comparative") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_color_manual(name = "", values = c("#999999", "#E69F00")) +
  facet_wrap(. ~ sigType) +
  geom_smooth(method = "lm") 


estim_comp <- mergeTB %>% 
  mutate(sigType = factor(sigType, levels = c("Shared", "Main", "Cell"))) %>%
  ggplot(aes(x = FC.y/10, y = FC.x/10, col = sigType)) +
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed") +
  scale_x_continuous(name = "Cell adjusted") + 
  scale_y_continuous("Main model") + 
  ggtitle("Estimates comparative") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_color_manual(name = "", values = c("#85C0F9", "#F5793A", "#0F2080")) +
  facet_wrap(. ~ sigType) +
  geom_smooth(method = "lm") 


png("paper/CompModelsEstimatesPvals.png", width = 2500, height = 2000, res = 300)
plot_grid(estim_comp, pvals_comp, nrow = 2, labels = c("A", "B"))
dev.off()

png("paper/CompModelsEstimates.png", width = 2500, height = 2500, res = 300)
plot_grid(top, bottom, nrow = 2)
dev.off()

filter(mergeTB, sigType == "Main") %>% lm(FC.y ~ FC.x, .) %>% summary()
filter(mergeTB, sigType == "Cell") %>% lm(FC.y ~ FC.x, .) %>% summary()
filter(mergeTB, sigType == "Both") %>% lm(FC.y ~ FC.x, .) %>% summary()


## Run enrichment on CpGs ####
mainCpGs <- setdiff(adjList$CpG, cellList$CpG)
comCpGs <- intersect(adjList$CpG, cellList$CpG)

methyAnnot %>% 
  filter(CpG %in% c(mainCpGs, comCpGs)) %>%
  mutate(GeneRel = ifelse(UCSC_RefGene_Name == "", "Intergenic", "Genic"),
         Type = ifelse(CpG %in% mainCpGs, "Main", "Common")) %>%
  dplyr::select(GeneRel, Type) %>%
  table()

chromStates <- c("TssA", "TssAFlnk", "TxFlnk", "TxWk", "Tx", "EnhG", "Enh",
                 "ZNF.Rpts", "Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC",
                 "ReprPCWk", "Quies")

getOR <- function(cols, df){
  t <- data.matrix(df[, cols])
  or <- t[1, 1]/(t[1, 2])/(t[2, 1])*t[2, 2]
  p.val <- chisq.test(t)$p.value
  ORl <- log(or)
  SEl <- sqrt(1/t[1, 1] + 1/t[1, 2] + 1/t[2, 1] + 1/t[2, 2])
  c(OR = or, p.val = p.val, ORm = exp(ORl - 1.96*SEl), 
    ORM = exp(ORl + 1.96*SEl))
}

g <- function(x, gpos, cols = c("in", "ou")){
  sapply(gpos, function(y) getOR(paste0(y, cols), df = x))
}

sum2 <- function(x) sum(!x)

methChromSt <- as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, eval(chromStates)) %>%
  mutate(Type = ifelse(CpG %in% mainCpGs, "Main", ifelse(CpG %in% comCpGs, "Shared", "non-eQTMs"))) %>%
  group_by(Type) %>%
  summarize_at(chromStates, list(sum = sum, sum2 = sum2, prop = mean))


chromSt_modComp_prop_plot <- methChromSt %>%
  select(Type, ends_with("prop")) %>%
  gather(categories, proportion, 2:(2+length(chromStates) - 1)) %>%
  ungroup() %>%
  mutate(Type = factor(Type,  levels = c("non-eQTMs", "Shared", "Main")),
         categories = gsub("_prop", "", categories),
         Group = factor(ifelse(categories %in% c("TssA", "TssAFlnk"), "TssProxProm",
                               ifelse(categories %in% c("Tx", "TxWk"), "ActTrans", 
                                      ifelse(categories %in% c("Enh", "EnhG"), "Enhancer", 
                                             ifelse(categories %in% c("TssBiv", "BivFlnk", "EnhBiv"), "BivReg", 
                                                    ifelse(categories %in% c("ReprPC", "ReprPCWk"), "ReprPoly", categories)
                                             )
                                      )
                               )
         ), 
         levels = c("TssProxProm", "TxFlnk", "ActTrans", "Enhancer", "ZNF.Rpts", "BivFlnk", "BivReg", "Het", "ReprPoly", "Quies")
         ),
         categories = factor(categories, levels = chromStates)) %>%
  ggplot(aes(x = categories, y = proportion*100, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_x_discrete(name = "ROADMAP chromatin states") +
  facet_wrap(~ Group, scales = "free_x") +  
  scale_y_continuous(name = "Proportion of CpGs (%)") +
  scale_fill_manual(name = "CpG Type", values = c("#000000", "#85C0F9", "#F5793A")) +
  theme_bw()

png("paper/chromStatesProp_model.png", width = 2500, height = 2000, res = 300)
chromSt_modComp_prop_plot
dev.off()

chromSt_enrich_compMods <- lapply(c("Common", "Main"), function(t){
  rbind(filter(methChromSt, Type == t),
        filter(methChromSt, Type == "non-eQTMs")) %>%
    select(ends_with("sum"), ends_with("sum2")) %>%
    g(chromStates, cols = c("_sum", "_sum2")) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:15) %>%
    mutate(Type = t)
})
chromSt_enrich_compMods <- Reduce(rbind, chromSt_enrich_compMods)



chrom_model <- as_tibble(methyAnnot) %>%
  filter(CpG %in% c(mainCpGs, comCpGs)) %>%
  mutate(Type = ifelse(CpG %in% mainCpGs, "Main", "Common")) %>%
  select(Type, eval(chromStates)) %>%
  group_by(Type) %>%
  summarize_at(chromStates, list(sum = sum, sum2 = sum2)) %>%
  arrange(desc(Type)) %>%
  g(chromStates, cols = c("_sum", "_sum2")) %>%
  as_tibble() %>%
  mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
  gather("Region", "Value", 1:15) %>%
  spread(par, Value) %>%
  mutate(Group = factor(ifelse(Region %in% c("TssA", "TssAFlnk"), "TssProxProm",
                               ifelse(Region %in% c("Tx", "TxWk"), "ActTrans", 
                                      ifelse(Region %in% c("Enh", "EnhG"), "Enhancer", 
                                             ifelse(Region %in% c("TssBiv", "BivFlnk", "EnhBiv"), "BivReg", 
                                                    ifelse(Region %in% c("ReprPC", "ReprPCWk"), "ReprPoly", Region)
                                             )
                                      )
                               )
  ), 
  levels = c("TssProxProm", "TxFlnk", "ActTrans", "Enhancer", "ZNF.Rpts", "BivFlnk", "BivReg", "Het", "ReprPoly", "Quies")
  )
  ) %>%
  ggplot(aes(x = Region, y = OR)) + 
  geom_bar(stat = "identity", position=position_dodge(), 
           fill = "#E69F00") + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "ROADMAP chromatin states") +
  facet_wrap(~ Group, scales = "free_x") +
  theme_bw() +
  ggtitle("Enrichment main vs common eQTMs") + 
  theme(plot.title = element_text(hjust = 0.5))



png("paper/enrich_chromStates_model.png", width = 2500, height = 1500, res = 300)
chrom_model
dev.off()



methyAnnot %>% 
  filter(CpG %in% c(mainCpGs, comCpGs)) %>%
  mutate(Type = factor(ifelse(CpG %in% mainCpGs, "Main", "Common"))) %>%
  glm(formula(paste("Type ~", paste(chromStates, collapse = "+"))), .,
     family = "binomial") %>%
  summary()

cot <- methyAnnot %>% 
  filter(CpG %in% c(mainCpGs, comCpGs)) %>%
  mutate(GeneRel = ifelse(UCSC_RefGene_Name == "", "Intergenic", "Genic"),
         Type = ifelse(CpG %in% mainCpGs, "Main", "Common")) 

methyAnnot %>% 
  filter(CpG %in% c(mainCpGs, comCpGs)) %>%
  mutate(GeneRel = ifelse(UCSC_RefGene_Name == "", "Intergenic", "Genic"),
         Type = ifelse(CpG %in% mainCpGs, "Main", "Common")) %>%
  dplyr::select(Enh, Type) %>%
  table()


mergeTB_Annot <- mergeTB %>%
  filter(sigPair.x) %>%
  dplyr::select(CpG, TC, sigPair.x, sigType) %>%
  left_join(dplyr::select(methyAnnot, CpG, UCSC_RefGene_Group, UCSC_RefGene_Name), by = "CpG") %>%
  left_join(dplyr::select(expAnnot, TC, GeneSymbol_Affy), by = "TC") %>%
  mutate(GeneAffy = strsplit(GeneSymbol_Affy, ";"))


mergeTB_Annot %>%
  mutate(match = sapply(seq_len(n()), function(x) any(GeneAffy[[x]] %in% UCSC_RefGene_Name[[x]]))) %>%
  group_by(sigType) %>%
  summarize(p = mean(match),
            n = sum(match))


## Run enrichment on specific TCs ####
### Define function to run enrichment
adjTCs <- setdiff(adjList$TC, cellList$TC)
adjGenes <- modU %>%
  group_by(TC) %>%
  summarize(sig = factor(
    ifelse(any(TC %in% adjTCs & sigPair), 1, 0)))
adjGOs <- computeGOs(adjGenes)

comTCs <- intersect(cellList$TC, adjList$TC)
comGenes <- modU %>%
  group_by(TC) %>%
  summarize(sig = factor(
    ifelse(any(TC %in% comTCs & sigPair), 1, 0)))
comGOs <- computeGOs(comGenes)

save(adjGOs, comGOs, file = "paper/GOobjectsCompModels.Rdata")

## It does not work on server. Run locally.
library(GOfuncR)
library(topGO)
library(dplyr)
server <- "//isg10174/data/WS_HELIX/HELIX_analyses/expr_met_SM/paper/"

load(paste0(server, "GOobjectsCompModels.Rdata"))

## Load categories and functions from eQTM_interpreation.R
## Select GOs with p-value < 0.001
mainMod <- addImmunityInfo(subset(adjGOs$table, as.numeric(w0) < 0.001))
comMod <- addImmunityInfo(subset(comGOs$table, as.numeric(w0) < 0.001))


tail(sort(sapply(top_GOs, function(x) sum(grepl(x, mainMod$parent)))))
tail(sort(sapply(top_GOs, function(x) sum(grepl(x, comMod$parent)))))

write.table(mainMod[, c("GO.ID", "GO_term", "w0", "parent", "immune")], 
            file = paste0(server, "/GOsMainModSp.txt"), 
            quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)
write.table(cellMod[, c("GO.ID", "GO_term", "w0", "parent", "immune")], 
            file = paste0(server, "/GOsCellModSp.txt"), 
            quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)




## Compare CpGs with Reinius ####
data(FlowSorted.Blood.450k.compTable)
adjCpGs <- setdiff(adjList$CpG, cellList$CpG)
comCpGs <- intersect(adjList$CpG, cellList$CpG)
cellCpGs <- setdiff(cellList$CpG, adjList$CpG)

data.frame(Fstat = FlowSorted.Blood.450k.compTable[c(adjCpGs, comCpGs), "Fstat"],
           Type = rep(c("Main", "Common"), c(length(adjCpGs), length(comCpGs)))) %>%
  lm(formula = log10(Fstat) ~ Type, data = .) %>%
  summary()
# Call:
#   lm(formula = log10(Fstat) ~ Type, data = .)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max
# -4.3185 -0.3776 -0.0170  0.4231  2.5337
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept) 0.841632   0.004882  172.40   <2e-16 ***
#   TypeMain    0.424345   0.007086   59.88   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6626 on 35063 degrees of freedom
# (163 observations deleted due to missingness)
# Multiple R-squared:  0.09279,   Adjusted R-squared:  0.09276
# F-statistic:  3586 on 1 and 35063 DF,  p-value: < 2.2e-16

png("paper/CompModelsCpGCellSpecific.png", width = 1500, height = 1000, res = 300)
methCell_compMods <- data.frame(Fstat = FlowSorted.Blood.450k.compTable[c(adjCpGs, comCpGs), "Fstat"],
           Type = factor(rep(c("Main", "Shared"), c(length(adjCpGs), length(comCpGs))),
                         levels = c("Shared", "Main"))) %>%
  ggplot(aes(y = log10(Fstat), x = Type, fill = Type)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(name = "") +
  scale_fill_manual(name = "", values = c("#85C0F9", "#F5793A")) +
  ggtitle("Methylation cell-type specificity") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")
methCell_compMods
dev.off()


## Compare TCs with Blueprint ####
bluep <- readRDS("data/lognormfpkm_allgenes.RDS")
bluepMap <- read.csv("data/blueprint_sampleID_mapping.csv", as.is = TRUE)

## Select cell types from venous blood and remove endothelial cells
cellMap <- subset(bluepMap, TISSUE_TYPE == "Venous blood" & GROUP != "ENDO")
comIDs <- intersect(cellMap$RUN_ID, colnames(bluep))
bluep_mat <- data.matrix(bluep)[, comIDs]
rownames(bluep_mat) <- bluep$sym
rownames(cellMap) <- cellMap$RUN_ID
modMat <- model.matrix(~ GROUP, cellMap[comIDs, ])

## Run linear model of genes vs cell type
fit <- lmFit(bluep_mat, modMat)
fit <- eBayes(fit)

## Select genes
mainTCs <- setdiff(adjList$TC, cellList$TC)
mainGenes <- filter(expAnnot, probeset_id %in% mainTCs)$GeneSymbol_Affy
mainGenes <- unlist(strsplit(mainGenes, ";"))

comTCs <- intersect(adjList$TC, cellList$TC)
comGenes <- filter(expAnnot, probeset_id %in% comTCs)$GeneSymbol_Affy
comGenes <- unlist(strsplit(comGenes, ";"))

## Remove genes regulated by different CpGs in both models
mainGenes.f <- setdiff(mainGenes, comGenes)
comGenes.f <- setdiff(comGenes, mainGenes)


### Make results table
bluep_res <- topTable(fit, n = Inf) %>% 
  mutate(gene = rownames(.),
         cat = ifelse(gene %in% mainGenes.f, "Main", 
                      ifelse(gene %in% comGenes.f, "Shared", "None")),
         cat = factor(cat, levels = c("Shared", "Main"))) %>%
  filter(cat != "None")

summary(lm(formula = log10(F) ~ cat, data = bluep_res))
  
png("paper/CompModelsGenesCellSpecific.png", width = 1500, height = 1000, res = 300)
exprCell_compMods <- bluep_res %>%
  ggplot(aes(y = log10(F), x = cat, fill = cat)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(name = "") +
  scale_y_continuous("log10(Fstat)") +
  scale_fill_manual(name = "", values = c("#85C0F9", "#F5793A")) +
  ggtitle("Gene expression cell-type specificity") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")
exprCell_compMods
dev.off()


png("paper/CompModelsCellSpecific.png", width = 2500, height = 1000, res = 300)
plot_grid(methCell_compMods, exprCell_compMods, nrow = 1, labels = "AUTO")
dev.off()


bluep_mat %>% group_by(mod) %>% summarize(m = median(log(v)))
## Los datos tienen muchos valores pequeños dp del log. Esto tira la media 
## muy para abajo en los genes especificos del main

png("paper/CompModelsTCCellSpecific.png", width = 1500, height = 1000, res = 300)
bluep_mat %>%
  ggplot(aes(x = mod, y = log(v), fill = mod)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(name = "")
dev.off()


## Compare with previous analysis
load("data/resultsList_filt.Rdata")
comRes <- datalist[intersect(comCpGs, names(datalist))]
intTCs <- intersect(comTCs, rownames(datalist[[1]]$crude) )

crudeTrans <- lapply(comRes, function(x) x$crude[intTCs, ]) %>%
  Reduce(f = rbind) %>%
  mutate(CpG = CPG, TC = trans) %>%
  left_join(dplyr::select(overDF, CpG, CpG_Pos) %>% distinct(), by = "CpG") %>% 
  left_join(dplyr::select(overDF, TC, TC_Pos, TC_strand) %>% distinct(), by = "TC") %>%
  mutate(Distance =  ifelse(TC_strand == "+", CpG_Pos - TC_Pos, TC_Pos - CpG_Pos))

mainRes <- datalist[intersect(mainCpGs, names(datalist))]
spTCs <- intersect(mainTCs, rownames(datalist[[1]]$crude) )

spTrans <- lapply(mainRes, function(x) x$crude[spTCs, ]) %>%
  Reduce(f = rbind) %>%
  mutate(CpG = CPG, TC = trans) %>%
  left_join(dplyr::select(overDF, CpG, CpG_Pos) %>% distinct(), by = "CpG") %>% 
  left_join(dplyr::select(overDF, TC, TC_Pos, TC_strand) %>% distinct(), by = "TC") %>%
  mutate(Distance =  ifelse(TC_strand == "+", CpG_Pos - TC_Pos, TC_Pos - CpG_Pos))


crudeTrans %>%
  mutate(type = "Common") %>%
  rbind(mutate(spTrans, type = "Main")) %>%
  ggplot(aes(x = Distance, y = -log10(P.Value))) +
  geom_point() +
  theme_bw() +
  facet_grid(type ~ .) +
#  ggtitle("Main model specific eQTMs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(-3e7, -2e7, -1e7, -1e6, 0, 1e6, 1e7, 2e7, 3e7), 
                     labels = c("-30Mb", "-20Mb", "-10Mb", "-1Mb", "0", "1Mb", "10Mb", "20Mb", "30Mb"))

cot <- spread(dplyr::select(spTrans, CpG, TC, P.Value), TC, P.Value)
pac <- cor(-log10(cot[, -1]))
h <- hclust(as.dist(1 - pac))

clas <- cutree(h, h = 0.1)
lapply(1:10, function(x){
  y <- pac[clas == x, clas == x]
  summary(y[upper.tri(y)])
})
## Select TCs with similar association with CpGS
tcs1 <- rownames(pac)[clas == 1]

filter(expAnnot, probeset_id %in% tcs1)$start

## Select CpGs associated with TCs with very low p-value
c1 <- filter(spTrans, TC %in% tcs1 & P.Value < 1e-10)$CpG
cpgs1 <- names(table(c1)[table(c1) == 5])
filter(spTrans,  TC %in% tcs1 & CpG %in% cpgs1)



## Compare with other eQTM studies (Bonder) PMID: 27918535 ####
bonder <- read.delim("data/Bonder_cis_eQTMs.txt")

## Remove CpGs and Genes not present in HELIX
bonder.f <- subset(bonder, SNPName %in% unique(modU$CpG))
bonder.f <- subset(bonder.f, HGNCName %in% expGenes)

## Add SYMBOL to significant pairs
sigDf_Annot <- sigDf %>%
  as_tibble() %>%
  left_join(dplyr::select(expAnnot, TC, GeneSymbol_Affy), by = "TC") %>%
  mutate(GeneAffy = strsplit(GeneSymbol_Affy, ";"))

## Bonder summary (after filtering)
length(unique(bonder.f$SNPName))
# [1] 10198
length(unique(bonder.f$HGNCName))
# [1] 3223

## Common CpGs
length(intersect(unique(bonder.f$SNPName), unique(sigDf_Annot$CpG)))
# [1] 4750
mean(unique(bonder.f$SNPName) %in% unique(sigDf_Annot$CpG))*100
# [1] 46.57776

## Common Genes
sigGenes <- unique(unlist(sigDf_Annot$GeneAffy))
length(intersect(unique(bonder.f$HGNCName), unique(sigDf_Annot$GeneAffy)))
# [1] 2003
mean(unique(bonder.f$HGNCName) %in% unique(sigDf_Annot$GeneAffy))*100
# [1] 62.14707


## Common Pairs
bonder.f$pair <- paste(bonder.f$SNPName, bonder.f$HGNCName)

sigPairs <- lapply(seq_len(nrow(sigDf_Annot)), function(x) {
  paste(sigDf_Annot[x, "CpG"], unlist(sigDf_Annot[x, "GeneAffy"]) )
})

length(intersect(bonder.f$pair, unlist(sigPairs)))
# [1] 4002
mean(bonder.f$pair %in% unlist(sigPairs))*100
# [1] 28.0783

## Restrict comparison to pairs present in HELIX dataset
allPairs <- lapply(seq_len(nrow(modU_Annot)), function(x) {
  paste(modU_Annot[x, "CpG"], unlist(modU_Annot[x, "GeneAffy"]) )
})

f <- function(x) mean(as.numeric(strsplit(as.character(x), ";")[[1]]))

bonder.merge <- bonder.f %>%
  mutate(CpG = SNPName, GeneSymbol_Affy = HGNCName,
         dir = sapply(IncludedDatasetsCorrelationCoefficient, f)) %>%
  left_join(sigDf_Annot, by = c("CpG", "GeneSymbol_Affy"))

table(bonder.merge$FC > 0, bonder.merge$dir > 0)
