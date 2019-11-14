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
library(topGO)
library(FlowSorted.Blood.450k)

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

## Distribution pairs per chromosomes ####
png("paper/eQTMs_Chr_distr.png", width = 2000, height = 1500, res = 300)
sigDf %>% 
  mutate(chr = substring(TC, 3, 4),
         chr = gsub("^0", "", chr)) %>%
  group_by(chr) %>%
  summarize(n = n()) %>%
  mutate(chr = factor(chr, levels = c(1:22, "X"))) %>%
  ggplot(aes(x = chr, y = n)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(name = "Chromosome") +
  scale_y_continuous(name = "Number Pairs") +
  theme_bw()
dev.off()


## Comparison distribution pairs per chromosome
chrDisteQTM <- sigDf %>% 
  mutate(chr = substring(TC, 3, 4),
         chr = gsub("^0", "", chr)) %>%
  group_by(chr) %>%
  summarize(n = n()) %>%
  mutate(chr = factor(chr, levels = c(1:22, "X")))

chrDistAll <- overDF %>% 
  mutate(chr = substring(TC, 3, 4),
         chr = gsub("^0", "", chr)) %>%
  group_by(chr) %>%
  summarize(n = n()) %>%
  mutate(chr = factor(chr, levels = c(1:22, "X", "Y")))

chrDistComb <- inner_join(chrDisteQTM, chrDistAll, by = "chr")

png("paper/chrDistr_all_eQTM.png", width = 2000, height = 1200, res = 300)
chrDistComb %>% mutate("All Pairs" = chrDistComb$n.y/sum(chrDistComb$n.y), 
                       "eQTMs" = chrDistComb$n.x/sum(chrDistComb$n.x)) %>% 
  gather(set, prop, 4:5) %>% 
  mutate(chr = factor(chr, levels = c(1:22, "X", "Y"))) %>%
  ggplot(aes( x = chr, y = prop, fill = set)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  theme_bw() +
  scale_fill_discrete(name = "") +
  scale_y_continuous(name = "Percentage of CpG-TC pairs") +
  scale_x_discrete(name = "Chromosome")
dev.off()
chisq.test(rbind(chrDistComb$n.x, chrDistComb$n.y))

chrDistComb %>% mutate(a = chrDistComb$n.y/sum(chrDistComb$n.y), 
                       e = chrDistComb$n.x/sum(chrDistComb$n.x),
                       da = e - a,
                       dr = e/a) %>%
  arrange(desc(abs(dr)))

chrDistComb %>% mutate(a = chrDistComb$n.y/sum(chrDistComb$n.y), 
                       e = chrDistComb$n.x/sum(chrDistComb$n.x),
                       da = e - a,
                       dr = e/a) %>%
  arrange(abs(dr))



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
  scale_y_continuous(name = "TC call rate")
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
t[1]/t[2]/t[3]*t[4]
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
  scale_color_discrete(name = "", labels = c("Non-significant", "Significant")) +
  scale_x_continuous(breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "", breaks = NULL) + 
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

# sigPair     m       l       h
# <lgl>   <int>   <dbl>   <dbl>
#   1 FALSE      14 -236288 236864
# 2 TRUE    -1333 -117011  84618.


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
         FC = abs(FC)/10) %>%
  group_by(Dir) %>%
  summarize(m = median(FC),
            l = quantile(FC, 0.25),
            h = quantile(FC, 0.75))

# Dir           m       l      h
# <chr>     <dbl>   <dbl>  <dbl>
#   1 Inverse  0.122 0.0655 0.240
# 2 Positive 0.112 0.0607 0.214

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
  geom_point() + 
  theme_bw() + 
  scale_color_manual(name = "", 
                     breaks = c("Inverse", "Positive"),
                     values = c("red", "blue")) +
  scale_x_continuous(breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "log2 FC/10% Methylation",
                     limits = c(-0.5, 0.5)) + 
  ggtitle("CpG-TC Distance vs Effect") + 
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
  left_join(dplyr::select(methyAnnot, CpG, UCSC_RefGene_Name), by = "CpG") %>%
  left_join(dplyr::select(expAnnot, TC, GeneSymbol_Affy), by = "TC") %>%
  mutate(GeneAffy = strsplit(GeneSymbol_Affy, ";"))

methGenes <- unique(unlist(modU_Annot$UCSC_RefGene_Name))
expGenes <- unique(unlist(modU_Annot$GeneAffy))
comGenes <- intersect(methGenes, expGenes)

modU_Annot_f <- modU_Annot %>%
  filter(sigPair)

methGenesf <- unique(unlist(modU_Annot_f$UCSC_RefGene_Name))
expGenesf <- unique(unlist(modU_Annot_f$GeneAffy))
comGenesf <- intersect(methGenesf, expGenesf)

modU_Annot_f %>%
  select(CpG, UCSC_RefGene_Name) %>%
  distinct() %>%
  summarize(n = n(),
            n_g = sum(UCSC_RefGene_Name != ""),
            n_GE = sum(sapply(UCSC_RefGene_Name, function(x) any(x %in% comGenes))))


modU_Annot_f2 <- modU_Annot_f %>%
  filter(sapply(UCSC_RefGene_Name, function(x) any(x %in% comGenes))) %>%
  mutate(match = GeneAffy %in% UCSC_RefGene_Name)


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

# Compare models ####
sigDf2 <- modC %>%
  as_tibble() %>%
  filter(sigPair)
## Statistics ####
### Significant pairs
nrow(sigDf2)
# [1] 39749

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
                    CellSp = unlist(Map(function(x, y) length(setdiff(x, y)), cellList, adjList)))
sink("paper/ModelsStatisticComparison.txt")
modCompDF[c(1, 2, 4, 3), ]
sink()

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
                        CellSp = unlist(Map(function(x, y) length(setdiff(x, y)), cellCpgs, adjCpgs)))
sink("paper/ModelsStatisticComparisonCpG.txt")
modCompDF2[c(1, 2, 4, 5, 3, 6), ]
sink()

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

## Volcano plot
png("paper/volcano_all_adjCells.png", width = 2000, height = 2000, res = 300)
volcano_plot(modC$p.value, modC$FC/100, paste(modC$CpG, modC$TC), 
             tPV = -log10(1e-8), tFC = 0.01, show.labels = FALSE) +
  geom_point(alpha = 0.1)
dev.off()



## Compare estimates ####
### Merge dataset
mergeTB <- modU %>%
  left_join(modC, by = c("CpG", "TC")) %>%
  as_tibble() %>%
  filter(sigPair.x == TRUE | sigPair.y == TRUE) %>%
  mutate(sigType = ifelse(sigPair.x == TRUE, ifelse(sigPair.y == TRUE, "Both", "Adjusted"), "Cell"))


png("paper/CompModelsP_values.png", width = 2500, height = 2500, res = 300)
ggplot(mergeTB, aes(x = -log10(p.value.y), y = -log10(p.value.x), col = sigType)) +
  geom_point() +
  scale_x_continuous(name = "Cell adjusted") + 
  scale_y_continuous("Adjusted") + 
  ggtitle("-log10 p-values comparative") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(name = "")
dev.off()

# mergeTB <- mergeTB %>%
#   mutate(Diff.pval = pnorm(abs(FC.x - FC.y)/sqrt(SD.x**2 + SD.y**2), lower.tail = FALSE), 
#          isDiff = ifelse(Diff.pval < 0.05, "Different", "Equal"))

png("paper/CompModelsEstimates.png", width = 2500, height = 2500, res = 300)
ggplot(mergeTB, aes(x = FC.y, y = FC.x, col = sigType)) +
  geom_point() +
  scale_x_continuous(name = "Cell adjusted") + 
  scale_y_continuous("Adjusted") + 
  ggtitle("Estimates comparative") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(name = "")
dev.off()


## Run enrichment on specific TCs ####
### Define function to run enrichment
computeGOs <- function(df){
  
  genesv <- df$sig
  names(genesv) <- df$TC
  
  Data <- new("topGOdata", 
              description = "GO analysis of TCs",
              ontology = "BP",
              allGenes = genesv,
              annot = annFUN.db,
              nodeSize = 10,
              affyLib = "hta20transcriptcluster.db")
  
  class <- runTest(Data, algorithm = "classic", statistic = "fisher")
  elim <- runTest(Data, algorithm = "elim", statistic = "fisher")
  weigth <- runTest(Data, algorithm = "weight", statistic = "fisher")
  mixed <- runTest(Data, algorithm = "weight01", statistic = "fisher")
  
  
  finTab <- GenTable(Data, 
                     classic = class, 
                     elim = elim, 
                     weight = weigth, 
                     w0 = mixed,
                     orderBy = "weight",
                     ranksOf = "elim",
                     topNodes = length(score(class)))
  gos <- list(classic = class, elim = elim,
              weigth = weigth, 
              mixed = mixed)
  list(gos = gos, table = finTab)
}
adjTCs <- setdiff(adjList$TC, cellList$TC)
adjGenes <- modU %>%
  group_by(TC) %>%
  summarize(sig = factor(
    ifelse(any(TC %in% adjTCs & sigPair), 1, 0)))
adjGOs <- computeGOs(adjGenes)

a <- adjGOs$table %>%
  filter(w0 < 0.01 & classic < 0.01) %>%
  arrange(w0)
write.table(a[, c("GO.ID", "Term", "w0", "classic")], file = "paper/GOsAdjSp.txt", 
            quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)

cellTCs <- setdiff(cellList$TC, adjList$TC)
cellGenes <- modC %>%
  group_by(TC) %>%
  summarize(sig = factor(
    ifelse(any(TC %in% cellTCs & sigPair), 1, 0)))
cellGOs <- computeGOs(cellGenes)

a <- cellGOs$table %>%
  filter(w0 < 0.01 & classic < 0.01) %>%
  arrange(w0)
write.table(a[, c("GO.ID", "Term", "w0", "classic")], file = "paper/GOsCellSp.txt", 
            quote = FALSE, col.names = TRUE,
            sep = "\t", row.names = FALSE)


## Compare CpGs with Reinius ####
data(FlowSorted.Blood.450k.compTable)
adjCpGs <- setdiff(adjList$CpG, cellList$CpG)
cellCpGs <- setdiff(cellList$CpG, adjList$CpG)

data.frame(Fstat = FlowSorted.Blood.450k.compTable[c(adjCpGs, cellCpGs), "Fstat"],
           Type = rep(c("Adj", "Cell"), c(length(adjCpGs), length(cellCpGs)))) %>%
  lm(formula = log10(Fstat) ~ Type, data = .) %>%
  summary()
  
# Call:
#   lm(formula = log10(Fstat) ~ Type, data = .)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max
# -3.6491 -0.3193 -0.0387  0.3510  2.1680
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  1.26598    0.00460  275.21   <2e-16 ***
#   TypeCell     0.16409    0.01114   14.73   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5934 on 20061 degrees of freedom
# (73 observations deleted due to missingness)
# Multiple R-squared:  0.0107,    Adjusted R-squared:  0.01065
# F-statistic: 216.9 on 1 and 20061 DF,  p-value: < 2.2e-16
# 

png("paper/CompModelsCpGCellSpecific.png", width = 1500, height = 1000, res = 300)
data.frame(Fstat = FlowSorted.Blood.450k.compTable[c(adjCpGs, cellCpGs), "Fstat"],
           Type = rep(c("Adjusted", "Cell"), c(length(adjCpGs), length(cellCpGs)))) %>%
  ggplot(aes(y = log10(Fstat), x = Type, fill = Type)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(name = "")
dev.off()


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
