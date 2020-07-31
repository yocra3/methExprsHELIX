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
library(eulerr)

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
  mutate(sigVar = ifelse(sig == "random", "CpGs not in eQTMs", "CpGs in eQTMs")) %>%
  ggplot(aes(x = sigVar, y = meth_range, fill = sigVar)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(name = "") +
  theme(legend.position = "none") +
  scale_y_continuous(name = "Methylation range") +
  scale_fill_manual(values = c("#999999", "#FFFFFF"))
dev.off()



### TC call rate
int_TC <- expAnnot %>%
  dplyr::select(probeset_id, Expressed, CallRate) %>%
  mutate(sig = ifelse(probeset_id %in% sigTCs, "TCs in eQTMs", "TCs not in eQTMs"))
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
  mutate(Direction = ifelse(!sigPair, "Non-eQTM", ifelse(FC > 0, "Positive", "Inverse")),
         Direction = factor(Direction, levels = c("Inverse", "Positive", "Non-eQTM"))) %>%
  ggplot(aes(x = Distance, color = Direction)) + geom_density() + 
  theme_bw() + 
  scale_color_manual(name = "eQTM type", values = c("#E69F00", "#009E73", "#000000")) +
  scale_x_continuous(name = "CpG-TC TSS distance", 
                     breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "Density") + 
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
  scale_color_manual(name = "eQTM type", 
                     breaks = c("Inverse", "Positive"),
                     values = c("#E69F00", "#009E73")) +
  scale_x_continuous(name = "CpG-TC TSS distance",
                     breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "log2 FC/0.1 Methylation") + 
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
  p.val = chisq.test(t)$p.value
  or <-  t[1]/t[2]/t[3]*t[4]
  ORl <- log(or)
  SEl <- sqrt(1/t[1, 1] + 1/t[1, 2] + 1/t[2, 1] + 1/t[2, 2])
  c(OR = or, p.val = p.val, ORm = exp(ORl - 1.96*SEl), 
    ORM = exp(ORl + 1.96*SEl))
  
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


ORsInv <- lapply(gpos, function(x) {
  cot <- filter(modU_gene_pairs, modU_gene_pairs$sigPair == FALSE | FC < 0)
  t <- table(cot$Gene_Group %in% x, cot$sigPair)
  p.val = chisq.test(t)$p.value
  or <-  t[1]/t[2]/t[3]*t[4]
  ORl <- log(or)
  SEl <- sqrt(1/t[1, 1] + 1/t[1, 2] + 1/t[2, 1] + 1/t[2, 2])
  c(OR = or, p.val = p.val, ORm = exp(ORl - 1.96*SEl), 
    ORM = exp(ORl + 1.96*SEl))
}) %>%
  Reduce(f = cbind)
colnames(ORsInv) <- gpos  

ORsPos <- lapply(gpos, function(x) {
  cot <- filter(modU_gene_pairs, modU_gene_pairs$sigPair == FALSE | FC > 0)
  t <- table(cot$Gene_Group %in% x, cot$sigPair)
  p.val = chisq.test(t)$p.value
  or <-  t[1]/t[2]/t[3]*t[4]
  ORl <- log(or)
  SEl <- sqrt(1/t[1, 1] + 1/t[1, 2] + 1/t[2, 1] + 1/t[2, 2])
  c(OR = or, p.val = p.val, ORm = exp(ORl - 1.96*SEl), 
    ORM = exp(ORl + 1.96*SEl))
}) %>%
  Reduce(f = cbind)
colnames(ORsPos) <- gpos  

convertORmat <- function(OR){
  OR %>%
    data.frame() %>%
    mutate(vars = c("OR", "p-value", "ORm", "ORM")) %>%
    gather(Region, value, 1:6) %>%
    mutate(Region = gsub("X", "", Region),
           Region = gsub(".", "'", Region, fixed = TRUE))
  
}

png("paper/CpGEnrichGenePosition.png", width = 3000, height = 2000, res = 300)
rbind(convertORmat(ORs) %>% mutate(type = "All eQTMs"),
      convertORmat(ORsInv) %>% mutate(type = "Inverse"),
      convertORmat(ORsPos) %>% mutate(type = "Positive")) %>%
  spread(vars, value) %>%
  mutate(Region = factor(Region, levels = gpos)) %>%
  ggplot(aes(x = type, y = OR, fill = Region)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORm, ymax = ORM)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "eQTM type") +
  scale_fill_discrete(name = "Gene relative position") +
  theme_bw() 
dev.off()

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

### Merge dataset ####
mergeTB <- modU %>%
  left_join(modC, by = c("CpG", "TC")) %>%
  as_tibble() %>%
  filter(sigPair.x == TRUE | sigPair.y == TRUE) %>%
  mutate(sigType = ifelse(sigPair.x == TRUE, ifelse(sigPair.y == TRUE, "Model-shared", "Unadj. cell-specific"), "Adj. cell-specific"))


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
  mutate(sigType = factor(sigType, levels = c("Model-shared", "Unadj. cell-specific", "Adj. cell-specific"))) %>%
  ggplot(aes(x = -log10(p.value.y), y = -log10(p.value.x), col = sigType)) +
  geom_point() +
  scale_x_continuous(name = "Cell adjusted") + 
  scale_y_continuous("Main model") + 
  ggtitle("Comparison of -log10 p-values") +
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
  mutate(sigType = factor(sigType, levels = c("Model-shared", "Unadj. cell-specific"))) %>%
  ggplot(aes(x = Distance, color = sigType)) +
  geom_density() +
  theme_bw() +
  scale_color_manual(name = "eQTM type", values = c("#85C0F9", "#F5793A")) +
  scale_x_continuous(name = "CpG-TC TSS distance", 
                     breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "Density") + 
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


tmp <- mergeTB %>%
  filter(sigPair.x) %>%
  inner_join(overDF, by = c("TC", "CpG")) %>%
  mutate(mod = ifelse(sigType == "Both", "Common", sigType))
  ks.test(subset(tmp, mod == "Shared")$Distance,
              subset(tmp, mod == "Main")$Distance)


estim_comp <- mergeTB %>% 
  mutate(sigType = factor(sigType, levels = c("Model-shared", "Unadj. cell-specific", "Adj. cell-specific"))) %>%
  ggplot(aes(x = FC.y/10, y = FC.x/10, col = sigType)) +
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed") +
  scale_x_continuous(name = "Cell adjusted") + 
  scale_y_continuous("Main model") + 
  ggtitle("Comparison of effect size stimates") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_color_manual(name = "", values = c("#85C0F9", "#F5793A", "#0F2080")) +
  facet_wrap(. ~ sigType) +
  geom_smooth(method = "lm") 


png("paper/CompModelsEstimatesPvals.png", width = 2500, height = 2000, res = 300)
plot_grid(estim_comp, pvals_comp, nrow = 2, labels = c("A", "B"))
dev.off()

filter(mergeTB, sigType == "Main") %>% lm(FC.y ~ FC.x, .) %>% summary()
filter(mergeTB, sigType == "Cell") %>% lm(FC.y ~ FC.x, .) %>% summary()
filter(mergeTB, sigType == "Both") %>% lm(FC.y ~ FC.x, .) %>% summary()

## Compare assocations strength ####
### FC
mergeTB %>%
  filter(sigPair.x) %>%
  group_by(sigType) %>%
  summarize(m = median(abs(FC.x)/10))

wilcox.test(abs(subset(mergeTB, sigPair.x & sigType == "Shared")$FC.x)/10,
        abs(subset(mergeTB, sigPair.x & sigType == "Unadj. Cell")$FC.x)/10, conf.int = TRUE)

mergeTB %>%
  filter(sigPair.x) %>%
  group_by(sigType) %>%
  summarize(m = median(SD.x))

wilcox.test(subset(mergeTB, sigPair.x & sigType == "Shared")$SD.x,
        subset(mergeTB, sigPair.x & sigType == "Unadj. Cell")$SD.x, conf.int = TRUE)


mergeTB %>%
  filter(sigPair.x) %>%
  group_by(sigType) %>%
  summarize(m = median(-log10(p.value.x)))

wilcox.test(-log10(subset(mergeTB, sigPair.x & sigType == "Shared")$p.value.x),
            -log10(subset(mergeTB, sigPair.x & sigType == "Unadj. Cell")$p.value.x), conf.int = TRUE)

## Run enrichment on CpGs ####
mainCpGs <- setdiff(adjList$CpG, cellList$CpG)
comCpGs <- intersect(adjList$CpG, cellList$CpG)

methyAnnot %>% 
  filter(CpG %in% c(mainCpGs, comCpGs)) %>%
  mutate(GeneRel = ifelse(UCSC_RefGene_Name == "", "Intergenic", "Genic"),
         Type = ifelse(CpG %in% mainCpGs, "Unadj. Cell", "Shared")) %>%
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
  mutate(Type = ifelse(CpG %in% mainCpGs, "CpGs in unadj. cell-specific eQTMs", ifelse(CpG %in% comCpGs, "CpGs in model-shared eQTMs", "CpGs not in eQTMs"))) %>%
  group_by(Type) %>%
  summarize_at(chromStates, list(sum = sum, sum2 = sum2, prop = mean))


chromSt_modComp_prop_plot <- methChromSt %>%
  select(Type, ends_with("prop")) %>%
  gather(categories, proportion, 2:(2+length(chromStates) - 1)) %>%
  ungroup() %>%
  mutate(Type = factor(Type,  levels = c("CpGs not in eQTMs", "CpGs in model-shared eQTMs", "CpGs in unadj. cell-specific eQTMs")),
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
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  scale_x_discrete(name = "ROADMAP chromatin states") +
  facet_wrap(~ Group, scales = "free_x") +  
  scale_y_continuous(name = "Proportion of CpGs (%)") +
  scale_fill_manual(name = "CpG type", values = c("#FFFFFF", "#85C0F9", "#F5793A")) +
  theme_bw()

png("paper/chromStatesProp_model.png", width = 2500, height = 2000, res = 300)
chromSt_modComp_prop_plot
dev.off()

chromSt_enrich_compMods <- lapply(c("CpGs in model-shared eQTMs", "CpGs in unadj. cell-specific eQTMs"), function(t){
  rbind(filter(methChromSt, Type == t),
        filter(methChromSt, Type == "CpGs not in eQTMs")) %>%
    select(ends_with("sum"), ends_with("sum2")) %>%
    g(chromStates, cols = c("_sum", "_sum2")) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:15) %>%
    mutate(Type = t)
})
chromSt_enrich_compMods <- Reduce(rbind, chromSt_enrich_compMods)


chrom_model <- chromSt_enrich_compMods %>%
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
  ),
  ) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge(), color = "black") + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "ROADMAP chromatin states") +
  scale_fill_manual(name = "CpG type", values = c("#85C0F9", "#F5793A")) +
  facet_wrap(~ Group, scales = "free_x") +
  theme_bw() +
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
           Type = factor(rep(c("CpGs in unadj. cell-specific eQTMs", "CpGs in model-shared eQTMs"), c(length(adjCpGs), length(comCpGs))),
                         levels = c("CpGs in model-shared eQTMs", "CpGs in unadj. cell-specific eQTMs"))) %>%
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
         cat = ifelse(gene %in% mainGenes.f, "TCs in unadj. cell-specific eQTMs", 
                      ifelse(gene %in% comGenes.f, "TCs in model-shared eQTMs", "None")),
         cat = factor(cat, levels = c("TCs in model-shared eQTMs", "TCs in unadj. cell-specific eQTMs"))) %>%
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


png("paper/CompModelsCellSpecific.png", width = 1500, height = 1500, res = 300)
plot_grid(methCell_compMods, exprCell_compMods, ncol = 1, labels = "AUTO")
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

## Age variability (MEDALL - Xu et al) ####
agedf <- read.csv2("data/DiffMethyAgeCpGs.csv", as.is = TRUE)
agedf <- agedf %>%
  as_tibble() %>%
  mutate(Dir = ifelse(as.numeric(beta8.avg) > as.numeric(beta0.avg), "Increased", "Decreased"), 
         CpG = ILMNID)

ageSum <- CpGsSum %>%
  dplyr::select(CpG, Direction) %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both", "Non-significant"))) %>%
  left_join(dplyr::select(agedf, Dir, CpG), by = "CpG") %>%
  mutate(Dir = ifelse(is.na(Dir), "Constant", Dir)) %>%
  group_by(Direction, Dir) %>%
  summarize(n = n()) %>%
  spread(Dir, n) %>%
  ungroup() %>%
  rbind(filter(., Direction != "Non-significant") %>% 
          summarize_all(fsum)) %>%
  mutate(Direction = as.character(Direction),
         Direction = ifelse(is.na(Direction), "eQTM", Direction),
         Variable = Decreased + Increased,
         tot = Variable + Constant)

ageG <- c("Variable", "Decreased", "Increased")


typesAge <- lapply(types, function(t){
  rbind(filter(ageSum, Direction == t),
        filter(ageSum, Direction == "Non-significant")) %>%
    g2(ageG) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:3) %>%
    mutate(Type = t)
})
combAge <- Reduce(rbind, typesAge)
combAge$dataset <- "MeDALL"

age_var <- combAge %>%
  spread(par, Value) %>%
  mutate(Type = ifelse(Type == "eQTM", "All", Type),
         Type = factor(Type, levels = c("All", "Inverse", "Positive"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge(), color = "black") + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "Methylation during childhood", limits = ageG) +
  scale_fill_manual(name = "CpG type", values = c("#999999", "#E69F00", "#009E73"),
                    labels = c("CpGs in eQTMs", "CpGs in inverse eQTMs", "CpGs in positive eQTMs")) +
  theme_bw() 

png("paper/CpGEnrichAge.png", width = 3000, height = 2000, res = 300)
age_var
dev.off()

tmp <- CpGsSum %>%
  dplyr::select(CpG, Direction) %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both", "Non-significant"))) %>%
  left_join(dplyr::select(agedf, Dir, CpG), by = "CpG") %>%
  mutate(Dir = ifelse(is.na(Dir), "Constant", Dir)) 

png("paper/CpGEnrichAge_methLevels.png", width = 3000, height = 2000, res = 300)
age_meth <-  as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, median) %>%
  right_join(tmp, by = "CpG") %>%
  mutate(Significant = ifelse(Direction == "Non-significant", "CpGs not in eQTMs", "CpGs in eQTMs")) %>%
  ggplot(aes(x = Dir, y = median, fill = Significant)) +
  geom_violin() +
  scale_x_discrete(name = "Methylation during childhood") +
  scale_y_continuous(name = "Median methylation") +
  scale_fill_manual(name = "CpG type", values = c("#999999", "#FFFFFF")) +
  theme_bw()
age_meth
dev.off()

## Age variability (Mulder) ####
epideltadf <- read.table("data/epidelta_2020-07-17.txt", as.is = TRUE)
epideltadf$CpG <- rownames(epideltadf)
epideltadf <- epideltadf %>%
  as_tibble() %>%
  mutate(Dir = ifelse(M1.change.p > 1e-7, "Constant",
                      ifelse(M1.change.estimate > 0, "Increased", "Decreased")))

age_cont_Sum <- CpGsSum %>%
  dplyr::select(CpG, Direction) %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both", "Non-significant"))) %>%
  left_join(dplyr::select(epideltadf, Dir, M1.change.estimate, CpG), by = "CpG") %>%
  group_by(Direction, Dir) %>%
  summarize(n = n()) %>%
  spread(Dir, n) %>%
  ungroup() %>%
  rbind(filter(., Direction != "Non-significant") %>% 
          summarize_all(fsum)) %>%
  mutate(Direction = as.character(Direction),
         Direction = ifelse(is.na(Direction), "eQTM", Direction),
         Variable = Decreased + Increased,
         tot = Variable + Constant)

ageG <- c("Variable", "Decreased", "Increased")


types_cont_Age <- lapply(types, function(t){
  rbind(filter(age_cont_Sum, Direction == t),
        filter(age_cont_Sum, Direction == "Non-significant")) %>%
    g2(ageG) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:3) %>%
    mutate(Type = t)
})
comb_cont_Age <- Reduce(rbind, types_cont_Age)
comb_cont_Age$dataset <- "Epidelta"


png("paper/CpGEnrichAge_Mulder.png", width = 3000, height = 2000, res = 300)
age_cont_var <- comb_cont_Age %>%
  spread(par, Value) %>%
  mutate(Type = ifelse(Type == "eQTM", "All", Type),
         Type = factor(Type, levels = c("All", "Inverse", "Positive"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge(), color = "black") + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "Methylation during childhood", limits = ageG) +
  scale_fill_manual(name = "CpG type", values = c("#999999", "#E69F00", "#009E73"),
                    labels = c("CpGs in eQTMs", "CpGs in inverse eQTMs", "CpGs in positive eQTMs")) +
  theme_bw() 
age_cont_var
dev.off()


tmp_cont <- CpGsSum %>%
  dplyr::select(CpG, Direction) %>%
  mutate(Direction = factor(Direction, levels = c("Inverse", "Positive", "Both", "Non-significant"))) %>%
  left_join(dplyr::select(epideltadf, Dir, CpG), by = "CpG") %>%
  mutate(Dir = ifelse(is.na(Dir), "Constant", Dir)) 

age_meth <-  as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, median) %>%
  right_join(tmp_cont, by = "CpG") %>%
  mutate(Significant = ifelse(Direction == "Non-significant", "Non-eQTM", "eQTM")) %>%
  ggplot(aes(x = Dir, y = median, fill = Significant)) +
  geom_violin() +
  scale_x_discrete(name = "Methylation during childhood") +
  scale_y_continuous(name = "Median methylation") +
  scale_fill_manual(name = "", values = c("#999999", "#FFFFFF")) +
  theme_bw()

png("paper/CpGEnrichAge_methLevels.png", width = 3000, height = 2000, res = 300)
age_meth
dev.off()

## Make combined plot
age_var_comb <- rbind(comb_cont_Age, combAge) %>%
  mutate(dataset = factor(dataset, levels = c("MeDALL", "Epidelta"))) %>%
  spread(par, Value) %>%
  mutate(Type = ifelse(Type == "eQTM", "All", Type),
         Type = factor(Type, levels = c("All", "Inverse", "Positive"))) %>%
  ggplot(aes(x = Region, y = OR, fill = Type)) + 
  geom_bar(stat = "identity", position=position_dodge(), color = "black") + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "Change in methylation", limits = ageG) +
  scale_fill_manual(name = "CpG type", values = c("#999999", "#E69F00", "#009E73"),
                    labels = c("CpGs in eQTMs", "CpGs in inverse eQTMs", "CpGs in positive eQTMs")) +
  theme_bw() +
  facet_grid(~ dataset)

## Overlap age variability datasets
merged_age <- left_join(select(agedf, CpG, Dir), 
                        select(epideltadf, CpG, Dir, M1.change.p), 
                        by = "CpG")
mean(merged_age$Dir.x == merged_age$Dir.y)     

png("paper/CpGEnrich_ageVar_methLevels.png", width = 3000, height = 1500, res = 300)
as_tibble(methyAnnot) %>%
  mutate(CpG = Name) %>%
  select(CpG, median) %>%
  left_join(select(CpGsSum, CpG, Type)) %>%
  mutate(Significant = ifelse(Type == "Non-significant", "CpGs not in eQTMs", "CpGs in eQTMs")) %>%
  left_join(dplyr::select(agedf, Dir, CpG), by = "CpG") %>%
  mutate(Xu = ifelse(is.na(Dir), "Constant", Dir))  %>%
  select(-Dir) %>%
  left_join(dplyr::select(epideltadf, Dir, CpG), by = "CpG") %>%
  mutate(Epidelta = Dir) %>%
  filter(!is.na(Significant) & !is.na(Epidelta)) %>%
  select(Xu, Epidelta, CpG, median, Significant) %>%
  gather(Database, Status, 1:2) %>%
  mutate(Database = ifelse(Database == "Xu", "MeDALL", "Epidelta"),
         Database = factor(Database, levels = c("MeDALL", "Epidelta"))) %>%
  ggplot(aes(x = Status, y = median, fill = Significant)) +
  geom_violin() +
  scale_x_discrete(name = "Change in methylation") +
  scale_y_continuous(name = "Median methylation") +
  scale_fill_manual(name = "CpG type", values = c("#999999", "#FFFFFF")) +
  facet_grid(~ Database) +
  theme_bw()
dev.off()




## Compare with other eQTM studies (Kennedy) PMID:  29914364 ####
eQTM_blood <- read.delim("data/Kennedy_eQTM_catalogue_wholeBlood.txt", as.is = TRUE)
eQTM_monocytes <- read.delim("data/Kennedy_eQTM_catalogue_Monocytes.txt", as.is = TRUE)

## Remove CpGs and Genes not present in HELIX
eQTM_blood.f <- subset(eQTM_blood, CpG.probe %in% unique(modU$CpG) & 
                         annot.gene %in% expGenes)
eQTM_monocytes.f <- subset(eQTM_monocytes, CpG.probe %in% unique(modU$CpG) & 
                         annot.gene %in% expGenes)

## Compare with main model
blood.merge <- eQTM_blood.f %>%
  mutate(CpG = CpG.probe, GeneSymbol_Affy = annot.gene) %>%
  inner_join(modU_Annot, by = c("CpG", "GeneSymbol_Affy")) %>%
  left_join(select(expAnnot, c("TC", "start", "stop")), by = "TC") %>%
  mutate(sigType = ifelse(sigPair, ifelse(p.val < 1e-11, "Both", "Main"), 
                          ifelse(p.val < 1e-11, "WholeBlood", "None")),
         sigType = factor(sigType, levels = c("Both", "None", "Main", "WholeBlood")),
         sigType2 = ifelse(sigPair, "Shared", "WholeBlood"))

monocytes.merge <- eQTM_monocytes.f %>%
  mutate(CpG = CpG.probe, GeneSymbol_Affy = annot.gene) %>%
  inner_join(modU_Annot, by = c("CpG", "GeneSymbol_Affy")) %>%
  mutate(sigType = ifelse(sigPair, ifelse(p.val < 1e-11, "Both", "Main"), 
                          ifelse(p.val < 1e-11, "Monocytes", "None")),
         sigType = factor(sigType, levels = c("Both", "None", "Main", "Monocytes")),
         sigType2 = ifelse(sigPair, "Shared", "Monocytes"))

blood.merge %>% group_by(sigType) %>% summarize(rFC = cor(FC, beta), 
                                                pFC = cor.test(FC, beta)$p.value,
                                                signConc = mean(sign(FC) == sign(beta)))
blood.merge %>% group_by(sigType2) %>% summarize(rFC = cor(FC, beta), 
                                                pFC = cor.test(FC, beta)$p.value,
                                                signConc = mean(sign(FC) == sign(beta)))

cor(monocytes.merge$FC, monocytes.merge$beta)
cor(-log10(monocytes.merge$p.val + 1e-100), -log10(monocytes.merge$p.value))
monocytes.merge %>% group_by(sigType) %>% summarize(rFC = cor(FC, beta), 
                                                pFC = cor.test(FC, beta)$p.value,
                                                signConc = mean(sign(FC) == sign(beta)))
monocytes.merge %>% group_by(sigType2) %>% summarize(rFC = cor(FC, beta), 
                                                 pFC = cor.test(FC, beta)$p.value,
                                                 signConc = mean(sign(FC) == sign(beta)))
table(blood.merge$sigType)
table(monocytes.merge$sigType)

table(blood.merge$sigType2)
table(monocytes.merge$sigType2)

## Map one pair of HELIX (lowest p-value) to each pair in GTP and MESA
blood.merge.uniq <- blood.merge %>% 
  group_by(exp.Probe, CpG.probe) %>%
  top_n(-1, p.value)

monocytes.merge.uniq <- monocytes.merge %>% 
  group_by(exp.Probe, CpG.probe) %>%
  top_n(-1, p.value)


blood.merge.uniq %>% group_by(sigType) %>% summarize(rFC = cor(FC, beta), 
                                                pFC = cor.test(FC, beta)$p.value,
                                                signConc = mean(sign(FC) == sign(beta)))
blood.merge.uniq %>% group_by(sigType2) %>% summarize(rFC = cor(FC, beta), 
                                                 pFC = cor.test(FC, beta)$p.value,
                                                 signConc = mean(sign(FC) == sign(beta)))


monocytes.merge.uniq %>% group_by(sigType) %>% summarize(rFC = cor(FC, beta), 
                                                    pFC = cor.test(FC, beta)$p.value,
                                                    signConc = mean(sign(FC) == sign(beta)))
monocytes.merge.uniq %>% group_by(sigType2) %>% summarize(rFC = cor(FC, beta), 
                                                     pFC = cor.test(FC, beta)$p.value,
                                                     signConc = mean(sign(FC) == sign(beta)))
table(blood.merge.uniq$sigType2)
table(monocytes.merge.uniq$sigType2)



#### Whole blood plots
pWBFC <- blood.merge %>%
  ggplot(aes(x = FC/10, y = beta, color = sigType)) +
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed") +
  scale_x_continuous(name = "Main model") + 
  scale_y_continuous("Kennedy (Whole Blood)") + 
  ggtitle("Estimates comparative") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
#  scale_color_manual(name = "", values = c("#85C0F9", "#F5793A", "#0F2080")) +
  facet_wrap(. ~ sigType) +
  geom_smooth(method = "lm") 


pWBpval <- blood.merge %>%
  ggplot(aes(x = -log10(p.value), y = -log10(p.val), color = sigType)) +
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed") +
  scale_x_continuous(name = "Main model") + 
  scale_y_continuous("Kennedy (Whole Blood)") + 
  ggtitle("-log10 p-value comparative") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  #  scale_color_manual(name = "", values = c("#85C0F9", "#F5793A", "#0F2080")) +
  facet_wrap(. ~ sigType)

#### Monocyte plots
pMONFC <- monocytes.merge %>%
  ggplot(aes(x = FC/10, y = beta, color = sigType)) +
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed") +
  scale_x_continuous(name = "Main model") + 
  scale_y_continuous("Kennedy (Monocytes)") + 
  ggtitle("Estimates comparative") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  #  scale_color_manual(name = "", values = c("#85C0F9", "#F5793A", "#0F2080")) +
  facet_wrap(. ~ sigType) +
  geom_smooth(method = "lm") 


pMONpval <- monocytes.merge %>%
  ggplot(aes(x = -log10(p.value), y = -log10(p.val), color = sigType)) +
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed") +
  scale_x_continuous(name = "Main model") + 
  scale_y_continuous("Kennedy (Monocytes)") + 
  ggtitle("-log10 p-value comparative") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  #  scale_color_manual(name = "", values = c("#85C0F9", "#F5793A", "#0F2080")) +
  facet_wrap(. ~ sigType)


png("paper/ComparisonKennedyEQTM.png", width = 2500, height = 2000, res = 300)
plot_grid(pWBFC, pWBpval, pMONFC, pMONpval, nrow = 2, ncol = 2, labels = c("A", "", "B", ""))
dev.off()




## Comparison catalogue pairs in Kennedy ####
monocytes.merge$pair <- paste(monocytes.merge$CpG, monocytes.merge$TC)
blood.merge$pair <- paste(blood.merge$CpG, blood.merge$TC)
modU_Annot$MESA <- ifelse(paste(modU_Annot$CpG, modU_Annot$TC) %in% monocytes.merge$pair, 
                          "MESA", "not in MESA")
modU_Annot$GTP <- ifelse(paste(modU_Annot$CpG, modU_Annot$TC) %in% blood.merge$pair, 
                          "GTP", "not in GTP")
modU_Annot$adults <- ifelse(modU_Annot$sigPair, 
                           ifelse(modU_Annot$MESA == "MESA" | modU_Annot$GTP == "GTP", "Shared", "Children"), 
                           "Other")

modU_Annot %>% 
  filter(adults != "Other") %>%
  group_by(adults) %>%
  summarize(nCpG = length(unique(CpG)),
            nTC = length(unique(TC)))


modU_Annot %>% filter(sigPair == TRUE) %>% group_by(GTP) %>% 
  summarize(m = median(-log10(p.value)))
modU_Annot %>% filter(sigPair == TRUE) %>% group_by(MESA) %>% 
  summarize(m = median(-log10(p.value)))
modU_Annot %>% filter(sigPair == TRUE) %>% group_by(adults) %>% 
  summarize(m = median(-log10(p.value)))

modU_Annot %>% filter(sigPair == TRUE) %>% group_by(MESA) %>% 
  summarize(m = median(p.value))
modU_Annot %>% filter(sigPair == TRUE) %>% group_by(GTP) %>% 
  summarize(m = median(p.value))
modU_Annot %>% filter(sigPair == TRUE) %>% group_by(adults) %>% 
  summarize(m = median(p.value))


modU_comp$MESA <- ifelse(paste(modU_comp$CpG, modU_comp$TC) %in% monocytes.merge$pair, 
                         "p-value < 1e-5", "p-value > 1e-5")
modU_comp$GTP <- ifelse(paste(modU_comp$CpG, modU_comp$TC) %in% blood.merge$pair, 
                         "p-value < 1e-5", "p-value > 1e-5")
modU_comp$adult <- ifelse(modU_comp$sigPair, 
                           ifelse(modU_comp$MESA == "p-value < 1e-5" | modU_comp$GTP == "p-value < 1e-5", "Shared", "Children"), 
                           "Other")


wilcox.test(abs(subset(modU_comp, sigPair & adult == "Shared")$FC)/10,
        abs(subset(modU_comp, sigPair & adult == "Children")$FC)/10)

modU_Annot %>% filter(sigPair == TRUE) %>% group_by(adult) %>% 
  summarize(m = median(abs(FC)/10))

modU_comp %>% filter(sigPair == TRUE) %>% group_by(adult) %>% 
  summarize(m = median(SD))


modU_comp %>% filter(sigPair == TRUE) %>% ggplot(aes(x = abs(FC), color = adult)) +
  geom_density() + scale_x_continuous(limits = c(0, 6))
modU_comp %>% filter(sigPair == TRUE) %>% ggplot(aes(x = SD, color = adult)) +
  geom_density() + scale_x_continuous(limits = c(0, 3))
wilcox.test(abs(subset(modU_comp, sigPair & adult == "Shared")$SD),
        abs(subset(modU_comp, sigPair & adult == "Children")$SD))


wilcox.test(-log10(subset(modU_comp, sigPair & adult == "Shared")$p.value),
            -log10(subset(modU_comp, sigPair & adult == "Children")$p.value))


plot_all <- modU_comp %>% filter(sigPair == TRUE) %>% 
  mutate(Type = ifelse(adult == "Children", "Child-specific", "Age-shared")) %>%
  ggplot(aes(x = Distance, color = Type)) +
  geom_density() +
  theme_bw() + 
  scale_x_continuous(name = "CpG-TC TSS distance",
                     breaks = c(-5e5, -2e5, 0, 2e5, 5e5), 
                     labels = c("-500Kb", "-250Kb", "0", "250Kb", "500Kb")) +
  scale_y_continuous(name = "Density") + 
  scale_color_discrete(name = "eQTM type") +
  theme(plot.title = element_text(hjust = 0.5))
png("paper/distance_distr_Kennedy.png", width = 2000, height = 1000, res = 300)
plot_all
dev.off()

modU_comp %>% filter(sigPair == TRUE) %>% 
  mutate(Type = ifelse(adult == "Children", "Children-specific", "Children and adult shared")) %>%
  group_by(Type) %>%
  summarize(m = median(Distance),
            q25 = quantile(Distance, 0.25),
            q75 = quantile(Distance, 0.75))
ks.test(subset(modU_comp, sigPair == TRUE & adult == "Children")$Distance,
        subset(modU_comp, sigPair == TRUE & adult != "Children")$Distance) 

### Plot overlap between datasets
adult_db <- full_join(select(blood.merge.uniq, exp.Probe, CpG.probe, sigType2),
                       select(monocytes.merge.uniq, exp.Probe, CpG.probe, sigType2),
                       by = c("exp.Probe", "CpG.probe")) %>%
  mutate(name = paste(exp.Probe, CpG.probe),
         sigType2.x = ifelse(is.na(sigType2.x), "None", sigType2.x),
         sigType2.y = ifelse(is.na(sigType2.y), "None", sigType2.y),
         GTP = sigType2.x != "None",
         MESA = sigType2.y != "None",
         HELIX = sigType2.x == "Shared" | sigType2.y == "Shared")
adult_over <- plot(euler(adult_db[, c("GTP", "MESA", "HELIX")], shape = "ellipse"), 
           fills = list(fill = c("cyan", "darksalmon", "darkseagreen"), alpha = 0.5),
     quantities = list(fontsize = 8), 
     main = "Adult eQTMs in HELIX")



helix_over <- modU_Annot %>% 
  filter(sigPair == TRUE) %>%
  mutate(HELIX = TRUE,
         GTP = GTP == "GTP",
         MESA = MESA == "MESA") %>%
  select(GTP, MESA, HELIX) %>%
  euler(shape = "ellipse") %>%
  plot(fills = list(fill = c("cyan", "darksalmon", "darkseagreen"), alpha = 0.5),
       quantities = list(fontsize = 8), 
       main = "Helix eQTMs in adult cohorts")

bottom <- plot_grid(adult_over, helix_over, labels = c("B", "C"))

png("paper/AgeEffectPlot.png", width = 3000, height = 2000, res = 300)
plot_grid(age_var_comb, bottom, ncol = 1, labels = c("A", ""))
dev.off()

## Compare with cell types  ####
#### Prepare data in previous sections
childCpGs <- unique(modU_Annot[modU_Annot$adults == "Children", ]$CpG)
sharedCpGs <- unique(modU_Annot[modU_Annot$adults == "Shared", ]$CpG)

length(intersect(childCpGs, sharedCpGs))
length(setdiff(sharedCpGs, childCpGs))
length(setdiff(childCpGs, sharedCpGs))

childCpGs.f <- setdiff(childCpGs, sharedCpGs)
sharedCpGs.f <- setdiff(sharedCpGs, childCpGs)

data.frame(Fstat = FlowSorted.Blood.450k.compTable[c(childCpGs, sharedCpGs), "Fstat"],
           Type = rep(c("Children", "Shared"), c(length(childCpGs), length(sharedCpGs)))) %>%
  lm(formula = log10(Fstat) ~ Type, data = .) %>%
  summary()

data.frame(Fstat = FlowSorted.Blood.450k.compTable[c(childCpGs.f, sharedCpGs), "Fstat"],
           Type = rep(c("Children", "Shared"), c(length(childCpGs.f), length(sharedCpGs)))) %>%
  lm(formula = log10(Fstat) ~ Type, data = .) %>%
  summary()


methCell_child_adult <- data.frame(Fstat = FlowSorted.Blood.450k.compTable[c(childCpGs, sharedCpGs), "Fstat"],
                                  Type = rep(c("CpGs in child-specific eQTMs", "CpGs in age-shared eQTMs"), c(length(childCpGs), length(sharedCpGs)))) %>%
  ggplot(aes(y = log10(Fstat), x = Type, fill = Type)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(name = "") +
  ggtitle("Methylation cell-type specificity") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

## Select genes
childTCs <- unique(modU_Annot[modU_Annot$adults == "Children", ]$TC)
childGenes <- filter(expAnnot, probeset_id %in% childTCs)$GeneSymbol_Affy
childGenes <- unique(unlist(strsplit(childGenes, ";")))

sharedTCs <- unique(modU_Annot[modU_Annot$adults == "Shared", ]$TC)
sharedGenes <- filter(expAnnot, probeset_id %in% sharedTCs)$GeneSymbol_Affy
sharedGenes <- unique(unlist(strsplit(sharedGenes, ";")))

length(intersect(sharedTCs, childTCs))
length(setdiff(sharedTCs, childTCs))
length(setdiff(childTCs, sharedTCs))

## Remove genes regulated by different CpGs in both models
childGenes.f <- setdiff(childGenes, sharedGenes)
sharedGenes.f <- setdiff(sharedGenes, childGenes)

### Make results table
bluep_res_adult <- topTable(fit, n = Inf) %>% 
  mutate(gene = rownames(.),
         cat = ifelse(gene %in% childGenes, "TCs in child-specific eQTMs", 
                      ifelse(gene %in% sharedGenes, "TCs in age-shared eQTMs", "None"))) %>%
  filter(cat != "None")

summary(lm(formula = log10(F) ~ cat, data = bluep_res_adult))

exprCell_adults <- bluep_res_adult %>%
  ggplot(aes(y = log10(F), x = cat, fill = cat)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(name = "") +
  scale_y_continuous("log10(Fstat)") +
  ggtitle("Gene expression cell-type specificity") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")



png("paper/CompModelsChildrenAdult.png", width = 1500, height = 1500, res = 300)
plot_grid(methCell_child_adult, exprCell_adults, ncol = 1, labels = "AUTO")
dev.off()


## Compare with age variability ####
### Epidelta
epideltadf <- read.table("data/epidelta_2020-07-17.txt", as.is = TRUE)
epideltadf$CpG <- rownames(epideltadf)
ageVarCpG <- subset(epideltadf, M1.change.p < 1e-7)$CpG
agevar_adult <-  data.frame(CpG = c(childCpGs, sharedCpGs), 
                            eQTM = rep(c("Child-specific", "Child and adult shared"), 
                                       c(length(childCpGs), length(sharedCpGs))),
                            stringsAsFactors = FALSE) %>%
  rbind(data.frame(CpG = filter(methyAnnot, !Name %in% .$CpG)$Name, 
                   eQTM = "None")) %>%
  mutate(Epidelta = ifelse(CpG %in% ageVarCpG, "Variable", "Constant")) %>%
  group_by(eQTM) %>%
  summarize(ageVarin = sum(Epidelta == "Variable"), ageVarout = sum(Epidelta == "Constant"))

sharedAgevar <- getOR(2:3, agevar_adult[-2, ])
childAgevar <- getOR(2:3, agevar_adult[-1, ])
  
### MeDALL
agedf <- read.csv2("data/DiffMethyAgeCpGs.csv", as.is = TRUE)
agedf <- agedf %>%
  as_tibble() %>%
  mutate(Dir = ifelse(as.numeric(beta8.avg) > as.numeric(beta0.avg), "Increased", "Decreased"), 
         CpG = ILMNID)

agevar_adult2 <-  data.frame(CpG = c(childCpGs, sharedCpGs), 
                            eQTM = rep(c("Child-specific", "Child and adult shared"), 
                                       c(length(childCpGs), length(sharedCpGs))),
                            stringsAsFactors = FALSE) %>%
  rbind(data.frame(CpG = filter(methyAnnot, !Name %in% .$CpG)$Name, 
                   eQTM = "None")) %>%
  mutate(Medall = ifelse(CpG %in% agedf$CpG, "Variable", "Constant")) %>%
  group_by(eQTM) %>%
  summarize(ageVarin = sum(Medall == "Variable"), ageVarout = sum(Medall == "Constant"))

sharedAgevar2 <- getOR(2:3, agevar_adult2[-2, ])
childAgevar2 <- getOR(2:3, agevar_adult2[-1, ])


## Compare with meQTLs ####
load("results/eQTLanalysis/comQTLs.Rdata") 
mQTLs <- read.table("data/ARIES_mQTLs.tab", header = TRUE, as.is = TRUE)
mQTLsH2 <- read.table("results/ARIES/mqtls.txt", header = TRUE, as.is = TRUE)

ARIESannot <- read.table("data/ariesmqtlsnps.bim", as.is = TRUE)
colnames(ARIESannot) <- c("chr", "SNP", "cm", "pos", "Ref", "Alt")

HELIXannot <- read.table("~/data/WS_HELIX/HELIX_preproc/gwas/Final_data_HRCimp_QC2/HELIX.impQC.rs.bim", as.is = TRUE)
colnames(HELIXannot) <- c("chr", "SNP", "cm", "pos", "Ref", "Alt")

comMQTLs <- mQTLsH2 %>%
  mutate(gene = CpG) %>%
  dplyr::select(SNP, gene, A1, A2, freq, b, se, p, N, r2) %>%
  semi_join(rbind(comCisQTL, comTransQTL), by = c("SNP", "gene")) %>%
  inner_join(left_join(mQTLs, dplyr::select(ARIESannot, SNP, Ref, Alt), by = "SNP"), 
             by = c("SNP", "gene")) %>%
  as_tibble()

comMQTLs.f <- comMQTLs %>%
  filter(!is.na(Ref)) %>%
  filter(!(sign(b) != sign(beta) & A1 == Ref)) %>%
  filter(!(sign(b) == sign(beta) & A1 == Alt))

meqtlCpGs <- unique(comMQTLs.f$gene)

meQTL_Sum <- data.frame(CpG = c(childCpGs, sharedCpGs), 
                        Type = rep(c("Children-specific", "Children and adult shared"), 
                                   c(length(childCpGs), length(sharedCpGs)))) %>%
  rbind(data.frame(CpG = filter(methyAnnot, !Name %in% .$CpG)$Name, 
                   Type = "None")) %>%
  mutate(meQTL = CpG %in% meqtlCpGs) %>%
  group_by(Type) %>%
  summarize(meQTLin = sum(meQTL), no_meQTL = sum(!meQTL))

sharedmeQTL <- getOR(2:3, meQTL_Sum[-2, ])
childmeQTL <- getOR(2:3, meQTL_Sum[-1, ])

## Compare ROADMAP chromatin states ####
chromSt_adult <-  data.frame(CpG = c(childCpGs, sharedCpGs), 
                             eQTM = rep(c("CpGs in child-specific eQTMs", "CpGs in age-shared eQTMs"), 
                                        c(length(childCpGs), length(sharedCpGs)))) %>%
  rbind(data.frame(CpG = filter(methyAnnot, !Name %in% .$CpG)$Name, 
                   eQTM = "None")) %>%
  left_join(mutate(methyAnnot, CpG = Name) %>%  select(CpG, eval(chromStates))) %>%
  group_by(eQTM) %>%
  summarize_at(chromStates, list(sum = sum, sum2 = sum2, prop = mean))


chromSt_enrich_adult <- lapply(c("CpGs in child-specific eQTMs", "CpGs in age-shared eQTMs"), function(t){
  rbind(filter(chromSt_adult, eQTM == t),
        filter(chromSt_adult, eQTM == "None")) %>%
    select(ends_with("sum"), ends_with("sum2")) %>%
    g(chromStates, cols = c("_sum", "_sum2")) %>%
    as_tibble() %>%
    mutate(par = c("OR", "p.val", "ORlow", "ORhigh")) %>%
    gather("Region", "Value", 1:15) %>%
    mutate(eQTM = t)
})
chromSt_enrich_adult <- Reduce(rbind, chromSt_enrich_adult)

chrom_adult <- chromSt_enrich_adult %>%
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
  ),
  ) %>%
  ggplot(aes(x = Region, y = OR, fill = eQTM)) + 
  geom_bar(stat = "identity", position=position_dodge(), color= "black") + 
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin = ORlow, ymax = ORhigh)) +
  scale_y_continuous(trans = "log2", 
                     breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  scale_fill_discrete(name = "CpG type") +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = "ROADMAP chromatin states") +
  facet_wrap(~ Group, scales = "free_x") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

png("paper/enrich_chromStates_adult.png", width = 2500, height = 1500, res = 300)
chrom_adult
dev.off()

chromSt_adult_prop <- chromSt_adult %>%
  select(eQTM, ends_with("prop")) %>%
  gather(categories, proportion, 2:(2+length(chromStates) - 1)) %>%
  ungroup() %>%
  mutate(Type = as.character(eQTM),
         Type = ifelse(Type == "None", "CpGs not in eQTMs", Type),
         Type = factor(Type, levels = c("CpGs not in eQTMs", "CpGs in age-shared eQTMs", "CpGs in child-specific eQTMs")),
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
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  scale_x_discrete(name = "ROADMAP chromatin states") +
  facet_wrap(~ Group, scales = "free_x") +  
  scale_y_continuous(name = "Proportion of CpGs (%)") +
  scale_fill_manual(name = "CpG type", values = c("#FFFFFF", "#F8766D", "#00BFC4")) +
  theme_bw()

png("paper/chromStatesProp_age.png", width = 2500, height = 2000, res = 300)
chromSt_adult_prop
dev.off()
