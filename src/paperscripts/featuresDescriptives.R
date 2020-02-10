###############################################################################
#' Create table with features descriptives
###############################################################################

library(dplyr)
library(S4Vectors)
library(minfi)
library(ggplot2)
library(cowplot)
library(matrixStats)

## Load data  ####
load("results/preprocessFiles/gexpAnnotation.Rdata")
load("results/preprocessFiles/methyAnnotation.Rdata")
load("results/preprocessFiles/allOverlaps.Rdata")
load("results/preprocessFiles/Methylation_GRSet.RData")
load("results/preprocessFiles/Expression_SE_residuals.RData")

## Total number CpGs
nrow(methyAnnot)
#[1] 396522

## Sexual CpGs
sum(methyAnnot$chr %in% c("chrX", "chrY"))
#[1] 10004

## Total number TCs
nrow(expAnnot)
#[1] 60692

## Sexual TCs
sum(expAnnot$seqname %in% c("chrX", "chrY"))
# [1] 2438


## Total number coding TCs
sum(expAnnot$Coding == "coding")
# [1] 24158

## Sexual coding TCs
sum(expAnnot$Coding == "coding" & expAnnot$seqname %in% c("chrX", "chrY"))
# [1] 1104

## Total number pairs
nrow(overDF)
#[1] 13868331

## Sexual pairs
length(grep("X|Y", overDF$TC))
# [1] 256200

## CpGs in pairs
length(unique(overDF$CpG))
# [1] 396414

## Sexual CpGs in pairs
sum(subset(methyAnnot, Row.names %in% unique(overDF$CpG))$chr %in% c("chrX", "chrY"))
#[1] 9996


## TCs in pairs
length(unique(overDF$TC))
# [1] 60363

## CpGs not in pairs
methyAnnot[!methyAnnot$Row.names %in% unique(overDF$CpG) , 1:10]

## TCs not in pairs
expAnnot[!expAnnot$probeset_id %in% overDF$TC, ]
           

## Coding TCs not in pairs
subset(expAnnot[!expAnnot$probeset_id %in% overDF$TC, ], Coding == "coding")


## Distribution CpGs per methylation level
prop.table(table(methyAnnot$median_cat))*100
# low   medium     high
# 42.05744 10.48441 47.45815

## Distribution CpGs per methylation variability
prop.table(table(methyAnnot$variability))*100
# invariant   variant
# 44.0306   55.9694
# 

## Proportion TCs with call rate == 100%
mean(expAnnot$CallRate > 90)
  # invariant   variant
# 44.0306   55.9694
# 


## Methylation distribution
betas <- as.vector(getBeta(gset))
png("paper/MethylationDistribution.png", width = 2000, height = 1500, res = 300)
betas %>%
  tibble::enframe(name = "a") %>%
  ggplot(aes(x = value)) + geom_density() + 
  theme_bw() + 
  scale_x_continuous("Methylation (betas)") +
  scale_y_continuous("", breaks = NULL) + 
  geom_vline(xintercept = c(0.3, 0.7), linetype = 2)
dev.off()

## CpGs/TC
CpG_plot <- overDF %>%
  group_by(CpG) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = n)) + geom_histogram(binwidth = 10) +
  scale_x_continuous("", limits = c(0, 1250)) +
  scale_y_continuous("")  +
  ggtitle("TCs paired with each CpG")

TC_plot <- overDF %>%
  group_by(TC) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = n)) + geom_histogram(binwidth = 10) +
  scale_x_continuous("", limits = c(0, 1250)) +
  scale_y_continuous("")  +
  ggtitle("CpGs paired with each TC")

png("paper/CpGs_TC_distr.png", width = 2000, height = 1500, res = 300)
plot_grid(CpG_plot, TC_plot, labels = "AUTO", nrow = 2)
dev.off()

## CpGs per TC
overDF %>%
  group_by(CpG) %>%
  summarize(n = n()) %>%
  summary()
# CpG               n
# cg00000029:     1   Min.   :  1.00
# cg00000108:     1   1st Qu.: 19.00
# cg00000109:     1   Median : 30.00
# cg00000165:     1   Mean   : 34.99
# cg00000236:     1   3rd Qu.: 46.00
# cg00000289:     1   Max.   :162.00
# (Other)   :396408

## TCs per CpG
overDF %>%
  group_by(TC) %>%
  summarize(n = n()) %>%
  summary()
# TC              n
# TC01000001.hg.1:    1   Min.   :   1.0
# TC01000002.hg.1:    1   1st Qu.:  89.0
# TC01000003.hg.1:    1   Median : 158.0
# TC01000004.hg.1:    1   Mean   : 229.8
# TC01000005.hg.1:    1   3rd Qu.: 291.0
# TC01000006.hg.1:    1   Max.   :3198.0
# (Other)        :60357

## Distribution pairs per chromosomes ####
png("paper/allPairs_Chr_distr.png", width = 2000, height = 1500, res = 300)
overDF %>% 
  mutate(chr = substring(TC, 3, 4),
         chr = gsub("^0", "", chr)) %>%
  group_by(chr) %>%
  summarize(n = n()) %>%
  mutate(chr = factor(chr, levels = c(1:22, "X", "Y"))) %>%
  ggplot(aes(x = chr, y = n)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(name = "Chromosome") +
  scale_y_continuous(name = "Number Pairs") +
  theme_bw()
dev.off()

## Distribution variability ####
png("paper/methy_Var.png", width = 2000, height = 1500, res = 300)
methyAnnot %>% 
  data.frame() %>%
  ggplot(aes(x = meth_range)) + 
  geom_histogram() +
  scale_x_continuous(name = "Methylation Range") +
  scale_y_continuous(name = "") +
  geom_vline(xintercept = 0.05, linetype = 2) +
  theme_bw()
dev.off()

## Correlation variability vs median methylation ####
png("paper/methy_Var_medianMeth.png", width = 2000, height = 1500, res = 300)
methyAnnot %>% 
  data.frame() %>%
  ggplot(aes(x = median_cat, y = meth_range)) + 
  geom_boxplot() +
  scale_x_discrete(name = "Median Methylation") +
  scale_y_continuous(name = "Methylation Range") +
  theme_bw()
dev.off()

summary(lm(meth_range ~ median_cat, methyAnnot))
# 
# Call:
#   lm(formula = meth_range ~ median_cat, data = methyAnnot)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max
# -0.17404 -0.04242 -0.02303  0.01670  0.86738
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)      0.0598260  0.0001875  319.00   <2e-16 ***
#   median_catmedium 0.1641175  0.0004198  390.91   <2e-16 ***
#   median_cathigh   0.0210212  0.0002576   81.61   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07659 on 396519 degrees of freedom
# Multiple R-squared:  0.2812,    Adjusted R-squared:  0.2812
# F-statistic: 7.757e+04 on 2 and 396519 DF,  p-value: < 2.2e-16

summary(lm(meth_range ~ factor(median_cat == "medium"), methyAnnot))

## Expression distribution
png("paper/ExpressionDistribution.png", width = 2000, height = 1500, res = 300)
se %>%
  assay() %>%
  rowMedians() %>%
  tibble::enframe(name = "a") %>%
  ggplot(aes(x = value)) + geom_density() + 
  theme_bw() + 
  scale_x_continuous("Median Expression per TC") +
  scale_y_continuous("", breaks = NULL) 
dev.off()


## Expression Call Rate distribution
png("paper/ExpressionCRDistribution.png", width = 2000, height = 1500, res = 300)
expAnnot %>%
  ggplot(aes(x = CallRate)) + 
  geom_histogram() + 
  theme_bw() + 
  scale_x_continuous("TCs Call Rate") +
  scale_y_continuous("", breaks = NULL) 
dev.off()

## Correlation median Expression and Call Rate
medExp <- rowMedians(assay(se))
names(medExp) <- rownames(se)

if (identical(expAnnot$probeset_id, names(medExp))) {
  png("paper/ExpressionMedianCRDistribution.png", width = 2000, height = 1500, res = 300)
  data.frame(median = medExp, CallRate = expAnnot$CallRate) %>%
    mutate(CallRateCat = ifelse(CallRate > 90, ">90%", "<90%")) %>%
    ggplot(aes(x = CallRateCat, y = median)) +
    geom_boxplot() +
    scale_x_discrete(name = "TC Call Rate") +
    scale_y_continuous(name = "Median Expression") +
    theme_bw()
  dev.off()
  
  data.frame(median = medExp, CallRate = expAnnot$CallRate) %>%
    mutate(CallRateCat = ifelse(CallRate > 90, ">90%", "<90%")) %>%
    lm(median ~ CallRateCat, .) %>%
    summary()
}
# Call:
#   lm(formula = median ~ CallRateCat, data = .)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max
# -3.5863 -0.6579 -0.2111  0.3318 15.8389
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)      3.63635    0.01003   362.5   <2e-16 ***
#   CallRateCat>90%  2.35810    0.01479   159.4   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.816 on 60690 degrees of freedom
# Multiple R-squared:  0.2951,    Adjusted R-squared:  0.2951
# F-statistic: 2.541e+04 on 1 and 60690 DF,  p-value: < 2.2e-16



