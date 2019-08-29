###############################################################################
#' Create table with features descriptives
###############################################################################

library(dplyr)
library(S4Vectors)
library(minfi)
library(ggplot2)
library(cowplot)

## Load data  ####
load("results/preprocessFiles/gexpAnnotation.Rdata")
load("results/preprocessFiles/methyAnnotation.Rdata")
load("results/preprocessFiles/allOverlaps.Rdata")
load("results/preprocessFiles/Methylation_GRSet.RData")

## Total number CpGs
nrow(methyAnnot)
#[1] 396233

## Sexual CpGs
sum(methyAnnot$chr %in% c("chrX", "chrY"))
#[1] 9715

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
# [1] 252449


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

