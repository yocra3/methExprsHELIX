#'##############################################################################
#' Genetic variants and eQTMs
#' This script contains the code for
#'  - Heritability
#'  - Association with genetics
#'##############################################################################

# Load data and libraries ####
library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyr)
library(GenomicRanges)

load("results/MethComBatExpResidualsNoCellAdj/allres_simP_cpgs.Rdata")
load("results/preprocessFiles/allOverlaps.Rdata")
load("results/preprocessFiles/methyAnnotation.Rdata")
load("results/preprocessFiles/gexpAnnotation.Rdata")

# Get useful variables ####
CpGsSum <- df %>%
  group_by(CpG) %>%
  summarise(Type = ifelse(sum(sigPair) == 0, "Non-significant",
                          ifelse(sum(sigPair) == 1, "Mono", "Multi")),
            Direction = ifelse(sum(sigPair) == 0, "Non-significant",
                               ifelse(all(FC[sigPair] > 0), "Positive", 
                                      ifelse(all(FC[sigPair] < 0), "Inverse", "Both")))) %>%
  mutate(Combined = ifelse(Type == "Non-significant", 
                           "Non-significant", 
                           paste(Type, Direction, sep = "_")))



## Heritability ####
## Create summary tibble
h2df <- read.csv("data/HeritabilityCpGs", as.is = TRUE)
herm <- h2df %>%
  as_tibble() %>%
  ## Remove CpGs with NA in h2
  filter(!is.na(h2_total)) %>%
  mutate(CpG = cgid) %>%
  inner_join(CpGsSum, h2df, by = "CpG") %>%
  mutate(Combined = factor(Combined, 
                           levels = c("Mono_Inverse", "Mono_Positive", "Multi_Inverse", "Multi_Positive", "Multi_Both", "Non-significant")))
  

h2tot <- herm %>%
  ggplot(aes(x = Combined, y = h2_total, fill = Combined)) + 
  geom_boxplot() + 
  geom_hline(yintercept = c(0.2, 0.5), linetype="dashed", colour = "blue") + 
  ggtitle("Total Heritability") +
  scale_x_discrete(name = "CpG Type") +
  scale_y_continuous(name = "h2", limits = c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_fill_manual(name = "CpG Type", values = c("#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73", "#FFFFFF")) 


h2SNP <- herm %>%
  ggplot(aes(x = Combined, y = h2_SNPs, fill = Combined)) + 
  geom_boxplot() + 
  geom_hline(yintercept = c(0.2, 0.5), linetype="dashed", colour = "blue") + 
  ggtitle("SNP Heritability") +
  scale_x_discrete(name = "CpG Type") +
  scale_y_continuous(name = "h2", limits = c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_fill_manual(name = "CpG Type", values = c("#56B4E9", "#E69F00", "#0072B2", "#D55E00", "#009E73", "#FFFFFF")) 

png("paper/eQTMsHerit.png", width = 3500, height = 2000, res = 300)
plot_grid(h2tot, h2SNP, ncol = 2)
dev.off()

## ARIES mQTLs ####
mQTLs <- read.table("data/ARIES_mQTLs.tab", header = TRUE, as.is = TRUE)

ariesCpGs <- unique(mQTLs$gene)
CpGsSum %>%
  mutate(mQTL = CpG %in% ariesCpGs) %>%
  group_by(Combined) %>%
  summarize(n = sum(mQTL),
            prop = mean(mQTL))

CpGsSum %>%
  mutate(mQTL = CpG %in% ariesCpGs,
         Type = ifelse(Combined != "Non-significant", "Significant", "Non-significant")) %>%
  group_by(Type) %>%
  summarize(n = sum(mQTL),
            prop = mean(mQTL))