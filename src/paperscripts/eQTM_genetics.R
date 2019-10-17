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

## Common mQTLs between ARIES and HELIX ####
load("results/eQTLanalysis/comQTLs.Rdata")

comCpGs <- unique(c(as.character(comCisQTL$gene), as.character(comTransQTL$gene)))

CpGsSum %>%
  mutate(mQTL = CpG %in% comCpGs,
         cisQTL = CpG %in% comCisQTL$gene,
         transQTL = CpG %in% comTransQTL$gene) %>%
  group_by(Combined) %>%
  summarize_if(is.logical, list(sum, mean))
# Combined     mQTL_fn1 cisQTL_fn1 transQTL_fn1 mQTL_fn2 cisQTL_fn2 transQTL_fn2
# <chr>           <int>      <int>        <int>    <dbl>      <dbl>        <dbl>
#  Mono_Inverse     3424       3394           77   0.271      0.269       0.00609
#  Mono_Positi…     2474       2438           71   0.276      0.272       0.00792
#  Multi_Both       1362       1360           41   0.344      0.343       0.0103
#  Multi_Inver…     2061       2037           68   0.336      0.332       0.0111
#  Multi_Posit…     1268       1259           28   0.359      0.357       0.00793
#  Non-signifi…    26129      24471         1944   0.0744     0.0697      0.00554

CpGsSum %>%
  mutate(mQTL = CpG %in% comCpGs,
         cisQTL = CpG %in% comCisQTL$gene,
         transQTL = CpG %in% comTransQTL$gene,
         Type = ifelse(Combined != "Non-significant", "Significant", "Non-significant")) %>%
  group_by(Type) %>%
  summarize_if(is.logical, list(sum, mean))
# Type         mQTL_fn1 cisQTL_fn1 transQTL_fn1 mQTL_fn2 cisQTL_fn2 transQTL_fn2
# <chr>           <int>      <int>        <int>    <dbl>      <dbl>        <dbl>
# Non-signifi…    26129      24471         1944   0.0744     0.0697      0.00554
# Significant     10589      10488          285   0.301      0.298       0.00809

