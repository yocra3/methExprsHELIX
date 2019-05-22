###############################################################################
# Code for used plots
###############################################################################

# Models comparison
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(UpSetR)

load("results/preprocessFiles/gexpAnnotation.Rdata")

load("results/MethComBatExpResidualsNoCellAdj/allres_simP_cpgs.Rdata")
modU <- df
featsU <- featStatsDF

load("results/MethComBatExpResidualsCellAdj/allres_simP_cpgs.Rdata")
modC <- df
featsC <- featStatsDF

## CpGs identified
sigFU <- subset(featsU, p.val.adj < 0.05)$feat
sigFC <- subset(featsC, p.val.adj < 0.05)$feat
sigCom <- intersect(sigFU, sigFC)

png("compareModelsCpG.png")
upset(fromList(list("Model A" = sigFU, "Model B" = sigFC)), order.by = "freq", 
      set_size.angles = 45, text.scale = 3.5)
dev.off()

# Multiple testing comparison
library(dplyr)
library(UpSetR)

load("results/preprocessFiles/gexpAnnotation.Rdata")
expAnnot$TC <- expAnnot$transcript_cluster_id

## No cell
load("results/MethComBatExpResidualsNoCellAdj/allres_simP_cpgs.Rdata")
modU <- df %>%
  as_tibble() %>%
  mutate(adj.p.value = p.adjust(p.value, "BH"),
         sigPair_BH = adj.p.value < 0.05,
         adj.p.value.BF = p.adjust(p.value, "bonferroni"),
         sigPair_BF = adj.p.value.BF < 0.05,
         pair = paste(CpG, TC))
featsU <- featStatsDF

## CpGs identified
permsCp <- unique(modU[modU$sigPair, ]$CpG)
BHCp <- unique(modU[modU$sigPair_BH, ]$CpG)
BFCp <- unique(modU[modU$sigPair_BF, ]$CpG)

png("CpGsComparisonPval.png", height = 600)
upset(fromList(list(Permutation = permsCp, BH = BHCp, Bonferroni = BFCp)), 
      order.by = "freq", set_size.angles = 45, text.scale = 3.5)
dev.off()