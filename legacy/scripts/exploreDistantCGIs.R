###############################################################################
# Test association with distal promoter (Did not work - not included in final version)
###############################################################################

library(dplyr)
library(ggplot2)
library(minfi)

load("results/preprocessFiles/Methylation_GRSet.RData")
load("results/preprocessFiles/Expression_SE_residuals.RData")
load("results/preprocessFiles/gexpAnnotation.Rdata")

load("results/MethComBatExpResidualsCellAdj/allres_simP_cpgs.Rdata")
load("results/preprocessFiles/methyAnnotation.Rdata")
load("results/preprocessFiles/allOverlaps.Rdata")

rownames(methyAnnot) <- methyAnnot$Name
gset_isles <- gset[rowData(gset)$HMM_Island != "", ]

df_sig <- df[df$sigPair, ]
df_comp <- inner_join(df_sig, overDF, by = c("CpG", "TC"))
df_comp_sel <- df_comp %>%
  as_tibble() %>%
  filter(Distance > -50e3 & Distance < -1.5e3 & CpG %in% rownames(gset_isles))


expAnnDF <- as_tibble(as.data.frame(expAnnot))
expAnnDF$TC <- expAnnDF$transcript_cluster_id

df_comp_sel <- df_comp_sel %>%
  select(CpG, TC, FC, p.value, Distance) %>%
  left_join(select(expAnnDF, TC, GeneSymbol_Affy), by = "TC") %>%
  mutate(GeneAffy = strsplit(GeneSymbol_Affy, ";"))

methGenes <- unique(unlist(rowData(gset_isles)$UCSC_RefGene_Name))
expGenes <- unique(unlist(df_comp_sel$GeneAffy))
comGenes <- intersect(methGenes, expGenes)

gset_end <- gset_isles[grepl(paste(comGenes, collapse = "|"), methGenes),]
gset_prom <- gset_end[grepl("TSS", rowData(gset_end)$UCSC_RefGene_Group),]

hist(as.vector(getBeta(gset_prom)))