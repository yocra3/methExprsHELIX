#'##############################################################################
#'##############################################################################
#' Add eQTMs to CpGs found in HELIX
#'##############################################################################
#'##############################################################################

## Load libraries ####
library(dplyr)
library(tidyr)

## Load datasets
load("results/MethComBatExpResidualsNoCellAdj/allres_simP_cpgs.Rdata")
cpgList <- read.delim("data/list_cpgs_helix.txt", as.is = TRUE)
load("results/preprocessFiles/gexpAnnotation.Rdata")

sigDf <- df %>%
  as_tibble() %>%
  filter(sigPair)

HELIXeQTMs <- cpgList %>%
  as_tibble() %>%
  mutate(CpG = Molecular.feature) %>%
  distinct() %>%
  left_join(sigDf, by = "CpG") %>%
  select(-Molecular.feature, -sigPair) %>%
  left_join(select(expAnnot, TC, GeneSymbol_Affy))

write.table(HELIXeQTMs, file = "paper/HELIXeQTMs_list.tab", 
            quote = FALSE, row.names = FALSE)