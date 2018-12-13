###############################################################################
#' Create SummarizedExperiment for analysis (Gene Expression)
###############################################################################

library(Biobase)
library(SummarizedExperiment)

## Load SummarizedExperiment
load("./data/transcriptome_subcohort_notfitr_inclsex_v3.RData")
gexp <- transcriptome_subcohort_notfitr_inclsex

## Load common IDs
load("results/preprocessFiles/comIds.Rdata")

# Select common IDs
gexp <- gexp[ , comIds]

## Remove probes with low call rate
load("./data/badProbes.Rdata")
gexp <- gexp[fData(gexp)$fil1 == "no", ]

# Making SE
se <- makeSummarizedExperimentFromExpressionSet(gexp)
save(se, file = "results/preprocessFiles/Expression_SE_raw.RData")

