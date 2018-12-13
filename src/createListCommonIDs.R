#'##############################################################################
#' Create list of common IDs between methylation and gene Expression
#'##############################################################################

# Load libraries
library(minfi)
library(SummarizedExperiment)

# Load Methylation data ####
load("./data/methylome_subcohort_ComBatSlide_6cells_notfitr_v3.Rdata")
meth <- methylome_subcohort_ComBatSlide_6cells_notfitr

# Load Gene Expression data ####
load("data/transcriptome_subcohort_notfitr_inclsex_v3.RData")
gexp <- transcriptome_subcohort_notfitr_inclsex

# Filter Gene Expression dataset
## Select European and discard panel 1B samples
gexp_sel <- gexp[, gexp$Period != "1B" & gexp$h_ethnicity_cauc == "yes"]

## Select common samples between Exprs and Methylation
comIds <- intersect(colnames(gexp_sel), colnames(meth))
save(comIds, file = "results/preprocessFiles/comIds.Rdata")