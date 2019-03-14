#!/bin/sh

###############################################################################
# Methylation vs Expression Analysis using:
#   - Methylation: ComBat
#   - Expression: SVA Residuals protecting age, sex, cohort and cell counts
#   - Model: adjust for age, sex, cohort, cell counts
#   - Stratified
#   - All chromosomes
###############################################################################

## Run lm models and simulations in server
resFolder="results/MethComBatExpResidualsCellAdjStrat"

# Male 
maleFolder="$resFolder/male"

## Generate beta dsitributions per CpG
Rscript src/getBetaDistributionsCpGs.R folder="'$maleFolder'" type="'male'"

## Generate beta dsitributions per gene
Rscript src/getBetaDistributionsGenes.R folder="'$maleFolder'" type="'male'"

# Female 
femaleFolder="$resFolder/female"

## Generate beta dsitributions per CpG
Rscript src/getBetaDistributionsCpGs.R folder="'$femaleFolder'" type="'female'"

## Generate beta dsitributions per gene
Rscript src/getBetaDistributionsGenes.R folder="'$femaleFolder'" type="'female'"

