#!/bin/sh

###############################################################################
# Methylation vs Expression Analysis using:
#   - Methylation: ComBat
#   - Expression: SVA Residuals protecting age, sex, cohort and cell counts
#   - Model: adjust for age, sex, cohort
#   - Stratified
#   - All chromosomes
###############################################################################

## Run lm models and simulations in server
resFolder="results/MethComBatExpResidualsNoCellAdjStrat"

# Male 
maleFolder="$resFolder/male"

## Generate beta dsitributions per CpG
Rscript src/getBetaDistributionsCpGs.R folder="'$maleFolder'" type="'male'"
Rscript src/getSignificantPairs.R resFolder="'$maleFolder'" base="'cpgs'" distribution="'$maleFolder/CpGsDistr.Rdata'"

## Generate beta dsitributions per gene
Rscript src/getBetaDistributionsGenes.R folder="'$maleFolder'" type="'male'"

# Female 
femaleFolder="$resFolder/female"

## Generate beta dsitributions per CpG
Rscript src/getBetaDistributionsCpGs.R folder="'$femaleFolder'" type="'female'"
Rscript src/getSignificantPairs.R resFolder="'$femaleFolder'" base="'cpgs'" distribution="'$femaleFolder/CpGsDistr.Rdata'"

## Generate beta dsitributions per gene
Rscript src/getBetaDistributionsGenes.R folder="'$femaleFolder'" type="'female'"

