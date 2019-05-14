#!/bin/sh

###############################################################################
# Methylation vs Expression Analysis using:
#   - Methylation: ComBat
#   - Expression: SVA Residuals protecting age, sex, cohort and cell counts
#   - Model: adjust for age, sex, cohort, cell counts
#   - All Samples
#   - Autosome chromosomes
###############################################################################

## Define paths to omic files
meth="results/preprocessFiles/Methylation_GRSet.RData"
gexp="results/preprocessFiles/Expression_SE_residuals.RData"
resFolder="results/MethComBatExpResidualsCellAdj"

methInput="${resFolder}/methyInput.Rdata"
gexpInput="${resFolder}/gexpInput.Rdata"

## Create folder for result
mkdir $resFolder

## Generate Input data for analysis
Rscript src/createInputFiles.R methy="'$meth'" exprs="'$gexp'" out_fold="'$resFolder'" "autosomes"

## Divide data for analysis
Rscript src/divide_data.R methy="'$methInput'" exprs="'$gexpInput'" out_fold="'$resFolder'"
  
## Run linear models
for i in {1..22}
do
  echo $i
  Rscript src/runLinearModelSubset.R data_fold="'$resFolder'" chr="'chr$i'" model="'cell'" out_fold="'$resFolder'"
done


## Generate beta dsitributions per CpG
Rscript src/getBetaDistributionsCpGs.R folder="'$resFolder'" type="'autosome'"
Rscript src/getSignificantPairs.R resFolder="'$resFolder'" base="'cpgs'" distribution="'$resFolder/CpGsDistr.Rdata'"

## Generate beta dsitributions per gene
Rscript src/getBetaDistributionsGenes.R folder="'$resFolder'" type="'autosome'"
