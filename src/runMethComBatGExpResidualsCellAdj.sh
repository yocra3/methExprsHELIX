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
Rscript src/createInputFiles.R $meth $gexp $resFolder autosomes

## Divide data for analysis
Rscript src/divide_data.R $methInput $gexpInput $resFolder

## Run linear models
for i in {1..22}
do
  echo $i
  Rscript src/runLinearModelSubset.R $resFolder chr$i cell $resFolder
done