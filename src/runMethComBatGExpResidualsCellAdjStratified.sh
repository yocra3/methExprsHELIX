#!/bin/sh

###############################################################################
# Methylation vs Expression Analysis using:
#   - Methylation: ComBat
#   - Expression: SVA Residuals protecting age, sex, cohort and cell counts
#   - Model: adjust for age, sex, cohort, cell counts
#   - Stratified
#   - All chromosomes
###############################################################################

## Define paths to omic files
meth="results/preprocessFiles/Methylation_GRSet.RData"
gexp="results/preprocessFiles/Expression_SE_residuals.RData"
resFolder="results/MethComBatExpResidualsCellAdjStrat"

methMaleInput="${resFolder}/methyInputMale.Rdata"
methFemaleInput="${resFolder}/methyInputFemale.Rdata"
gexpMaleInput="${resFolder}/gexpInputMale.Rdata"
gexpFemaleInput="${resFolder}/gexpInputFemale.Rdata"

# Male 
maleFolder="$resFolder/male"
femaleFolder="$resFolder/female"

## Create folder for result
mkdir $resFolder

## Create folder for result per sex
mkdir $maleFolder
mkdir $femaleFolder

## Generate Input data for analysis
Rscript src/createInputFiles.R methy="'$meth'" exprs="'$gexp'" out_fold="'$resFolder'" "sex-stratify"

## Divide data for analysis
### male
Rscript src/divide_data.R methy="'$methMaleInput'" exprs="'$gexpMaleInput'" out_fold="'$maleFolder'"
  
### female
Rscript src/divide_data.R methy="'$methFemaleInput'" exprs="'$gexpFemaleInput'" out_fold="'$femaleFolder'"

## Run linear models
### Male
for i in {1..22}
do
  echo $i
  Rscript src/runLinearModelSubset.R data_fold="'$maleFolder'" chr="'chr$i'" model="'cellStrat'" out_fold="'$maleFolder'"
done
Rscript src/runLinearModelSubset.R data_fold="'$maleFolder'" chr="'chrX'" model="'cellStrat'" out_fold="'$maleFolder'"
Rscript src/runLinearModelSubset.R data_fold="'$maleFolder'" chr="'chrY'" model="'cellStrat'" out_fold="'$maleFolder'"

### Female
for i in {1..22}
do
  echo $i
  Rscript src/runLinearModelSubset.R data_fold="'$femaleFolder'" chr="'chr$i'" model="'cellStrat'" out_fold="'$femaleFolder'"
done
Rscript src/runLinearModelSubset.R data_fold="'$femaleFolder'" chr="'chrX'" model="'cellStrat'" out_fold="'$femaleFolder'"


## Run simulations in server
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

