#!/bin/sh

###############################################################################
# Methylation vs Expression Analysis using data from previous project
#   - Model: adjust for age, sex, cohort
#   - All Samples
#   - Autosome chromosomes
###############################################################################

## Define paths to omic files
methInput="archive/results/data/gset_sm_1.RData"
gexpInput="archive/results/data/summ_exp_sm_res2.RData"
resFolder="results/SergiObjectsNoCellAdj"

## Create folder for result
mkdir $resFolder

## Change name to object
echo 'load("'$gexpInput'")
se <- se_res
save(se, file = "'$resFolder/gexp.Rdata'")
' | Rscript --vanilla -

## Divide data for analysis
Rscript src/divide_data.R methy="'$methInput'" exprs="'$resFolder/gexp.Rdata'" out_fold="'$resFolder'"
  
## Run linear models
for i in {1..22}
do
  echo $i
  Rscript src/runLinearModelSubset.R data_fold="'$resFolder'" chr="'chr$i'" model="'nocell'" out_fold="'$resFolder'"
done

## Generate beta dsitributions per CpG
Rscript src/getBetaDistributionsCpGs.R folder="'$resFolder'" type="'autosome'"
Rscript src/getSignificantPairs.R resFolder="'$resFolder'" base="'cpgs'" distribution="'$resFolder/CpGsDistr.Rdata'"

## Generate beta dsitributions per gene
Rscript src/getBetaDistributionsGenes.R folder="'$resFolder'" type="'autosome'"
