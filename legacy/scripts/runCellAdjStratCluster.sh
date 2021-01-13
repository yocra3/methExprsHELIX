#!/bin/sh

###############################################################################
# Methylation vs Expression Analysis using:
#   - Methylation: ComBat
#   - Expression: SVA Residuals protecting age, sex, cohort and cell counts
#   - Model: adjust for age, sex, cohort, and cell counts
#   - Stratified males and females
#   - ALl chromosomes
###############################################################################

#set the job name
#SBATCH --job-name=cellAdjStrat

#set the number of CPUS per task
#SBATCH --ntasks-per-node=1

#set the memory
#SBATCH --mem=15000

# job output file information
#SBATCH -o cellAdjStratPrepare.out

# job errors file
#SBATCH -e cellAdjStratPrepare.err

# set the partition where the job will run
#SBATCH --partition=normal
 
# set the number of nodes
#SBATCH --nodes=1
 
# set max wallclock time
#SBATCH --time=1:00:00
 
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=END
 
# send mail to this address
#SBATCH --mail-user= carlos.ruiz@isglobal.org

module load R/3.5.1-foss-2018b

## Define paths to omic files
meth="results/preprocessFiles/Methylation_GRSet.RData"
gexp="results/preprocessFiles/Expression_SE_residuals.RData"
resFolder="results/MethComBatExpResidualsCellAdjStrat"

methMaleInput="${resFolder}/methyInputMale.Rdata"
methFemaleInput="${resFolder}/methyInputFemale.Rdata"
gexpMaleInput="${resFolder}/gexpInputMale.Rdata"
gexpFemaleInput="${resFolder}/gexpInputFemale.Rdata"

## Create folder for result
mkdir $resFolder

## Create folder for result per sex
mkdir $resFolder/male
mkdir $resFolder/female


## Generate Input data for analysis
R CMD BATCH --vanilla '--args methy="'$meth'" exprs="'$gexp'" out_fold="'$resFolder'" sex-stratify' src/createInputFiles.R  $resFolder/input.out

## Divide data for analysis
### Male
R CMD BATCH --vanilla '--args methy="'$methMaleInput'" exprs="'$gexpMaleInput'" out_fold="'$resFolder/male'"' src/divide_data.R  $resFolder/divideMale.out

### Female
R CMD BATCH --vanilla '--args methy="'$methFemaleInput'" exprs="'$gexpFemaleInput'" out_fold="'$resFolder/female'"' src/divide_data.R  $resFolder/divideFemale.out