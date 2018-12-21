#!/bin/sh

###############################################################################
# Run linear models for all the different models
###############################################################################

#set the job name
#SBATCH --job-name=allModels

#set the number of CPUS per task
#SBATCH --ntasks-per-node=16

#set the memory
##SBATCH --mem=15000

# job output file information
#SBATCH -o allModels.out

# job errors file
#SBATCH -e allModels.err

# set the partition where the job will run
#SBATCH --partition=normal
 
# set the number of nodes
#SBATCH --nodes=1
 
# set max wallclock time
#SBATCH --time=4:00:00
 
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=END
 
# send mail to this address
#SBATCH --mail-user= carlos.ruiz@isglobal.org

module load R/3.5.1-foss-2018b

## Cell adjusted
resFolder="results/MethComBatExpResidualsCellAdj"
for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'"' src/runLinearModelSubset.R $resFolder/modchr$i.out
done


## No cell adjusted
resFolder="results/MethComBatExpResidualsNoCellAdj"
for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'"' src/runLinearModelSubset.R $resFolder/modchr$i.out
done

## Cell adjusted Male
resFolder="results/MethComBatExpResidualsCellAdjStrat/male"
model="cellStrat"
for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'"' src/runLinearModelSubset.R $resFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'"' src/runLinearModelSubset.R $resFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'"' src/runLinearModelSubset.R $resFolder/modchrY.out

## Cell adjusted Female
resFolder="results/MethComBatExpResidualsCellAdjStrat/female"
for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'"' src/runLinearModelSubset.R $resFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'"' src/runLinearModelSubset.R $resFolder/modchrX.out

## No cell adjusted Male
resFolder="results/MethComBatExpResidualsNoCellAdjStrat/male"
model="nocellStrat"
for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'"' src/runLinearModelSubset.R $resFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'"' src/runLinearModelSubset.R $resFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'"' src/runLinearModelSubset.R $resFolder/modchrY.out

## Cell adjusted Female
resFolder="results/MethComBatExpResidualsCellAdjStrat/female"
for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'"' src/runLinearModelSubset.R $resFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'"' src/runLinearModelSubset.R $resFolder/modchrX.out

