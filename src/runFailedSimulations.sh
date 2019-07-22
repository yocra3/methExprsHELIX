#!/bin/sh

###############################################################################
# Run linear models for all the different models
###############################################################################

#set the job name
#SBATCH --job-name=failedSims

#set the number of CPUS per task
#SBATCH --ntasks-per-node=16

#set the memory
#SBATCH --mem=35000

# job output file information
#SBATCH -o failedSims.out

# job errors file
#SBATCH -e failedSims.err

# set the partition where the job will run
#SBATCH --partition=normal
 
# set the number of nodes
#SBATCH --nodes=1
 
# set max wallclock time
#SBATCH --time=30:00:00
 
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=END
 
# send mail to this address
#SBATCH --mail-user= carlos.ruiz@isglobal.org

module load R/3.5.1-foss-2018b

## Cell adjusted
echo Cell Adjusted
resFolder="results/MethComBatExpResidualsCellAdjStrat"

## Male
maleFolder="$resFolder/male"

# sim1
num=1
simFolder="$maleFolder/sim$num"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$maleFolder'" chr="chr'$i'" model="cellStrat" out_fold="'$maleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$maleFolder'" chr="chrX" model="cellStrat" out_fold="'$maleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
mv $simFolder/*.* /gpfs42/projects/lab_helix_omics/shared_data/methExprsHELIX/$simFolder/

# sim49
num=49
simFolder="$maleFolder/sim$num"
mkdir $simFolder

R CMD BATCH '--args data_fold="'$maleFolder'" chr="chr22" model="cellStrat" out_fold="'$maleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr22.out
R CMD BATCH '--args data_fold="'$maleFolder'" chr="chrX" model="cellStrat" out_fold="'$maleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
mv $simFolder/*.* /gpfs42/projects/lab_helix_omics/shared_data/methExprsHELIX/$simFolder

# sim50
num=50
simFolder="$maleFolder/sim$num"
mkdir $simFolder

for i in {21..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$maleFolder'" chr="chr'$i'" model="cellStrat" out_fold="'$maleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$maleFolder'" chr="chrX" model="cellStrat" out_fold="'$maleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
mv $simFolder/*.* /gpfs42/projects/lab_helix_omics/shared_data/methExprsHELIX/$simFolder

## Female
femaleFolder="$resFolder/female"

# sim31
num=31
simFolder="$femaleFolder/sim$num"
mkdir $simFolder

for i in {11..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$femaleFolder'" chr="chr'$i'" model="cellStrat" out_fold="'$femaleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$femaleFolder'" chr="chrX" model="cellStrat" out_fold="'$femaleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
mv $simFolder/*.* /gpfs42/projects/lab_helix_omics/shared_data/methExprsHELIX/$simFolder

# sim34
num=34
simFolder="$femaleFolder/sim$num"
mkdir $simFolder

for i in {4..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$femaleFolder'" chr="chr'$i'" model="cellStrat" out_fold="'$femaleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$femaleFolder'" chr="chrX" model="cellStrat" out_fold="'$femaleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
mv $simFolder/*.* /gpfs42/projects/lab_helix_omics/shared_data/methExprsHELIX/$simFolder



## No cell adjusted
echo No Cell Adjusted
resFolder="results/MethComBatExpResidualsNoCellAdjStrat"

## Male
maleFolder="$resFolder/male"

# sim1
num=1
simFolder="$maleFolder/sim$num"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$maleFolder'" chr="chr'$i'" model="nocellStrat" out_fold="'$maleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$maleFolder'" chr="chrX" model="nocellStrat" out_fold="'$maleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
mv $simFolder /gpfs42/projects/lab_helix_omics/shared_data/methExprsHELIX/$resFolder

# sim45
num=45
simFolder="$maleFolder/sim$num"
mkdir $simFolder

R CMD BATCH '--args data_fold="'$maleFolder'" chr="chr3" model="nocellStrat" out_fold="'$maleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr3.out
mv $simFolder/*.* /gpfs42/projects/lab_helix_omics/shared_data/methExprsHELIX/$simFolder

## Female
femaleFolder="$resFolder/female"

# sim1
num=1
simFolder="$femaleFolder/sim$num"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$femaleFolder'" chr="chr'$i'" model="nocellStrat" out_fold="'$femaleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$femaleFolder'" chr="chrX" model="nocellStrat" out_fold="'$femaleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
mv $simFolder /gpfs42/projects/lab_helix_omics/shared_data/methExprsHELIX/$resFolder

# sim72
num=72
simFolder="$femaleFolder/sim$num"
mkdir $simFolder

for i in {14..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$femaleFolder'" chr="chr'$i'" model="nocellStrat" out_fold="'$femaleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$femaleFolder'" chr="chrX" model="nocellStrat" out_fold="'$femaleFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
mv $simFolder/*.* /gpfs42/projects/lab_helix_omics/shared_data/methExprsHELIX/$simFolder
