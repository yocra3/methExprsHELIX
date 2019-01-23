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
 
## cell adjusted
echo Cell Adjusted
resFolder="results/MethComBatExpResidualsCellAdj"


### sim54
num=54
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {19..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


### sim55
num=55
simFolder="$resFolder/sim$num"
mkdir $simFolder

i=22
R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out

### sim56
num=56
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {14..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


### sim57
num=57
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {14..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


### sim58
num=58
simFolder="$resFolder/sim$num"
mkdir $simFolder

i=11
R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out

for i in {13..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


### sim59
num=59
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {9..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


### sim60
num=60
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {7..18}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


### sim61
num=61
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {6..18}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


### sim62
num=62
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {5..15}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


### sim63
num=63
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {2..9}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


### sim64
num=64
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..5}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


### sim65
num=65
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..4}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


### sim66
num=66
simFolder="$resFolder/sim$num"
mkdir $simFolder

i=1
R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out


## No cell adjusted
echo No Cell Adjusted
resFolder="results/MethComBatExpResidualsNoCellAdj"

# sim44
num=44
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..4}
do
 echo $i
 R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

## No cell adjusted Male
echo No Cell Adjusted Male
resFolder="results/MethComBatExpResidualsNoCellAdjStrat/male"
model="nocellStrat"


## sim51
num=51
simFolder="$resFolder/sim$num"
mkdir $simFolder

i=22
R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out