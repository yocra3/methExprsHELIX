#!/bin/sh

###############################################################################
# Run linear models for all the different models
###############################################################################

#set the job name
#SBATCH --job-name=nocellLM

#set the number of CPUS per task
#SBATCH --ntasks-per-node=16

#set the memory
##SBATCH --mem=3000

# job output file information
#SBATCH -o nocellLM.out

# job errors file
#SBATCH -e nocellLM.err

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

## No cell adjusted
echo No Cell Adjusted
resFolder="results/MethComBatExpResidualsNoCellAdj"
for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'"' src/runLinearModelSubset.R $resFolder/modchr$i.out
done