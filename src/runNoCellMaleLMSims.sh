#!/bin/sh

###############################################################################
# Run linear models for all the different models
###############################################################################

#set the job name
#SBATCH --job-name=noCellMaleSims

#set the number of CPUS per task
#SBATCH --ntasks-per-node=16

#set the memory
#SBATCH --mem=35000

# job output file information
#SBATCH -o noCellMaleSims.out

# job errors file
#SBATCH -e noCellMaleSims.err

# set the partition where the job will run
#SBATCH --partition=normal
 
# set the number of nodes
#SBATCH --nodes=1
 
# set max wallclock time
#SBATCH --time=7:00:00
 
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=END
 
# send mail to this address
#SBATCH --mail-user= carlos.ruiz@isglobal.org

module load R/3.5.1-foss-2018b


## No cell adjusted Male
echo No Cell Adjusted Male
resFolder="results/MethComBatExpResidualsNoCellAdjStrat/male"
model="nocellStrat"

simFolder="$resFolder/sim${SLURM_ARRAY_TASK_ID}"
mkdir $simFolder


for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" ${SLURM_ARRAY_TASK_ID}' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" ${SLURM_ARRAY_TASK_ID}' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" ${SLURM_ARRAY_TASK_ID}' src/runLinearModelSubset.R $simFolder/modchrY.out
mv $simFolder /gpfs42/projects/lab_helix_omics/shared_data/methExprsHELIX/$resFolder