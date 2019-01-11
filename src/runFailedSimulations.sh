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

### sim79
simFolder="$resFolder/sim79"
mkdir $simFolder

for i in {19..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

### sim87
simFolder="$resFolder/sim87"
mkdir $simFolder

for i in {1..2}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


## cell adjusted Male
echo Cell Adjusted Male
resFolder="results/MethComBatExpResidualsCellAdjStrat/male"
model="cellStrat"

# sim90
simFolder="$resFolder/sim90"
mkdir $simFolder


for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrY.out

# sim89
simFolder="$resFolder/sim89"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrY.out


# sim88
simFolder="$resFolder/sim88"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrY.out


# sim84
simFolder="$resFolder/sim84"
mkdir $simFolder

for i in {7..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrY.out


# sim93
simFolder="$resFolder/sim93"
mkdir $simFolder

for i in {1..10}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim91
simFolder="$resFolder/sim91"
mkdir $simFolder

for i in {1..11}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim92
simFolder="$resFolder/sim92"
mkdir $simFolder

for i in {1..11}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


# sim87
simFolder="$resFolder/sim87"
mkdir $simFolder

for i in {2..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrY.out


## No cell adjusted
echo No Cell Adjusted
resFolder="results/MethComBatExpResidualsNoCellAdj"

# sim42
simFolder="$resFolder/sim42"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim40
simFolder="$resFolder/sim40"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim41
simFolder="$resFolder/sim41"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim64
simFolder="$resFolder/sim64"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim63
simFolder="$resFolder/sim63"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim23
simFolder="$resFolder/sim23"
mkdir $simFolder

for i in {11..16}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim59
simFolder="$resFolder/sim59"
mkdir $simFolder

for i in {5..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim37
simFolder="$resFolder/sim37"
mkdir $simFolder

for i in {5..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim30
simFolder="$resFolder/sim30"
mkdir $simFolder

for i in {2..3}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim43
simFolder="$resFolder/sim43"
mkdir $simFolder

for i in {1..5}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim55
simFolder="$resFolder/sim55"
mkdir $simFolder

for i in {9..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


# sim66
simFolder="$resFolder/sim66"
mkdir $simFolder

for i in {1..21}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim65
simFolder="$resFolder/sim65"
mkdir $simFolder

for i in {1..18}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim33
simFolder="$resFolder/sim33"
mkdir $simFolder

for i in {10..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

## No cell adjusted Male
echo No Cell Adjusted Male
resFolder="results/MethComBatExpResidualsNoCellAdjStrat/male"
model="nocellStrat"

## sim62
simFolder="$resFolder/sim62"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrY.out

## sim60
simFolder="$resFolder/sim60"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrY.out

## sim61
simFolder="$resFolder/sim61"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrY.out

## sim59
simFolder="$resFolder/sim59"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrY.out

## sim65
simFolder="$resFolder/sim65"
mkdir $simFolder

for i in {1..11}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

## sim51
simFolder="$resFolder/sim51"
mkdir $simFolder

for i in {22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrY.out

## sim52
simFolder="$resFolder/sim52"
mkdir $simFolder

for i in {16..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrY.out

## sim53
simFolder="$resFolder/sim53"
mkdir $simFolder

R CMD BATCH '--args data_fold="'$resFolder'" chr="chr17" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchr$i.out

## sim50
simFolder="$resFolder/sim50"
mkdir $simFolder

R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'${SLURM_ARRAY_TASK_ID}'"' src/runLinearModelSubset.R $simFolder/modchrY.out
