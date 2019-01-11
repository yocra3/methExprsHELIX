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
num=79
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {19..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

### sim87
num=87
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..2}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="cell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


## cell adjusted Male
echo Cell Adjusted Male
resFolder="results/MethComBatExpResidualsCellAdjStrat/male"
model="cellStrat"

# sim90
num=90
simFolder="$resFolder/sim$num"
mkdir $simFolder


for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrY.out

# sim89
num=89
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrY.out


# sim88
num=88
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrY.out


# sim84
num=84
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {7..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrY.out


# sim93
num=93
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..10}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim91
num=91
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..11}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim92
num=92
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..11}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


# sim87
num=87
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {2..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrY.out


## No cell adjusted
echo No Cell Adjusted
resFolder="results/MethComBatExpResidualsNoCellAdj"

# sim42
num=42
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim40
num=40
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim41
num=41
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim64
num=64
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim63
num=63
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim23
num=23
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {11..16}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim59
num=59
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {5..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim37
num=37
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {5..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim30
num=30
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {2..3}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim43
num=43
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..5}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim55
num=55
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {9..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done


# sim66
num=66
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..21}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim65
num=65
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..18}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

# sim33
num=33
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {10..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="nocell" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

## No cell adjusted Male
echo No Cell Adjusted Male
resFolder="results/MethComBatExpResidualsNoCellAdjStrat/male"
model="nocellStrat"

## sim62
num=62
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrY.out

## sim60
num=60
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrY.out

## sim61
num=61
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrY.out

## sim59
num=59
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrY.out

## sim65
num=65
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {1..11}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done

## sim51
num=51
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrY.out

## sim52
num=52
simFolder="$resFolder/sim$num"
mkdir $simFolder

for i in {16..22}
do
  echo $i
  R CMD BATCH '--args data_fold="'$resFolder'" chr="chr'$i'" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out
done
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrY.out

## sim53
num=53
simFolder="$resFolder/sim$num"
mkdir $simFolder

R CMD BATCH '--args data_fold="'$resFolder'" chr="chr17" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchr$i.out

## sim50
num=50
simFolder="$resFolder/sim$num"
mkdir $simFolder

R CMD BATCH '--args data_fold="'$resFolder'" chr="chrX" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrX.out
R CMD BATCH '--args data_fold="'$resFolder'" chr="chrY" model="'$model'" out_fold="'$resFolder'" "'$num'"' src/runLinearModelSubset.R $simFolder/modchrY.out
