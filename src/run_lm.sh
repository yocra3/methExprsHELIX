# Script for running the linear regressions between
# methylation and expression data using an R script

# Pasted directly to the terminal

# Linear reg 1

for i in {1..22}
do
echo $i
Rscript linear_m_subsets.R $i
done 2> x.log

# Linear reg 2

for i in {1..22}
do
echo $i
Rscript linear_m_subsets2.R $i
done 2> x2.log

# Linear reg 3

for i in {1..22}
do
echo $i
Rscript linear_m_subsets3.R $i
done 2> x3.log
