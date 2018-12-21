###############################################################################
# Commands to run analyses in cluster
#############################################################################

# Prepare objects

## No cell Adjusted
sbatch src/runNoCellAdjCluster.sh

## Cell Adjusted
sbatch src/runCellAdjCluster.sh
