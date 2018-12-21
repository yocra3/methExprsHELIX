###############################################################################
# Commands to run analyses in cluster
#############################################################################

# Prepare objects

## No cell Adjusted
sbatch src/runNoCellAdjCluster.sh

### Stratified
sbatch src/runNoCellAdjStratCluster.sh

## Cell Adjusted
sbatch src/runCellAdjCluster.sh

### Stratified
sbatch src/runCellAdjStratCluster.sh
