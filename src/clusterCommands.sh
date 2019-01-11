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

# Run models

## Run base models
sbatch src/runNoCellLM.sh
sbatch src/runNoCellFemaleLM.sh
sbatch src/runNoCellMaleLM.sh
sbatch src/runCellAdjLM.sh
sbatch src/runCellAdjMaleLM.sh
sbatch src/runCellAdjFemaleLM.sh

## Simulations
sbatch --array=1-100 src/runNoCellLMSims.sh
sbatch --array=1-100 src/runNoCellMaleLMSims.sh
sbatch --array=1-100 src/runNoCellFemaleSims.sh
sbatch --array=1-100 src/runCellLMSims.sh
sbatch --array=1-100 src/runCellMaleLMSims.sh
sbatch --array=1-100 src/runCellFemaleSims.sh

## Failed simulations
sbatch src/runFailedSimulations.sh
