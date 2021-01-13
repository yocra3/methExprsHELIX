###############################################################################
# Commands to run analyses in cluster
#############################################################################

# Prepare objects
## Send files (repeat per folder)
rsync -azvh --progress -e 'ssh -o "ProxyCommand ssh -A cruizh@sit-web.upf.edu -W %h:%p"' *chr*.RData cruizh@marvin.s.upf.edu:/scratch/lab_helix_omics/cruizh/methExprsHELIX/results/MethComBatExpResidualsCellAdj
rsync -azvh --progress -e 'ssh -o "ProxyCommand ssh -A cruizh@sit-web.upf.edu -W %h:%p"' pheno.Rdata cruizh@marvin.s.upf.edu:/scratch/lab_helix_omics/cruizh/methExprsHELIX/results/MethComBatExpResidualsCellAdj

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

## Compress results files in a folder
for i in `(find -name outputchr*.txt)`
do
 gzip $i
done

## Failed simulations
sbatch src/runFailedSimulations.sh
sbatch src/runFailedSimulations2.sh

## Send files
rsync -azvh --progress -e 'ssh -o "ProxyCommand ssh -A cruizh@sit-web.upf.edu -W %h:%p"' . cruizh@marvin.s.upf.edu:/scratch/lab_helix_omics/cruizh/methExprsHELIX/results/
## Receive files
rsync -azvh --progress -e 'ssh -o "ProxyCommand ssh -A cruizh@sit-web.upf.edu -W %h:%p"' cruizh@marvin.s.upf.edu:/gpfs42/projects/lab_helix_omics/shared_data/methExprsHELIX/results .