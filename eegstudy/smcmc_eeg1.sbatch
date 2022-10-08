#!/bin/bash

#SBATCH --account=nn8050k
#SBATCH --job-name=EEG1
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --time=8:00:00

set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## module restore
module --quiet purge
module load R/4.2.1-foss-2022a
## module load intel/2020b

cd /cluster/projects/nn8050k/eeg

Rscript mcmc_eeg1.R $SLURM_ARRAY_TASK_ID > mcmc_eeg1-$SLURM_ARRAY_TASK_ID.Rout 2>&1
