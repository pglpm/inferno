#!/bin/bash

#SBATCH --account=nn0000k
#SBATCH --job-name=MCMCold11
#SBATCH --partition=normal
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=48:00:00

set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## module restore
module --quiet purge
module load R/4.1.0-foss-2021a
## module load intel/2020b

cd ~/bayes_regression

Rscript rmsd_regression_11_C21D21-sbatch.R > mcmc_old11_C21D21.Rout 2>&1
