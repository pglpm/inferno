#!/bin/bash

#SBATCH --account=nn0000k
#SBATCH --job-name=R_script_test
#SBATCH --partition=normal
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --time=28:00:00

set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## module restore
module --quiet purge
module load R/4.1.0-foss-2021a
## module load intel/2020b

cd ~/bayes_regression

Rscript rmsd_regression_4-sbatch.R > output_4.Rout 2>&1
