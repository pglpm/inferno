#!/bin/bash

#SBATCH --account=nn0000k
#SBATCH --job-name=R_script_test
#SBATCH --partition=devel
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --time=20:00

set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

module restore
module load R/4.1.0-foss-2021a
module load intel/2020b

cd ~/bayes_regression

Rscript rmsd_regression_sbatch_1.R > output.Rout 2>&1