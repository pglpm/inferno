#!/bin/bash

#SBATCH --account=nn0000k
#SBATCH --job-name=Reval
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:55:00

set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## module restore
module --quiet purge
module load R/4.1.0-foss-2021a
## module load intel/2020b

cd ~/bayes_regression/finished2

Rscript test_distr.R > evaluation.Rout 2>&1
