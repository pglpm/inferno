#!/bin/bash

echo 'directory: '$1
echo 'parallel chains: '$2

workdir=$1
mkdir -p $workdir
cp mcmc_eeg.R $workdir/
cp smcmc_eeg.sbatch $workdir/
cd $workdir
sleep 3

sbatch --array=1-$2 smcmc_eeg.sbatch
