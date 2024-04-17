#!/bin/bash
## 1st arg: R script
## 2nd arg: output directory
## 3rd arg: first chain number or number of chains
## 4th arg: last chain number

if [ $# -lt 4 ]; then
    firstcore=1
    endcore=$3
else
    firstcore=$3
    endcore=$4
fi

echo 'script: '$1
echo 'directory: '$2
echo 'parallel chains: '$firstcore'-'$endcore
echo 'total '$((endcore - firstcore + 1))' chains'

workdir=$2
script=$1
mkdir -p $workdir
cp $script $workdir/basemcmcscript_$firstcore.R
cd $workdir
sleep 3

for ((core=$firstcore; core<=$endcore; core++))
do
   Rscript basemcmcscript_$firstcore.R $core > mcmc_out-$core.Rout 2>&1 &
done
