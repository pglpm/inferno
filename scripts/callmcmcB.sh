#!/bin/bash

if [ $# -lt 3 ]; then
    firstcore=1
    endcore=$2
else
    firstcore=$2
    endcore=$3
fi

echo 'directory: '$1
echo 'parallel chains: '$firstcore'-'$endcore
echo 'total '$((endcore - firstcore + 1))' chains'

workdir=$1
mkdir -p $workdir
cp mcmc_2212120902.R $workdir/basemcmcscript_$firstcore.R
cd $workdir
sleep 3

for ((core=$firstcore; core<=$endcore; core++))
do
   Rscript basemcmcscript_$firstcore.R $core > mcmc_out-$core.Rout 2>&1 &
done
