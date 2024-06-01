startdir <- getwd()

if(basename(startdir) == 'tests'){
    setwd('../package/')
    cat('\nSwitching to "package" directory\n')
} else if(basename(startdir) == 'bayes_nonparametric_inference'){
    setwd('package/')
    cat('\nSwitching to "package" directory\n')
} else if(basename(startdir) != 'package'){
    stop('Please run this script from directory "package" or "tests" or base')
}

source('bnpi.R')

testdir <- '../tests/'

refFdistributionFile <- paste0(testdir,
                               'reference_packagetest_seed16_240601T0912-V8-D15-K64-S120/Fdistribution.rds')

seed <- 16

outputdirPrefix <- paste0(paste0(testdir,'__packagetest-'),
                          strftime(as.POSIXlt(Sys.time()), '%y%m%dT%H%M'))

test <- inferpopulation(data=paste0(testdir,'testdata.csv'),
                        metadata=paste0(testdir,'metatestdata.csv'),
                        outputdir=outputdirPrefix,
                        nsamples=120, nchains=12, cleanup=F, parallel=4,
                        seed=seed)

currentFdistribution <- readRDS(paste0(outputdirPrefix ,
                                       '-V8-D15-K64-S120/Fdistribution.rds'))

cat('\nResult of package test (TRUE = passed):\n')
print(identical(currentFdistribution, readRDS(refFdistributionFile)))

## return to original directory, in case called with 'source'
cat('\nSwitching to original directory\n')
setwd(startdir)

