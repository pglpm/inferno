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

print('Result of package test (TRUE = passed):')
print(identical(currentFdistribution, readRDS(refFdistributionFile)))
