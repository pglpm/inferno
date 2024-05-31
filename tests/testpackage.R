source('bnpi.R')
testdir <- '../tests/'
outputdirPrefix <- paste0(paste0(testdir,'_testpackage_'), strftime(as.POSIXlt(Sys.time()), '%y%m%dT%H%M'))
test <- inferpopulation(data=paste0(testdir,'testdata.csv'), 
                        metadata=paste0(testdir,'metatestdata.csv'), 
                        outputdir=outputdirPrefix, 
                        nsamples=120, nchains=12, cleanup=F, parallel=4)
