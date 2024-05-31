source('bnpi.R')
testdir <- '../tests/'
test <- inferpopulation(data=paste0(testdir,'testdata.csv'), metadata=paste0(testdir,'metatestdata.csv'), outputdir=paste0(paste0(testdir,'testpackage_'),strftime(as.POSIXlt(Sys.time()), '%y%m%dT%H%M')), nsamples=120, nchains=12, cleanup=F, parallel=4)
