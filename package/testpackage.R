source('bnpi.R')
test <- inferpopulation(data='testdata.csv', metadata='metatestdata.csv', outputdir=paste0('testpackage_',strftime(as.POSIXlt(Sys.time()), '%y%m%dT%H%M')), nsamples=120, nchains=12, cleanup=F, parallel=4)
