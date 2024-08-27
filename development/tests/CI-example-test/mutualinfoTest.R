library('inferno')

seed <- 16

#outputdirPrefix <- file.path('../tests/BN_functdependence')
pathToLearnt <- file.path('checkResults')

mi <- mutualinfo(
    Y1names = c('N2vrt'),
    Y2names = c('Rvrt'),
    X = cbind(Bvrt = 'no'),
    learnt = pathToLearnt,
    nsamples = 3600,
    parallel = 2,
    silent = FALSE
)

print(mi)