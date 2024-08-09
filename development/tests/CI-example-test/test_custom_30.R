library('modelfreeinference')

seed <- 16

outputdirPrefix <- file.path('_deletepackagetest')

## ncores <- 4
## library('doParallel')
## mycluster <- makeCluster(ncores)
## registerDoParallel(mycluster)

currenttestdir <- inferpopulation(
    data = 'data_test_custom_30.csv',
    metadata = 'metadata_test_custom.csv',
    outputdir = outputdirPrefix,
    output = 'directory',
    appendtimestamp = F,
    appendinfo = TRUE,
    nsamples = 120,
    nchains = 12,
    parallel = 4,
    maxhours = 1/60,
    ## relerror = 0.062,
    ncheckpoints = NULL,
    cleanup = FALSE,
    ## startupMCiterations = 1200,
    prior = FALSE,
    showKtraces = TRUE,
    showAlphatraces = TRUE,
    seed = seed
)

warnings()

cat('\nTest calculation of mutual information:\n')

mi <- mutualinfo(
    Y1names = c('N2vrt'),
    Y2names = c('Rvrt'),
    X = cbind(Bvrt = 'no'),
    mcoutput = currenttestdir,
    nsamples = 3600,
    parallel = 4
)

print(mi)

warnings()

cat('\nEnd\n')