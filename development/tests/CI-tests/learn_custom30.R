library('inferno')

seed <- 16

outputdirPrefix <- file.path('_test_custom')

currenttestdir <- learn(
    data = 'dataset_custom30.csv',
    metadata = 'metadata_custom.csv',
    outputdir = outputdirPrefix,
    output = 'directory',
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    nsamples = 12,
    nchains = 2,
    parallel = 2,
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

mi <- mutualinfo(
    Y1names = c('N2vrt'),
    Y2names = c('Rvrt'),
    X = cbind(Bvrt = 'no'),
    learnt = currenttestdir,
    nsamples = 3600,
    parallel = 4
)

print(mi)
warnings()

cat('\nEnd\n')
