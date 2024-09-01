library('inferno')
# Checking that everything works together

seed <- 16

outputdir <- file.path('_fullchain')
datafile <- 'dataset_custom30.csv'
metadatafile <- 'metadata_custom.csv'
dataset <- read.csv('dataset_custom30.csv', na.strings = '')
silent <- TRUE

currenttestdir <- learn(
    data = datafile,
    metadata = metadatafile,
    outputdir = outputdir,
    output = 'directory',
    nsamples = 12,
    nchains = 2,
    parallel = 2,
    appendtimestamp = FALSE,
    appendinfo = FALSE,
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

mi <- mutualinfo(
    Y1names = c('N2vrt'),
    Y2names = c('Rvrt'),
    X = cbind(Bvrt = 'no'),
    learnt = outputdir,
    nsamples = 3600,
    parallel = 4,
    silent = silent
)

nv <- round(ncol(dataset) / 2)
probs <- Pr(
    Y = dataset[1:20, 1:nv, drop = FALSE],
    X = dataset[1:20, (nv + 1):ncol(dataset), drop = FALSE],
    learnt = outputdir,
    parallel = 4,
    silent = silent
)

warnings()
