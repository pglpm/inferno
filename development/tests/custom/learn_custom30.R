startdir <- getwd()

#### Check and change working directory if necessary
if(basename(startdir) != 'custom'){
    cat('\nAre you in the correct folder?\n')
}

cat('\nLoading package "inferno"\n\n')

if(!(Sys.getenv("R_LIBS_USER") %in% .libPaths())) {
    stop('Make sure your local installation directory,\n',
        Sys.getenv("R_LIBS_USER"),
        '\nexists.\n')
}

# devtools::install()
library('inferno')

seed <- 16

outputdirPrefix <- file.path('__test_custom')

## ncores <- 4
## library('doParallel')
## mycluster <- makeCluster(ncores)
## registerDoParallel(mycluster)

currenttestdir <- learn(
    data = 'dataset_custom30.csv',
    metadata = 'metadata_custom.csv',
    outputdir = outputdirPrefix,
    output = 'directory',
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    nsamples = 120,
    nchains = 8,
    parallel = 4,
    maxhours = 1/60,
    ## relerror = 0.062,
    ncheckpoints = NULL,
    cleanup = FALSE,
    startupMCiterations = 1200,
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

dataset <- read.csv('dataset_custom30.csv', na.strings='')
nv <- round(ncol(dataset)/2)
probs <- Pr(
    Y = dataset[1:20,1:nv,drop=F],
    X = dataset[1:20,(nv+1):ncol(dataset), drop=F],
    learnt = currenttestdir,
    parallel = 4
)

print(probs$values)

print(probs$quantiles)

warnings()


cat('\nEnd\n')
