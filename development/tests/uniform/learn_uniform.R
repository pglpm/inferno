startdir <- getwd()

#### Check and change working directory if necessary
if(basename(startdir) != 'uniform'){
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

outputdirPrefix <- file.path('test_uniform')

## ncores <- 4
## library('doParallel')
## mycluster <- makeCluster(ncores)
## registerDoParallel(mycluster)

currenttestdir <- learn(
    data = 'uniformdata.csv',
    metadata = 'uniformmetadata.csv',
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

## mi <- mutualinfo(
##     Y1names = c('N2vrt'),
##     Y2names = c('Rvrt'),
##     X = cbind(Bvrt = 'no'),
##     learnt = currenttestdir,
##     nsamples = 3600,
##     parallel = 4
## )
## 
## print(mi)
## 
## warnings()

dataset <- read.csv('uniformdata.csv', na.strings='')
probs <- Pr(
    Y = dataset[1:20,1,drop=F],
    X = NULL,
    learnt = currenttestdir,
    parallel = 4
)

print(probs$values)

print(probs$quantiles)

warnings()


cat('\nEnd\n')
