startdir <- getwd()

#### Check and change working directory if necessary
if(basename(startdir) != 'custom'){
    cat('\nAre you in the correct folder?\n')
}

cat('\nInstalling local package "predict" in local library\n\n')

if(!(Sys.getenv("R_LIBS_USER") %in% .libPaths())) {
    stop('Make sure your local installation directory,\n',
        Sys.getenv("R_LIBS_USER"),
        '\nexists.\n')
}

devtools::install()
library('predict')

seed <- 16

outputdirPrefix <- file.path('_newdeletepackagetest')

## ncores <- 4
## library('doParallel')
## mycluster <- makeCluster(ncores)
## registerDoParallel(mycluster)

currenttestdir <- learn(
    data = 'dataset_custom500.csv',
    metadata = 'metadata_custom.csv',
    outputdir = outputdirPrefix,
    output = 'directory',
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    nsamples = 3600,
    nchains = 60,
    parallel = 8,
    maxhours = +Inf,
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
