startdir <- getwd()

#### Check and change working directory if necessary
if(basename(startdir) != 'custom'){
    cat('\nAre you in the correct folder?\n')
}

cat('\nInstalling local package "modelfreeinference" in local library\n\n')

if(!(Sys.getenv("R_LIBS_USER") %in% .libPaths())) {
    stop('Make sure your local installation directory,\n',
        Sys.getenv("R_LIBS_USER"),
        '\nexists.\n')
}

# devtools::install()
library('modelfreeinference')

seed <- 16

outputdirPrefix <- file.path('test_custom')

## ncores <- 4
## library('doParallel')
## mycluster <- makeCluster(ncores)
## registerDoParallel(mycluster)

currenttestdir <- inferpopulation(
    data = 'dataset_custom30.csv',
    metadata = 'metadata_custom.csv',
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
