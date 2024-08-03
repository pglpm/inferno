startdir <- getwd()

#### Check and change working directory if necessary
if(basename(startdir) != 'test_custom'){
    cat('\nAre you in the correct folder?\n')
}

cat('\nInstalling local package "modelfreeinference" in local library\n\n')

if(!(Sys.getenv("R_LIBS_USER") %in% .libPaths())) {
    stop('Make sure your local installation directory,\n',
        Sys.getenv("R_LIBS_USER"),
        '\nexists.\n')
}

devtools::install()
library('modelfreeinference')

seed <- 16

outputdirPrefix <- file.path('_oneVpackagetest')

currenttestdir <- inferpopulation(data = 'data_test_custom_30.csv',
    metadata = 'metadata_test_R.csv',
    outputdir = outputdirPrefix,
    output = 'directory',
    appendtimestamp = FALSE,
    appendinfo = TRUE,
    nsamples = 30,
    nchains = 3,
    parallel = 1,
    relerror = 0.05,
    ncheckpoints = NULL,
    cleanup = FALSE,
    minMCiterations = 1200,
    ## prior = TRUE,
    showKtraces = T,
    showAlphatraces = T,
    seed = seed
)
