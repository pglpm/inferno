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

outputdirPrefix <- file.path('_deletetest_packagetest')

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
    relerror = 0.062,
    ncheckpoints = NULL,
    cleanup = FALSE,
    miniter = 1200,
    ## prior = TRUE,
    showKtraces = T,
    showAlphatraces = T,
    seed = seed
)

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

cat('\nEnd\n')



if(FALSE){
    ## The maths of the new package version has changed a little,
    ## so a comparison with old-version results are not meaningful
    refdir <- 'reference_seed16-vrt9_dat15_smp120'

#### Test whether Fdistribution output is identical
    cat('\nVerifying equality of "Fdistribution.rds" (TRUE = passed):\n')
    print(identical(
        readRDS(file.path(currenttestdir,'Fdistribution.rds')),
        readRDS(file.path(refdir,'Fdistribution.rds'))
    ))

#### Test whether MCtraces output is identical
    cat('\nVerifying equality of "MCtraces.rds" (TRUE = passed):\n')
    print(identical(
        readRDS(file.path(currenttestdir,'MCtraces.rds')),
        readRDS(file.path(refdir,'MCtraces.rds'))
    ))

#### Test whether computation log-1 is identical
    ## remove lines containing diagnostic times, as they can vary
    cat('\nVerifying equality of "log-1" (TRUE = passed):\n')
    reffile <- readLines(file.path(refdir,'log-1.txt'))
    reffile <- reffile[!grepl('time', reffile, fixed=TRUE)]
    currentfile <- readLines(file.path(currenttestdir,'log-1.log'))
    currentfile <- currentfile[!grepl('time', currentfile, fixed=TRUE)]

    print(identical(currentfile, reffile))
}
