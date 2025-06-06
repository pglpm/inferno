## devtools::install()
library('inferno')

data(mtcars)

seed <- 1000

outputdir <- learn(
    data = mtcars,
    metadata = 'mtcars_metadata.csv',
    outputdir = 'mtcars_output',
    output = 'directory',
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    nsamples = 512,
    nchains = 4,
    parallel = 4,
    startupMCiterations = 1024,
    maxhours = 0/60,
    ## relerror = 0.02,
    ncheckpoints = NULL,
    cleanup = FALSE,
    ## miniter = 1200,
    prior = FALSE,
    showKtraces = TRUE,
    showAlphatraces = TRUE,
    seed = seed
    ## ,hyperparams = list(
    ##     ncomponents = 64,
    ##     minalpha = -4,
    ##     maxalpha = 4,
    ##     byalpha = 1,
    ##     Rshapelo = 0.5,
    ##     Rshapehi = 0.5,
    ##     Rvarm1 = 3^2,
    ##     Cshapelo = 0.5,
    ##     Cshapehi = 0.5,
    ##     Cvarm1 = 3^2,
    ##     Dshapelo = 0.5,
    ##     Dshapehi = 0.5,
    ##     Dvarm1 = 3^2,
    ##     Lshapelo = 0.5,
    ##     Lshapehi = 0.5,
    ##     Lvarm1 = 3^2,
    ##     Bshapelo = 1,
    ##     Bshapehi = 1,
    ##     Dthreshold = 1
    ## )
)
