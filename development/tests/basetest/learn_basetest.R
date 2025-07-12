library('inferno')

set.seed(16)
parallel <- 2

outputdir <- '__testbase_ess100'
learntdir <- learn(
    data = 'data_basetest.csv',
    metadata = 'metadata_basetest.csv',
    nsamples = 3600,
    nchains = parallel,
    ## minMCiterations = 3600 * 3,
    prior = FALSE,
    outputdir = outputdir,
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    cleanup = FALSE,
    output = 'directory',
    parallel = parallel,
    maxrelMCSE = +Inf,
    minESS = 100,
    ## ncheckpoints = 12,
    ##
    ## ## parameters for short test run:
    ## subsampledata = 10,
    ## maxhours = 0,
    ## nsamplesperchain = 60,
    ## nchains = parallel + 1,
    ##
    hyperparams = list(
        ## ncomponents = 64,
        ## minalpha = -4,
        ## maxalpha = 4,
        ## byalpha = 1,
        ## Rshapelo = 0.5,
        ## Rshapehi = 0.5,
        ## Rvarm1 = 8^2,
        ## Cshapelo = 0.5,
        ## Cshapehi = 0.5,
        ## Cvarm1 = 8^2,
        ## Dshapelo = 0.5,
        ## Dshapehi = 0.5,
        ## Dvarm1 = 8^2,
        ## Bshapelo = 1,
        ## Bshapehi = 1,
        ## Dthreshold = 1,
        ## tscalefactor = 1,
        ## initmethod = 'allinone',
        ## avoidzeroW = TRUE
        ## precluster, prior
    )
)
