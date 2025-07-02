library('inferno')

set.seed(16)
parallel <- 6

outputdir <- '__gauss2B_nDN_cp12'
learntdir <- learn(
    data = 'data_gaussians2B.csv',
    metadata = 'metadata_gaussians2B.csv',
    prior = FALSE,
    maxrelMCSE = +Inf,
    minESS = NULL,
    outputdir = outputdir,
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    output = 'directory',
    cleanup = TRUE,
    parallel = parallel,
    ## subsampledata = 30,
    ncheckpoints = 12,
    ## parameters for short test run:
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
        ## Rvarm1 = 3^2,
        ## Cshapelo = 0.5,
        ## Cshapehi = 0.5,
        ## Cvarm1 = 3^2,
        ## Dshapelo = 0.5,
        ## Dshapehi = 0.5,
        ## Dvarm1 = 3^2,
        ## Bshapelo = 1,
        ## Bshapehi = 1,
        ## Dthreshold = 1,
        ## tscalefactor = 1.35,
        ## avoidzeroW = NULL
        ## initmethod = 'precluster'
    )
)


file.copy('../../../R/learn.R', learntdir)

