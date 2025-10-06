library('inferno')

dataset <- data.frame(V = 1:10)
metadata <- data.frame(name = 'V', type = 'continuous')

set.seed(16)
parallel <- 1

outputdir <- 'minitest'
learntdir <- learn(
    data = dataset,
    metadata = metadata,
    nsamples = 200,
    nchains = 2,
    minMCiterations = 200,
    outputdir = outputdir,
    output = 'directory',
    parallel = parallel,
    maxrelMCSE = +Inf,
    minESS = 100,
    maxhours = 0,
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
