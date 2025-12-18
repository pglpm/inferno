library('inferno')

dataset <- data.frame(V = letters[1:5])
metadata <- data.frame(name = 'V', type = 'nominal', V1 = letters[1], V2 = letters[2], V3 = letters[3], V4 = letters[4], V5 = letters[5])

set.seed(16)
parallel <- 2

outputdir <- '__minitest_nom'
learntdir <- learn(
    data = dataset,
    metadata = metadata,
    nsamples = 200,
    nchains = 2,
    minMCiterations = 200,
    outputdir = outputdir,
    outputvalue = 'outputdir',
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
        ## avoidzeroW = FALSE
        ## precluster, prior
    )
)
