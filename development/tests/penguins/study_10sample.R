library('inferno')

## metadatatemplate(penguins, file='meta_penguins')

seed <- 16
parallel <- 4

outputdir <- '__learn_penguin_10'
learntdir <- learn(
    data = 'penguin_data10.csv',
    prior = FALSE,
    metadata = 'penguin_metadata.csv',
    outputdir = outputdir,
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    outputvalue = 'directory',
    parallel = parallel,
    relerror = 0.038,
    ## parameters for short test run:
    ## subsampledata = 10,
    ## maxhours = 0,
    ## nsamplesperchain = 60,
    ## nchains = parallel + 1,
    ##
    ## hyperparams = list(
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
    ##     Dthreshold = 1,
    ##     tscalefactor = 2,
    ##     initmethod = 'prior'
    ##     ## precluster, prior
    ## ),
    seed = seed
)
