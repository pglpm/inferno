
metadatatemplate(penguins, file='meta_penguins')

seed <- 16
parallel <- 8

outputdir <- 'learn_penguins'
learntdir <- learn(
    data = penguins,
    prior = FALSE,
    metadata = 'meta_penguins.csv',
    outputdir = outputdir,
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    output = 'directory',
    parallel = parallel,
    ## parameters for short test run:
    ## subsampledata = 10,
    ## maxhours = 0,
    ## nsamples = 60 * parallel,
    ## nchains = parallel,
    ##
    seed = seed
)

