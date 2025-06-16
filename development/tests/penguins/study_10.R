library('inferno')

## metadatatemplate(penguins, file='meta_penguins')

seed <- 16
parallel <- 15

outputdir <- '__learn_penguins_10'
learntdir <- learn(
    data = 'penguin_data10.csv',
    prior = FALSE,
    metadata = 'penguin_metadata.csv',
    outputdir = outputdir,
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    output = 'directory',
    parallel = parallel,
    ## parameters for short test run:
    ## subsampledata = 10,
    ## maxhours = 0,
    ## nsamplesperchain = 60,
    ## nchains = parallel + 1,
    ##
    seed = seed
)
