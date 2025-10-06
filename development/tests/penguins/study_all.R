library('inferno')

## metadatatemplate(penguins, file='meta_penguins')

seed <- 16
parallel <- 4

outputdir <- '__learn_penguins_all'
learntdir <- learn(
    data = 'penguin_data.csv',
    prior = FALSE,
    metadata = 'penguin_metadata.csv',
    outputdir = outputdir,
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    outputvalue = 'directory',
    parallel = parallel,
    ## parameters for short test run:
    ## subsampledata = 10,
    ## maxhours = 0,
    ## nsamplesperchain = 60,
    ## nchains = parallel + 1,
    ##
    seed = seed
)
