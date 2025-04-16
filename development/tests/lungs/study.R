unloadNamespace('inferno')
library('inferno')

seed <- 16
parallel <- 6

outputdir <- '__learn_lungs'
learntdir <- learn(
    data = 'calibration_train_cleaned.csv',
    prior = FALSE,
    metadata = 'meta_calibration_train.csv',
    outputdir = outputdir,
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    output = 'directory',
    parallel = parallel,
    ## parameters for short test run:
    subsampledata = 10,
    maxhours = 0,
    nsamplesperchain = 60,
    nchains = parallel,
    ##
    seed = seed
)
