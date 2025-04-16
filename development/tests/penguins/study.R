unloadNamespace('inferno')
library('inferno')


## metadatatemplate(penguins, file='meta_penguins')

seed <- 16
parallel <- 6

outputdir <- '__learn_penguins'
learntdir <- learn(
    data = penguins,
    prior = FALSE,
    metadata = 'meta_penguins2.csv',
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


testf <- function(){
cl <- parallel::makeCluster(3)
doParallel::registerDoParallel(cl)
suppressMessages(
    foreach(i = 1:3) %dopar% {
    library('nimble')
    nimble::messageIfVerbose('test')
    i
    }
)
foreach::registerDoSEQ()
parallel::stopCluster(cl)
}

testf()
