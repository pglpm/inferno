#library('prova')
devtools::load_all()

set.seed(16)
dataset <- data.frame(V = rnorm(n = 3))
metadata <- data.frame(name = 'V', type = 'continuous')
learnt <- learn(
    data = dataset, metadata = metadata,
    ## the following parameters are unrealistic
    ## only used to reduce computation time for this example
    nsamples = 10, nchains = 1,
    startupMCiterations = 10, maxMCiterations = 10,
    minESS = 0, initES = 0
)

