library('inferno')

seed <- 16

#outputdirPrefix <- file.path('../tests/BN_functdependence')
pathToLearnt <- file.path('checkResults')

dataset <- read.csv('dataset_custom30.csv', na.strings='')
nv <- round(ncol(dataset)/2)
probs <- Pr(
    Y = dataset[1:20,1:nv,drop=F],
    X = dataset[1:20,(nv+1):ncol(dataset), drop=F],
    learnt = pathToLearnt,
    parallel = 4,
    silent = FALSE
)

print(probs$values[0])

print(probs$quantiles[0])