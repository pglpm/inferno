library('nimble')

## uncomment below for problematic W
testW <- unlist(read.csv('testW24.csv'))

## uncomment below for W that doesn't cause problems
## testW <- unlist(read.csv('testW23.csv'))

nclusters <- 64L
datapoints <- list(W = testW)
##
constants <- list(nclusters=nclusters)
##
initsFunction <- function(){
    catalpha0 <- rep(0, 11)
    catalpha0[4:7] <- rep(1/4, 4) # only values 4:7 have nonzero prob.
    dirchalpha0 <- matrix(2^((-5):5)/nclusters, nrow=11, ncol=nclusters)
    list(
        catalpha0 = catalpha0,
        dirchalpha0 = dirchalpha0,
        Alphaindex = 6
    )
}

print(initsFunction()$catalpha0)
## [1] 0.00 0.00 0.00 0.25 0.25 0.25 0.25 0.00 0.00 0.00 0.00


finitemix <- nimbleCode({
    Alphaindex ~ dcat(prob=catalpha0[1:11])
    W[1:nclusters] ~ ddirch(alpha=dirchalpha0[Alphaindex, 1:nclusters])
})

finitemixnimble <- nimbleModel(code=finitemix, name='finitemixnimble1',
                               constants=constants,
                               inits=initsFunction(),
                               data=datapoints
                               )

Cfinitemixnimble <- compileNimble(finitemixnimble, showCompilerOutput=FALSE)

confnimble <- configureMCMC(Cfinitemixnimble, monitors=c('Alphaindex', 'W'))
print(confnimble)

mcsampler <- buildMCMC(confnimble)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

## try different seeds or niter if problem doesn't appear
set.seed(890)
Cfinitemixnimble$setInits(initsFunction())
todelete <- Cmcsampler$run(niter=40, thin=1, thin2=1, nburnin=0, time=TRUE, reset=TRUE, resetMV=TRUE)
mcsamples <- as.matrix(Cmcsampler$mvSamples)


## using 'testW24.csv'
print(dcat(x=mcsamples[,'Alphaindex'], prob=initsFunction()$catalpha0))
## [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
