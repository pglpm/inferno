library('nimble')

rm(finitemixnimble,Cfinitemixnimble,confnimble,mcsampler,Cmcsampler)
##
finitemix <- nimbleCode({
    for(d in 1:ndata){
        Laux[d] ~ dconstraint(Lcont[d] >= Lleft[d] & Lcont[d] <= Lright[d])
        Lcont[d] ~ dnorm(mean=0, var=1)
    }
})
##
datapoints <- list(
    Laux = rep(c(1L, 1L, 1L, NA),10),
    Lleft = rep(c(-Inf, -Inf, -Inf, -Inf),10),
    Lright = rep(c(+Inf, -1, -1, +Inf),10)
    ,Lcont = rep(c(0, NA, NA, NA),10)
)
##
constants <- list(
    ndata=length(datapoints$Laux)
)
##
initsFunction <- function(){
    list(
        Lcont = rep(c(NA, -2, -100, 0),10)
    )
}
##  
finitemixnimble <- nimbleModel(code=finitemix, name='test',
                               constants=constants,
                               inits=initsFunction(),
                               data=datapoints)

Cfinitemixnimble <- compileNimble(finitemixnimble, showCompilerOutput=FALSE)

confnimble <- configureMCMC(Cfinitemixnimble, monitors=c('Lcont','Laux','Lleft','Lright'))
## confnimble$removeSamplers('Lcont[4]')
## confnimble$addSampler(target='Lcont[4]', type='posterior_predictive')
##
print(confnimble)
##
mcsampler <- buildMCMC(confnimble)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

set.seed(890)
Cfinitemixnimble$setInits(initsFunction())
timecount <- Sys.time()
todelete <- Cmcsampler$run(niter=1000000, thin=1000, thin2=1, nburnin=0, reset=TRUE, resetMV=TRUE)
print(Sys.time()-timecount)
rm(todelete)
mcsamples <- as.matrix(Cmcsampler$mvSamples)


tplot(y=mcsamples[,'Lcont[1]'])


> Time difference of 3.27896 secs


## uncomment below for problematic W
testW <- unlist(read.csv('testWbad.csv'))

## uncomment below for W that doesn't cause problems
## testW <- unlist(read.csv('testWgood.csv'))

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

confnimble$removeSamplers('Alphaindex')
confnimble$addSampler(target='Alphaindex', type='slice')

print(confnimble)

mcsampler <- buildMCMC(confnimble)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

## try different seeds or niter if problem doesn't appear
set.seed(890)
Cfinitemixnimble$setInits(initsFunction())
todelete <- Cmcsampler$run(niter=40, thin=1, thin2=1, nburnin=0, time=TRUE, reset=TRUE, resetMV=TRUE)
mcsamples <- as.matrix(Cmcsampler$mvSamples)


## using 'testWbad.csv'
print(dcat(x=mcsamples[,'Alphaindex'], prob=initsFunction()$catalpha0))
## [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

plot(x=1:nrow(mcsamples),y=mcsamples[,'Alphaindex'],type='l')
