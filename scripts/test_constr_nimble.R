library('nimble')

rm(finitemixnimble,Cfinitemixnimble,confnimble,mcsampler,Cmcsampler)
##
finitemix <- nimbleCode({
    vmean ~ dnorm(mean=0, var=1)
    vvar ~ dinvgamma(shape=1, rate=1)
    for(d in 1:ndata){
        Laux[d] ~ dconstraint(Lcont[d] >= Lleft[d] & Lcont[d] <= Lright[d])
        Lcont[d] ~ dnorm(mean=vmean, var=vvar)
    }
})
##
datapoints <- list(
    Laux = rep(c(1L, 1L, 1L, NA),1)
   ,Lcont = rep(c(0, NA, NA, NA),1)
)
##
constants <- list(
    ndata=length(datapoints$Laux),
    Lleft = rep(c(-Inf, -Inf, -Inf, -Inf),1),
    Lright = rep(c(+Inf, -1, -1, +Inf),1)
)
##
initsFunction <- function(){
    list(
        Lcont = rep(c(NA, -2, -100, 0),1)
    )
}
##  
finitemixnimble <- nimbleModel(code=finitemix, name='test',
                               constants=constants,
                               inits=initsFunction(),
                               data=datapoints)

Cfinitemixnimble <- compileNimble(finitemixnimble, showCompilerOutput=FALSE)

confnimble <- configureMCMC(Cfinitemixnimble, monitors=c('Lcont','Laux','vmean','vvar'))
## confnimble$removeSamplers('Lcont[4]')
## confnimble$addSampler(target='Lcont[4]', type='posterior_predictive')
##
print(confnimble)
#3
targetslist <- sapply(confnimble$getSamplers(), function(xx)xx$target)
samplerorder <- c('vmean','vvar')
## neworder <- c(
##     foreach(var=setdiff(targetslist,samplerorder), .combine=c)%do%{grep(paste0('^',var,'$'),targetslist)},
##     foreach(var=samplerorder, .combine=c)%do%{grep(paste0('^',var,'(\\[.+\\])*$'),targetslist)}
neworder <- foreach(var=samplerorder, .combine=c)%do%{grep(paste0('^',var,'(\\[.+\\])*$'),targetslist)}
neworder <- c(setdiff(confnimble$getSamplerExecutionOrder(), neworder), neworder)
confnimble$setSamplerExecutionOrder(neworder)

##
mcsampler <- buildMCMC(confnimble)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

set.seed(890)
Cfinitemixnimble$setInits(initsFunction())
timecount <- Sys.time()
todelete <- Cmcsampler$run(niter=100000, thin=100, thin2=1, nburnin=0, reset=TRUE, resetMV=TRUE)
print(Sys.time()-timecount)
rm(todelete)
mcsamples <- as.matrix(Cmcsampler$mvSamples)


tplot(y=mcsamples[,'Lcont[4]'])

## bounded
## leftright as data
## > Time difference of 1.43506 mins
## leftright as init
## > Time difference of 1.43888 mins
## leftright as const
## > Time difference of 1.41759 mins

## ordinal
## cont as just init
## > Time difference of 43.3734 secs
## cont as data
## > Time difference of 45.2004 secs

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
