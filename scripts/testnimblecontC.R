library('data.table')
library('png')
library('foreach')
library('doFuture')
library('doRNG')
registerDoFuture()
cat('\navailableCores: ')
cat(availableCores())
cat('\navailableCores-multicore: ')
cat(availableCores('multicore'))
if(Sys.info()['nodename']=='luca-HP-Z2-G9'){
    ncores <- 20}else{
    ncores <- 6}
cat(paste0('\nusing ',ncores,' cores\n'))
if(ncores>1){
    if(.Platform$OS.type=='unix'){
        plan(multicore, workers=ncores)
    }else{
        plan(multisession, workers=ncores)
    }
}else{
    plan(sequential)
}
library('nimble')

extract <- function(x){
    grep(paste0('^',x,'(\\[.+\\])*$'), colnames(mcsamples))
}

set.seed(345)
npoints <- 200
## testp <- rbind(100+(1:npoints),rep(0,npoints))
testp <- matrix(rnorm(npoints,mean=c(20,0),sd=c(0.5,2)),nrow=2)
testp <- cbind(testp,-testp)
ang <- pi/6
testp <- matrix(c(cos(ang),sin(ang),sin(-ang),cos(ang)),nrow=2) %*% testp
tplot(x=testp[1,],y=testp[2,],type='p',pch='.',
      xlim=c(-1,1)*max(testp),
      ylim=c(-1,1)*max(testp)
      )

locdis <- t(apply(testp,1,function(xx){c(median(xx), mad(xx))}))
mtestp <- t(apply(testp,1,function(xx){(xx-median(xx))/mad(xx)}))

nclusters <- 64L
maxalpha <- 4
alpha0 <- 2^((-maxalpha+1):maxalpha)
nalpha <- length(alpha0)
npoints <- ncol(testp)
nvar <- nrow(testp)
shapehi0 <- 1.5
shapelo0 <- 2.5

Rmean1 <- 0#NULL
Rrate1 <- 1#NULL
Rvar1 <- 1#NULL

datapoints <- list(datapoints = mtestp)
##
constants <- list(npoints=npoints, nvar=nvar, nclusters=nclusters, nalpha=nalpha)
##
initsFunction <- function(){
    probalpha0 <- rep(1/nalpha, nalpha)
    walpha0 <- matrix(alpha0/nclusters, nrow=nalpha, ncol=nclusters)
    Alphaindex <- extraDistr::rcat(n=1, prob=probalpha0)
    W <- rdirch(n=1, alpha=walpha0[Alphaindex,])
    ##
    if((exists('Rmean1') && !is.null(Rmean1))){
        Rmean1 <- rep(Rmean1, nvar)
    }else{
        Rmean1 <- rnorm(n=nvar, mean=0, sd=1)
    }
    if((exists('Rrate1') && !is.null(Rrate1))){
        Rrate1 <- rep(Rrate1, nvar)
    }else{
        Rrate1 <- rinvgamma(n=nvar, shape=shapehi0, rate=1)
    }
    if((exists('Rvar1') && !is.null(Rvar1))){
        Rvar1 <- rep(Rvar1, nvar)
        }else{
            Rvar1 <- rinvgamma(n=nvar, shape=shapelo0, rate=Rrate1)
        }
    Rmean <- matrix(rnorm(n=nvar*nclusters, mean=Rmean1, sd=1), nrow=nvar, ncol=nclusters)
    Rrate <- matrix(rinvgamma(n=nvar*nclusters, shape=shapehi0, rate=Rvar1), nrow=nvar, ncol=nclusters)
    Rvar <- matrix(rinvgamma(n=nvar*nclusters, shape=shapelo0, rate=Rrate), nrow=nvar, ncol=nclusters)
    K <- rep(1, npoints)
    list(
        shapehi0 = shapehi0,
        shapelo0 = shapelo0,
        probalpha0 = probalpha0,
        walpha0 = walpha0,
        Alphaindex = round(nalpha/2),
        W = W,
        K = K,
        ##
        Rmean1 = Rmean1,
        Rrate1 = Rrate1,
        Rvar1 = Rvar1,
        Rmean = Rmean,
        Rrate = Rrate,
        Rvar = Rvar
    )
}

finitemix <- nimbleCode({
    Alphaindex ~ dcat(prob=probalpha0[1:nalpha])
    W[1:nclusters] ~ ddirch(alpha=walpha0[Alphaindex, 1:nclusters])
    ##
    for(v in 1:nvar){
        if(!(exists('Rmean1') && !is.null(Rmean1))){ Rmean1[v] ~ dnorm(mean=0, var=1) }
        if(!(exists('Rrate1') && !is.null(Rrate1))){ Rrate1[v] ~ dinvgamma(shape=shapehi0, rate=1) }
        if(!(exists('Rvar1') && !is.null(Rvar1))){ Rvar1[v] ~ dinvgamma(shape=shapelo0, rate=Rrate1[v]) }
    }
    for(k in 1:nclusters){
        for(v in 1:nvar){
            Rmean[v, k] ~ dnorm(mean=Rmean1[v], var=1)
            Rrate[v, k] ~ dinvgamma(shape=shapehi0, rate=Rvar1[v])
            Rvar[v, k] ~ dinvgamma(shape=shapelo0, rate=Rrate[v, k])
        }
    }
    for(d in 1:npoints){
        K[d] ~ dcat(prob=W[1:nclusters])
        ##
        for(v in 1:nvar){
            datapoints[v,d] ~ dnorm(mean=Rmean[v, K[d]], var=Rvar[v, K[d]])
        }
    }
})

finitemixnimble <- nimbleModel(code=finitemix, name='finitemixnimble1',
                               constants=constants,
                               inits=initsFunction(),
                               data=datapoints
                               )

Cfinitemixnimble <- compileNimble(finitemixnimble, showCompilerOutput=FALSE)

confnimble <- configureMCMC(Cfinitemixnimble,
                            monitors=c('Alphaindex', 'K',  'W', 'Rmean', 'Rvar', 'Rrate')
                            )
print(confnimble)

## confnimble$printSamplers(executionOrder=TRUE)
## confnimble$getSamplerExecutionOrder()
confnimble$setSamplerExecutionOrder(rev(confnimble$getSamplerExecutionOrder()))

mcsampler <- buildMCMC(confnimble)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

shapehi0 <- 2
shapelo0 <- 4
Cfinitemixnimble$setInits(initsFunction())
todelete <- Cmcsampler$run(niter=16000, thin=1, thin2=1, nburnin=15000, time=FALSE, reset=TRUE, resetMV=TRUE)
mcsamples <- as.matrix(Cmcsampler$mvSamples)

uk <- apply(mcsamples,1,function(rr){length(unique(rr[extract('K')]))})
table(uk)
table(alpha0[mcsamples[,extract('Alphaindex')]])


## pdff('testCnimble')
## tplot(y=mcsamples[,'Alphaindex'], ylab='Alphaindex')
## dev.off()


draws <- function(n,i){
    kk <- rcat(n=n, prob=mcsamples[i, extract('W')])
    means <- t(matrix(mcsamples[i, extract('Rmean')], nrow=nvar, ncol=nclusters)[,kk,drop=F])
    sds <- sqrt(t(matrix(mcsamples[i, extract('Rvar')], nrow=nvar, ncol=nclusters)[,kk,drop=F]))
    t(matrix(rnorm(n=2*n, mean=means, sd=sds), nrow=n, ncol=nvar))
}
##
pdff(paste0('testmcnohyper_lo',shapelo0,'_hi',shapehi0))
for(i in 1000:990){
testpoints <- draws(1e4,i)*locdis[,2]+locdis[,1]
xm <- tquant(testpoints[1,],c(2.5,97.5)/100)
ym <- tquant(testpoints[2,],c(2.5,97.5)/100)
tplot(x=list(testpoints[1,],testp[1,]), y=list(testpoints[2,],testp[2,]), type='p', pch=c(46,16), cex=c(1,0.5), alpha=c(0.9,0.8),
      xlim=xm, ylim=xm)
}
dev.off()
