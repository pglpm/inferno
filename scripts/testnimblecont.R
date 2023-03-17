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
npoints <- 100
## testp <- rbind(100+(1:npoints),rep(0,npoints))
testp <- matrix(rnorm(npoints,mean=c(20,0),sd=c(0.5,2)),nrow=2)
testp <- cbind(testp,-testp)
ang <- pi/6
testp <- matrix(c(cos(ang),sin(ang),sin(-ang),cos(ang)),nrow=2) %*% testp
tplot(x=testp[1,],y=testp[2,],type='p',pch='.',
      xlim=c(-1,1)*max(testp),
      ylim=c(-1,1)*max(testp)
      )

nclusters <- 64L
maxalpha <- 4
alpha0 <- 2^((-maxalpha+1):maxalpha)
nalpha <- length(alpha0)
npoints <- ncol(testp)
nvar <- nrow(testp)
shapehi0 <- 1
shapelo0 <- 1

datapoints <- list(datapoints = testp)
##
constants <- list(npoints=npoints, nvar=nvar, nclusters=nclusters, nalpha=nalpha)
##
initsFunction <- function(){
    probalpha0 <- rep(1/nalpha, nalpha)
    walpha0 <- matrix(alpha0/nclusters, nrow=nalpha, ncol=nclusters)
    Alphaindex <- extraDistr::rcat(n=1, prob=probalpha0)
    W <- rdirch(n=1, alpha=walpha0[Alphaindex,])
    ##
    Rmean <- matrix(rnorm(n=nvar*nclusters, mean=0, sd=1), nrow=nvar, ncol=nclusters)
    Rrate <- rinvgamma(n=nclusters, shape=shapehi0, rate=1)
    Rvar <- rinvgamma(n=nclusters, shape=shapelo0, rate=Rrate)
    K <- rep(1, npoints)
    list(
        shapehi0 = shapehi0,
        shapelo0 = shapelo0,
        probalpha0 = probalpha0,
        walpha0 = walpha0,
        Alphaindex = Alphaindex,
        W = W,
        K = K,
        ##
        Rmean = Rmean,
        Rrate = Rrate,
        Rvar = Rvar
    )
}

finitemix <- nimbleCode({
    Alphaindex ~ dcat(prob=probalpha0[1:nalpha])
    W[1:nclusters] ~ ddirch(alpha=walpha0[Alphaindex, 1:nclusters])
    ##
    for(k in 1:nclusters){
        for(v in 1:nvar){
            Rmean[v, k] ~ dnorm(mean=0, var=1)
        }
        Rrate[k] ~ dinvgamma(shape=shapehi0, rate=1)
        Rvar[k] ~ dinvgamma(shape=shapelo0, rate=Rrate[k])
    }
    for(d in 1:npoints){
        K[d] ~ dcat(prob=W[1:nclusters])
        ##
        for(v in 1:nvar){
            datapoints[v,d] ~ dnorm(mean=Rmean[v, K[d]], var=Rvar[K[d]])
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

Cfinitemixnimble$setInits(initsFunction())
todelete <- Cmcsampler$run(niter=10000, thin=1, thin2=1, nburnin=9000, time=FALSE, reset=TRUE, resetMV=TRUE)
mcsamples <- as.matrix(Cmcsampler$mvSamples)

uk <- apply(mcsamples,1,function(rr){length(unique(rr[extract('K')]))})
summary(uk)


pdff('testCnimble')
tplot(y=mcsamples[,'Alphaindex'], ylab='Alphaindex')
dev.off()


draws <- function(n,i){
    kk <- rcat(n=n, prob=mcsamples[i, extract('W')])
    means <- t(matrix(mcsamples[i, extract('Rmean')], nrow=nvar, ncol=nclusters)[,kk,drop=F])
    sds <- sqrt(mcsamples[i, extract('Rvar')][kk])
    matrix(rnorm(n=2*n, mean=means, sd=sds), nrow=n, ncol=nvar)
}

pdff()
for(i in 1000:990){
testpoints <- draws(1e4,i)
xm <- c(-1,1)*max(abs(testpoints))
tplot(x=list(testpoints[,1],testp[1,]), y=list(testpoints[,2],testp[2,]), type='p', pch=c(46,16), cex=1, alpha=0.9,
      xlim=xm, ylim=xm)
}
dev.off()
