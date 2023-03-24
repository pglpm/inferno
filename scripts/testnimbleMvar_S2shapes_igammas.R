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
## nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)

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

locdis <- t(apply(testp,1,function(xx){c(median(xx), mad(xx))}))
mtestp <- t(apply(testp,1,function(xx){(xx-median(xx))/mad(xx)}))

nclusters <- 64L
nalpha2 <- 5
minalpha <- 2^-3
maxalpha <- 2^3
nshape2 <- 5
minshape <- 1
maxshape <- 1
##
npoints <- ncol(testp)
nvar <- nrow(testp)
nalpha <- 2*nalpha2+1
alpha0 <- 2^(minalpha:maxalpha)
nshape <- 2*nshape2+1
shape0 <- 2^(minshape:maxshape)
shapehi0 <- 1
shapelo0 <- 1


datapoints <- list(datapoints = mtestp)
##
constants <- list(npoints=npoints, nvar=nvar, nclusters=nclusters, nalpha=nalpha, nshape=nshape)
##
initsFunction <- function(){
    probalpha0 <- rep(1/nclusters, nclusters)
    Alpha <- runif(1, minalpha, maxalpha)
    walpha0 <- probalpha0 * Alpha
    W <- rdirch(n=1, alpha=walpha0)
    ##
    ## Rmean1 <- rnorm(n=nvar, mean=0, sd=1)
    Rmean1 <- rep(0, nvar)
    Rratem1 <- rinvgamma(n=nvar, shape=shapehi0, rate=1)
    Rvarm1 <- rinvgamma(n=nvar, shape=shapelo0, rate=Rratem1)
    ##
    ## Rrate1 <- rinvgamma(n=nvar, shape=shapehi0, rate=1)
    ## Rvar1 <- rinvgamma(n=nvar, shape=shapelo0, rate=Rrate1)
    Rvar1 <- rep(1, nvar)
    probshape0 <- numeric(nshape)
    ##probshape0[(minshape:maxshape)+nshape2+1] <- rep(1/(maxshape-minshape+1),maxshape-minshape+1)
    ## wshape0 <- 2^((-nshape2):nshape2)
    Shapeindexlo <- 1+0*runif(nvar, minshape, maxshape)
    Shapeindexhi <- 1+0*runif(nvar, minshape, maxshape)
    ##
    Rmean <- matrix(rnorm(n=nvar*nclusters, mean=Rmean1, sd=sqrt(Rvarm1)), nrow=nvar, ncol=nclusters)
    Rrate <- matrix(rinvgamma(n=nvar*nclusters, shape=2^Shapeindexhi, rate=Rvar1), nrow=nvar, ncol=nclusters)
    Rvar <- matrix(rinvgamma(n=nvar*nclusters, shape=2^Shapeindexlo, rate=Rrate), nrow=nvar, ncol=nclusters)
    ## Rrate <- matrix(rinvgamma(n=nvar*nclusters, shape=shapehi0, rate=Rvar1), nrow=nvar, ncol=nclusters)
    ## Rvar <- matrix(rinvgamma(n=nvar*nclusters, shape=shapelo0, rate=Rrate), nrow=nvar, ncol=nclusters)
    K <- rep(1, npoints)
    list(
        minalpha = minalpha,
        maxalpha = maxalpha,
        probalpha0 = probalpha0,
        walpha0 = walpha0,
        Alpha = 2^0,
        W = W,
        K = K,
        ##
        Rmean1 = Rmean1,
        Rvarm1 = Rvarm1,
        Rratem1 = Rratem1,
        ## Rrate1 = Rrate1,
        Rvar1 = Rvar1,
        Rmean = Rmean,
        Rrate = Rrate,
        Rvar = Rvar,
        ##
        minshape = minshape,
        maxshape = maxshape,
        ##wshape0 = wshape0,
        probshape0 = probshape0,
        Shapeindexlo = Shapeindexlo,
        Shapeindexhi = Shapeindexhi,
        shapehi0 = shapehi0,
        shapelo0 = shapelo0
    )
}

finitemix <- nimbleCode({
    Alpha ~ dunif(minalpha, maxalpha)
    walpha0[1:nclusters] <- probalpha0[1:nclusters] * Alpha
    W[1:nclusters] ~ ddirch(alpha=walpha0[1:nclusters])
    ##
    for(v in 1:nvar){
        ## Rmean1[v] ~ dnorm(mean=0, var=1)
        Rratem1[v] ~ dinvgamma(shape=shapehi0, rate=1)
        Rvarm1[v] ~ dinvgamma(shape=shapelo0, rate=Rratem1[v])
        ##
        Shapeindexlo[v] ~ dinvgamma(shape=minshape, scale=maxshape)
        Shapeindexhi[v] ~ dinvgamma(shape=minshape, scale=maxshape)
        ## Rrate1[v] ~ dinvgamma(shape=shapehi0, rate=1)
        ## Rvar1[v] ~ dinvgamma(shape=shapelo0, rate=Rrate1[v])
    }
    for(k in 1:nclusters){
        for(v in 1:nvar){
            Rmean[v, k] ~ dnorm(mean=Rmean1[v], var=Rvarm1[v])
            Rrate[v, k] ~ dinvgamma(shape=2^Shapeindexhi[v], rate=Rvar1[v])
            Rvar[v, k] ~ dinvgamma(shape=2^Shapeindexlo[v], rate=Rrate[v, k])
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
                            monitors=c('Alpha', 'K',  'W', 'Rmean', 'Rvar',
                                       'Rvarm1',
                                       ##'Rrate',
                                       'Shapeindexhi',
                                       'Shapeindexlo',
                                       'Rvar1'
                                       )
                            )
confnimble$removeSamplers(c('Alpha', 'Shapeindexlo','Shapeindexhi'))
for(no in c('Alpha','Shapeindexhi[1]','Shapeindexhi[2]','Shapeindexlo[1]','Shapeindexlo[2]')){confnimble$addSampler(target=no,type='slice')}
print(confnimble)

## confnimble$printSamplers(executionOrder=TRUE)
## confnimble$getSamplerExecutionOrder()
confnimble$setSamplerExecutionOrder(rev(confnimble$getSamplerExecutionOrder()))

mcsampler <- buildMCMC(confnimble)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

nclusters <- 64L
minalpha <- 2^-3
maxalpha <- 2^3
minshape <- 1
maxshape <- 1
shapehi0 <- 1
shapelo0 <- 1
##
alpha0 <- 2^(minalpha:maxalpha)
walpha0 <- 2^((-nalpha2):nalpha2)
shape0 <- 2^(minshape:maxshape)
wshape0 <- 2^((-nshape2):nshape2)
##
set.seed(890)
Cfinitemixnimble$setInits(initsFunction())
thistime <- Sys.time()
todelete <- Cmcsampler$run(niter=1, thin=1, thin2=1, nburnin=0, time=TRUE, reset=TRUE, resetMV=TRUE)
todelete <- Cmcsampler$run(niter=15000*3, thin=150, thin2=1, nburnin=0, time=TRUE, reset=FALSE, resetMV=FALSE)
thistime <- Sys.time()-thistime
print(thistime)
mcsamples <- as.matrix(Cmcsampler$mvSamples)
##
## tplot(y=log2(mcsamples[,'Alpha']))
##
samplertimes <- Cmcsampler$getTimes()
names(samplertimes) <- sapply(confnimble$getSamplers(),function(x)x$target)
sprefixes <- unique(sub('^([^[]+)(\\[.*\\])', '\\1', names(samplertimes)))
cat(paste0('\nSampler times:\n'))
print(sort(sapply(sprefixes, function(x)sum(samplertimes[grepl(x,names(samplertimes))])),decreasing=T))
## > Sampler times:
## >            K         Rvar        Rmean Shapeindexhi Shapeindexlo        Rrate 
##    57.624822    31.890464    25.770928     4.600648     3.966075     3.226554 
##            W        Alpha       Rvarm1      Rratem1 
##     1.358030     0.629642     0.243421     0.037783 
## Sampler times:
## >            K         Rvar        Rmean        Rrate            W Shapeindexlo 
##    59.860321    36.723186    27.529190     1.881085     1.390121     1.382306 
## Shapeindexhi   Alphaindex       Rvarm1      Rratem1 
##     1.322861     0.400863     0.181544     0.024994 
##
## Sampler times:
## >            K         Rvar        Rmean        Rrate            W Shapeindexhi 
##    57.077974    34.620131    30.770944     1.969160     1.319746     0.477963 
## Shapeindexlo       Rvarm1        Alpha      Rratem1 
##     0.425648     0.173787     0.076777     0.026015 
## print('alpha')
## mean(mcsamples[,extract('Alpha')])
## print('shapehi')
## mean(mcsamples[,extract('Shapeindexhi')])
## print('shapelo')
## mean(mcsamples[,extract('Shapeindexlo')])
## print('Rvar1')
## (mean(log10(mcsamples[,extract('Rvar1')])/2))
## #tplot(y=list(log10(mcsamples[,'Rvar1[1]'])/2,log10(mcsamples[,'Rvar1[2]'])/2))
## print('Rvarm1')
## (mean(log10(mcsamples[,extract('Rvarm1')])))
## #tplot(y=list(log10(mcsamples[,'Rvarm1[1]'])/2,log10(mcsamples[,'Rvarm1[2]'])/2))

## tplot(y=apply(mcsamples,1,function(rr){length(unique(rr[extract('K')]))}))
## table(apply(mcsamples,1,function(rr){length(unique(rr[extract('K')]))}))


## tplot(y=(mcsamples[,extract('Shapeindexlo')]))

## tplot(y=(mcsamples[,extract('Shapeindexhi')]))

## test <- apply(mcsamples,1,function(ss){
##     probs <- sapply(1:11,function(ii){
##         ddirch(ss[extract('W')], alpha=rep(walpha0[ii]/nclusters, nclusters), log=T)
##     })
##     dcat(1:11, prob=exp(probs-max(probs)))
## })





## pdff('testCnimble')
## tplot(y=mcsamples[,'Alphaindex'], ylab='Alphaindex')
## dev.off()


#### Check resulting densities for various choices of hyperpriors - cont variate
incl <- 95
nsamples <- 2^12
rowcol <- c(24,34)
prc <- prod(rowcol)
drawf <- function(n, q, means, sds){
    ks <- sample(1:nclusters, size=n, prob=q, replace=T)
    cbind(
        rnorm(n,
              mean=means[1,ks],
              sd=sds[1,ks]
              )
       ,
        rnorm(n,
              mean=means[2,ks],
              sd=sds[2,ks]
              )
    )
}
plotpoints2d <- function(nsamples, q, means, sds){
    xl <-- numeric(prc)
    ## 2D
    for(i in 1:prc){
        points <- drawf(n=nsamples, q=q[i,], means=means[i,,], sds=sds[i,,])
        xl[i] <- max(abs(tquant(c(points[,1]), c((100-incl)/2,(100+incl)/2)/100)), 1)
        yl <- max(xl[i],abs(tquant(c(points[,2]), c((100-incl)/2,(100+incl)/2)/100)), 1)
        tplot(x=points[,1], y=points[,2], type='p', pch='.',
              alpha=0.9,
              xlab=NA, ylab=NA, xticks=NA, yticks=NA,
              xlabels=NA, ylabels=NA,
              xlim=c(-yl,yl), ylim=c(-yl,yl),
              mar=rep(0.2,4))
        abline(h=c(-1,1),col=alpha2hex2(0.5,2))
        abline(v=c(-1,1),col=alpha2hex2(0.5,2))
        abline(v=par('usr')[1:2],col='black',lwd=0.5)
        abline(h=par('usr')[3:4],col='black',lwd=0.5)
    }
    xl <<- xl
    ## 1D
    for(i in 1:prc){
        xgrid <- seq(-xl[i], xl[i], length.out=256)
        ygrid <- c(sapply(1:nclusters, function(cc){
            dnorm(xgrid, mean=means[i,1,cc], sd=sds[i,1,cc])
            }) %*% q[i,])
        tplot(x=xgrid, y=ygrid, lwd=1,
              xlab=NA, ylab=NA, xticks=NA, yticks=NA,
              xlabels=NA, ylabels=NA,
              ylim=c(0,NA),
              mar=rep(0.25,4))
        abline(h=c(0),col=alpha2hex2(0.5,7))
        abline(v=c(-1,1),col=alpha2hex2(0.5,2))
        ## abline(v=par('usr')[1:2],col='black',lwd=0.5)
        ## abline(h=par('usr')[3:4],col='black',lwd=0.5)
    }
}
##
set.seed(111)
alphas <- runif(prc, minalpha, maxalpha)
q <- extraDistr::rdirichlet(n=prc,alpha=matrix(alphas/nclusters,nrow=prc,ncol=nclusters))
##
##meansm <- rnorm(prc*2, mean=0, sd=1)
meanss <- sqrt(nimble::rinvgamma(prc*2, shape=1, rate= nimble::rinvgamma(prc*2, shape=1, rate=1)))
sdssl <- rinvgamma(prc*2, shape=minshape, scale=maxshape)
sdssh <- rinvgamma(prc*2, shape=minshape, scale=maxshape)
sdsr <- 1#nimble::rinvgamma(prc*2, shape=baseshape1, rate= nimble::rinvgamma(prc*2, shape=baseshape1, rate=1))
##
set.seed(987)
means <- array(rnorm(2*prc*nclusters, mean=0, sd=meanss), dim=c(prc,2,nclusters))
sds <- array(sqrt(nimble::rinvgamma(prc*2*nclusters, shape=sdssl, rate=
                                                           nimble::rinvgamma(prc*2*nclusters, shape=sdssh, rate=sdsr))), dim=c(prc,2,nclusters))
##
pdff(paste0('_xigsprior2D_Mvar_S2shapes_A',minalpha,'_',maxalpha,'_S',minshape,'_',maxshape,'_N',npoints),apaper=3)
##
tplot(y=log2(mcsamples[,extract('Alpha')]), main=paste0('lb-alpha ',mean(mcsamples[,extract('Alpha')])), xlab=NA, ylab=NA)
tplot(y=apply(mcsamples,1,function(rr){length(unique(rr[extract('K')]))}), main=paste0('K  ',mean(apply(mcsamples,1,function(rr){length(unique(rr[extract('K')]))}))), xlab=NA, ylab=NA)
tplot(y=log2(mcsamples[,extract('Shapeindexlo')]), main=paste0('lb-shape-lo ',mean(mcsamples[,extract('Shapeindexlo')])), xlab=NA, ylab=NA)
if(length(extract('Shapeindexhi'))){
    tplot(y=log2(mcsamples[,extract('Shapeindexhi')]), main=paste0('lb-shape-hi ',mean(mcsamples[,extract('Shapeindexhi')])), xlab=NA, ylab=NA)
}
tplot(y=0.5*log10(mcsamples[,extract('Rvar1')]), main=paste0('hlg-Rvar1 ',mean(0.5*log10(mcsamples[,extract('Rvar1')]))), xlab=NA, ylab=NA)
tplot(y=0.5*log10(mcsamples[,extract('Rvarm1')]), main=paste0('hlg-Rvarm1 ',mean(0.5*log10(mcsamples[,extract('Rvarm1')]))), xlab=NA, ylab=NA)
##
par(mfrow=rowcol,mar = c(0,0,0,0))
## 
plotpoints2d(nsamples=nsamples, q=q, means=means, sds=sds)
##dev.off()
##
draws <- function(n,i){
    kk <- rcat(n=n, prob=mcsamples[i, extract('W')])
    means <- t(matrix(mcsamples[i, extract('Rmean')], nrow=nvar, ncol=nclusters)[,kk,drop=F])
    sds <- sqrt(t(matrix(mcsamples[i, extract('Rvar')], nrow=nvar, ncol=nclusters)[,kk,drop=F]))
    t(matrix(rnorm(n=2*n, mean=means, sd=sds), nrow=n, ncol=nvar))
}
##
##pdff(paste0('testmchyper-mean_var-sd_2shapes'),apaper=3)
##
par(mfrow=c(1,1))
for(i in round(seq(nrow(mcsamples),1,length.out=10))){
testpoints <- draws(1e5,i)*locdis[,2]+locdis[,1]
xm <- tquant(testpoints[1,],c(2.5,97.5)/100)
ym <- tquant(testpoints[2,],c(2.5,97.5)/100)
tplot(x=list(testpoints[1,],testp[1,]), y=list(testpoints[2,],testp[2,]), type='p', pch=c(46,16), cex=c(1,1), alpha=c(0.9,0.75),
      xlim=xm, ylim=xm, main=thistime)
}
dev.off()
