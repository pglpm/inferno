library('data.table')
library('png')
library('foreach')
library('doRNG')
## library('doFuture')
## registerDoFuture()
## cat('\navailableCores: ')
## cat(availableCores())
## cat('\navailableCores-multicore: ')
## cat(availableCores('multicore'),'\n')
## if(Sys.info()['nodename']=='luca-HP-Z2-G9'){
##     ncores <- 20}else{
##     ncores <- 6}
## cat(paste0('\nusing ',ncores,' cores\n'))
## if(ncores>1){
##     if(.Platform$OS.type=='unix'){
##         plan(multicore, workers=ncores)
##     }else{
##         plan(multisession, workers=ncores)
##     }
## }else{
##     plan(sequential)
## }
##library('nimble')
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
minshape <- -1
maxshape <- 1
##
npoints <- ncol(testp)
nvar <- nrow(testp)
nalpha <- 2*nalpha2+1
alpha0 <- 2^(minalpha:maxalpha)
nshape <- length(minshape:maxshape)
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
    W <- nimble::rdirch(n=1, alpha=walpha0)
    ##
    ## Rmean1 <- rnorm(n=nvar, mean=0, sd=1)
    Rmean1 <- rep(0, nvar)
    Rratem1 <- nimble::rinvgamma(n=nvar, shape=shapehi0, rate=1)
    Rvarm1 <- 1+0*nimble::rinvgamma(n=nvar, shape=shapelo0, rate=Rratem1)
    ##
    ## Rrate1 <- nimble::rinvgamma(n=nvar, shape=shapehi0, rate=1)
    ## Rvar1 <- nimble::rinvgamma(n=nvar, shape=shapelo0, rate=Rrate1)
    Rvar1 <- rep(1, nvar)
    probshape0 <- rep(1/nshape, nshape)
    ##probshape0[(minshape:maxshape)+nshape2+1] <- rep(1/(maxshape-minshape+1),maxshape-minshape+1)
    ## wshape0 <- 2^((-nshape2):nshape2)
    Shapelo <- sample(1:nshape, nvar, replace=T)
    Shapehi <- sample(1:nshape, nvar, replace=T)
    ##
    Rmean <- matrix(rnorm(n=nvar*nclusters, mean=Rmean1, sd=sqrt(Rvarm1)), nrow=nvar, ncol=nclusters)
    Rrate <- matrix(nimble::rinvgamma(n=nvar*nclusters, shape=2^(Shapehi+minshape-1), rate=Rvar1), nrow=nvar, ncol=nclusters)
    Rvar <- matrix(nimble::rinvgamma(n=nvar*nclusters, shape=2^(Shapelo+minshape-1), rate=Rrate), nrow=nvar, ncol=nclusters)
    ## Rrate <- matrix(nimble::rinvgamma(n=nvar*nclusters, shape=shapehi0, rate=Rvar1), nrow=nvar, ncol=nclusters)
    ## Rvar <- matrix(nimble::rinvgamma(n=nvar*nclusters, shape=shapelo0, rate=Rrate), nrow=nvar, ncol=nclusters)
    K <- rep(1, npoints)
    list(
        minalpha = minalpha,
        maxalpha = maxalpha,
        probalpha0 = probalpha0,
        walpha0 = walpha0,
        Alpha = Alpha,
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
        Shapelo = Shapelo,
        Shapehi = Shapehi,
        shapehi0 = shapehi0,
        shapelo0 = shapelo0
    )
}
##
    finitemix <- nimble::nimbleCode({
    Alpha ~ dunif(minalpha, maxalpha)
    walpha0[1:nclusters] <- probalpha0[1:nclusters] * Alpha
    W[1:nclusters] ~ ddirch(alpha=walpha0[1:nclusters])
    ##
if(TRUE){
    for(v in 1:nvar){
        ## Rmean1[v] ~ dnorm(mean=0, var=1)
        ## Rratem1[v] ~ dinvgamma(shape=shapehi0, rate=1)
        ##Rvarm1[v] ~ dinvgamma(shape=shapelo0, rate=Rratem1[v])
        ##
        Shapelo[v] ~ dcat(prob=probshape0[1:nshape])
        Shapehi[v] ~ dcat(prob=probshape0[1:nshape])
        ## Rrate1[v] ~ dinvgamma(shape=shapehi0, rate=1)
        ## Rvar1[v] ~ dinvgamma(shape=shapelo0, rate=Rrate1[v])
    }
    }
    for(k in 1:nclusters){
        for(v in 1:nvar){
            Rmean[v, k] ~ dnorm(mean=Rmean1[v], var=Rvarm1[v])
            Rrate[v, k] ~ dinvgamma(shape=2^(Shapehi[v]+minshape-1), rate=Rvar1[v])
            Rvar[v, k] ~ dinvgamma(shape=2^(Shapelo[v]+minshape-1), rate=Rrate[v, k])
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
##
nclusters <- 64L
minalpha <- 2^-3
maxalpha <- 2^3
minshape <- -1
maxshape <- 1
shapehi0 <- 1
shapelo0 <- 1
##
alpha0 <- 2^(minalpha:maxalpha)
walpha0 <- 2^((-nalpha2):nalpha2)
shape0 <- 2^(minshape:maxshape)
##

ncores <- 10
stopCluster(cluster)
stopImplicitCluster()
registerDoSEQ()
## cl <- makePSOCKcluster(ncores)
cluster <- makeCluster(ncores, outfile='')
registerDoParallel(cluster)
## if(ncores>1){
##     if(.Platform$OS.type=='unix'){
##         plan(multicore, workers=ncores)
##     }else{
##         plan(multisession, workers=ncores)
##     }
## }else{
##     plan(sequential)
## }
##

thistime <- Sys.time()
mcsamples <- foreach(cor=1:ncores, .combine=rbind, .packages='nimble', .inorder=FALSE)%dorng%{
    thisseed <- 890+cor
##  
finitemixnimble <- nimbleModel(code=finitemix, name='finitemixnimble1',
                               constants=constants,
                               inits=initsFunction(),
                               data=datapoints
                               )
##
Cfinitemixnimble <- compileNimble(finitemixnimble, showCompilerOutput=FALSE)
##
confnimble <- configureMCMC(Cfinitemixnimble,
                            monitors=c('Alpha', 'K',  'W', 'Rmean', 'Rvar',
                                       'Rvarm1',
                                       ##'Rrate',
                                       'Shapehi',
                                       'Shapelo',
                                       'Rvar1'
                                       )
                            )
confnimble$removeSamplers(c('Alpha', 'Shapelo','Shapehi'))
for(no in c('Alpha','Shapehi[1]','Shapehi[2]','Shapelo[1]','Shapelo[2]')){confnimble$addSampler(target=no,type='slice')}
print(confnimble)
##
## confnimble$printSamplers(executionOrder=TRUE)
## confnimble$getSamplerExecutionOrder()
confnimble$setSamplerExecutionOrder(rev(confnimble$getSamplerExecutionOrder()))
##
mcsampler <- buildMCMC(confnimble)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)
    ##
    gc()
    ##
set.seed(thisseed)
Cfinitemixnimble$setInits(initsFunction())
todelete <- Cmcsampler$run(niter=1, thin=1, thin2=1, nburnin=0, time=TRUE, reset=TRUE, resetMV=TRUE)
todelete <- Cmcsampler$run(niter=10000, thin=10000, thin2=1, nburnin=0, time=TRUE, reset=FALSE, resetMV=FALSE)
todelete <- Cmcsampler$run(niter=100*100, thin=100, thin2=1, nburnin=0, time=TRUE, reset=FALSE, resetMV=FALSE)
    t(as.matrix(Cmcsampler$mvSamples))
    }
thistime <- Sys.time()-thistime
print(thistime)
mcsamples <- t(matrix(mcsamples, nrow=nrow(mcsamples)/ncores, dimnames=list(rownames(mcsamples)[1:(nrow(mcsamples)/ncores)],NULL)))

stopCluster(cluster)
stopImplicitCluster()
registerDoSEQ()
gc()



##
## tplot(y=log2(mcsamples[,'Alpha']))
##
samplertimes <- Cmcsampler$getTimes()
names(samplertimes) <- sapply(confnimble$getSamplers(),function(x)x$target)
sprefixes <- unique(sub('^([^[]+)(\\[.*\\])', '\\1', names(samplertimes)))
cat(paste0('\nSampler times:\n'))
print(sort(sapply(sprefixes, function(x)sum(samplertimes[grepl(x,names(samplertimes))])),decreasing=T))
## > Sampler times:
## >            K         Rvar        Rmean Shapehi Shapelo        Rrate 
##    57.624822    31.890464    25.770928     4.600648     3.966075     3.226554 
##            W        Alpha       Rvarm1      Rratem1 
##     1.358030     0.629642     0.243421     0.037783 
## Sampler times:
## >            K         Rvar        Rmean        Rrate            W Shapelo 
##    59.860321    36.723186    27.529190     1.881085     1.390121     1.382306 
## Shapehi   Alphaindex       Rvarm1      Rratem1 
##     1.322861     0.400863     0.181544     0.024994 
##
## Sampler times:
## >            K         Rvar        Rmean        Rrate            W Shapehi 
##    57.077974    34.620131    30.770944     1.969160     1.319746     0.477963 
## Shapelo       Rvarm1        Alpha      Rratem1 
##     0.425648     0.173787     0.076777     0.026015 
## print('alpha')
## mean(mcsamples[,extract('Alpha')])
## print('shapehi')
## mean(mcsamples[,extract('Shapehi')])
## print('shapelo')
## mean(mcsamples[,extract('Shapelo')])
## print('Rvar1')
## (mean(log10(mcsamples[,extract('Rvar1')])/2))
## #tplot(y=list(log10(mcsamples[,'Rvar1[1]'])/2,log10(mcsamples[,'Rvar1[2]'])/2))
## print('Rvarm1')
## (mean(log10(mcsamples[,extract('Rvarm1')])))
## #tplot(y=list(log10(mcsamples[,'Rvarm1[1]'])/2,log10(mcsamples[,'Rvarm1[2]'])/2))

## tplot(y=apply(mcsamples,1,function(rr){length(unique(rr[extract('K')]))}))
## table(apply(mcsamples,1,function(rr){length(unique(rr[extract('K')]))}))


## tplot(y=(mcsamples[,extract('Shapelo')]))

## tplot(y=(mcsamples[,extract('Shapehi')]))

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
pdfname <- paste0('_tnp-2D_Mn_S2shapesDiscr_A',minalpha,'_',maxalpha,'_S',minshape,'_',maxshape,'_N',npoints)
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
sdssl <- 2^(sample(1:nshape,prc*2,replace=T)+minshape-1)
sdssh <- 2^(sample(1:nshape,prc*2,replace=T)+minshape-1)
sdsr <- 1#nimble::rinvgamma(prc*2, shape=baseshape1, rate= nimble::rinvgamma(prc*2, shape=baseshape1, rate=1))
##
set.seed(987)
means <- array(rnorm(2*prc*nclusters, mean=0, sd=meanss), dim=c(prc,2,nclusters))
sds <- array(sqrt(nimble::rinvgamma(prc*2*nclusters, shape=sdssl, rate=
                                                           nimble::rinvgamma(prc*2*nclusters, shape=sdssh, rate=sdsr))), dim=c(prc,2,nclusters))
##
pdff(pdfname,apaper=3)
##
tplot(y=log2(mcsamples[,extract('Alpha')]), main=paste0('lb-alpha ',mean(mcsamples[,extract('Alpha')])), xlab=NA, ylab=NA)
tplot(y=apply(mcsamples,1,function(rr){length(unique(rr[extract('K')]))}), main=paste0('K  ',mean(apply(mcsamples,1,function(rr){length(unique(rr[extract('K')]))}))), xlab=NA, ylab=NA, ylim=c(0,NA))
tplot(y=-1+minshape+(mcsamples[,extract('Shapelo')]), main=paste0('lb-shape-lo ',minshape-1+mean(mcsamples[,extract('Shapelo')])), xlab=NA, ylab=NA)
if(length(extract('Shapehi'))){
    tplot(y=-1+minshape+(mcsamples[,extract('Shapehi')]), main=paste0('lb-shape-hi ',minshape-1+mean(mcsamples[,extract('Shapehi')])), xlab=NA, ylab=NA)
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
    kk <- nimble::rcat(n=n, prob=mcsamples[i, extract('W')])
    means <- t(matrix(mcsamples[i, extract('Rmean')], nrow=nvar, ncol=nclusters)[,kk,drop=F])
    sds <- sqrt(t(matrix(mcsamples[i, extract('Rvar')], nrow=nvar, ncol=nclusters)[,kk,drop=F]))
    t(matrix(rnorm(n=2*n, mean=means, sd=sds), nrow=n, ncol=nvar))
}
##
##pdff(paste0('testmchyper-mean_var-sd_2shapes'),apaper=3)
##
par(mfrow=c(1,1))
for(i in round(seq(nrow(mcsamples),1,length.out=min(20,nrow(mcsamples))))){
testpoints <- draws(1e5,i)*locdis[,2]+locdis[,1]
xm <- tquant(testpoints[1,],c(2.5,97.5)/100)
ym <- tquant(testpoints[2,],c(2.5,97.5)/100)
tplot(x=list(testpoints[1,],testp[1,]), y=list(testpoints[2,],testp[2,]), type='p', pch=c(46,16), cex=c(1,1), alpha=c(0.9,0.75),
      xlim=xm, ylim=xm, main=thistime)
}
dev.off()




## nsam <- 1e6
## xgrid <- seq(-4,4,length.out=64)
## his1 <- thist(0.5*log10(nimble::rinvgamma(nsam, shape=1, rate=nimble::rinvgamma(nsam, shape=1, rate=1))), n=xgrid)
## his2 <- thist(0.5*log10(nimble::rinvgamma(nsam, shape=0.5, rate=nimble::rinvgamma(nsam, shape=0.5, rate=1))), n=xgrid)
## his3 <- thist(0.5*log10(nimble::rinvgamma(nsam, shape=2^sample((-1):1,nsam,replace=T), rate=nimble::rinvgamma(nsam,shape=2^sample((-1):1,nsam,replace=T), rate=1))), n=xgrid)
## ##
## tplot(x=his1$mids,y=lapply(list(his1$density,his2$density,his3$density),log10), ylim=c(NA,NA))

nsam <- 1e5
xgrid <- seq(-4,4,length.out=64)
sam2 <- matrix(0.5*log10(nimble::rinvgamma(nsam*2, shape=0.5, rate=nimble::rinvgamma(nsam*2, shape=0.5, rate=1))), ncol=2)
sam3 <- matrix(0.5*log10(nimble::rinvgamma(nsam*2, shape=2^sample((-1):1,nsam,replace=T), rate=nimble::rinvgamma(nsam*2, shape=2^sample((-1):1,nsam,replace=T), rate=1))), ncol=2)
xm <- max(abs(tquant(c(sam2,sam3),c(0.5,99.5)/100)))
tplot(x=list(sam2[,1], sam3[,1]), y=list(sam2[,2], sam3[,2]), type='p', cex=0.25, alpha=0.9, xlim=c(-1,1)*xm, ylim=c(-1,1)*xm)


his1 <- thist(0.5*log10(nimble::rinvgamma(nsam, shape=1, rate=nimble::rinvgamma(nsam, shape=1, rate=1))), n=xgrid)
his2 <- thist(0.5*log10(nimble::rinvgamma(nsam, shape=0.5, rate=nimble::rinvgamma(nsam, shape=0.5, rate=1))), n=xgrid)
his3 <- thist(0.5*log10(nimble::rinvgamma(nsam, shape=2^sample((-1):1,nsam,replace=T), rate=nimble::rinvgamma(nsam,shape=2^sample((-1):1,nsam,replace=T), rate=1))), n=xgrid)
##
tplot(x=his1$mids,y=lapply(list(his1$density,his2$density,his3$density),log10), ylim=c(NA,NA))
