mystartup()
stopCluster(cluster)
stopImplicitCluster()
registerDoSEQ()

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
minshape <- 0.5
maxshape <- 0.5
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
    Rvarm1 <- 1+0*rinvgamma(n=nvar, shape=shapelo0, rate=Rratem1)
    ##
    ## Rrate1 <- rinvgamma(n=nvar, shape=shapehi0, rate=1)
    ## Rvar1 <- rinvgamma(n=nvar, shape=shapelo0, rate=Rrate1)
    Rvar1 <- rep(1, nvar)
    probshape0 <- numeric(nshape)
    ##probshape0[(minshape:maxshape)+nshape2+1] <- rep(1/(maxshape-minshape+1),maxshape-minshape+1)
    ## wshape0 <- 2^((-nshape2):nshape2)
    Shapelo <- minshape + 0*rinvgamma(nvar, shape=minshape, scale=maxshape)
    Shapehi <- maxshape + 0*rinvgamma(nvar, shape=minshape, scale=maxshape)
    ##
    Rmean <- matrix(rnorm(n=nvar*nclusters, mean=Rmean1, sd=sqrt(Rvarm1)), nrow=nvar, ncol=nclusters)
    Rrate <- matrix(rinvgamma(n=nvar*nclusters, shape=Shapehi, rate=Rvar1), nrow=nvar, ncol=nclusters)
    Rvar <- matrix(rinvgamma(n=nvar*nclusters, shape=Shapelo, rate=Rrate), nrow=nvar, ncol=nclusters)
    ## Rrate <- matrix(rinvgamma(n=nvar*nclusters, shape=shapehi0, rate=Rvar1), nrow=nvar, ncol=nclusters)
    ## Rvar <- matrix(rinvgamma(n=nvar*nclusters, shape=shapelo0, rate=Rrate), nrow=nvar, ncol=nclusters)
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
if(FALSE){
    for(v in 1:nvar){
        ## Rmean1[v] ~ dnorm(mean=0, var=1)
        ## Rratem1[v] ~ dinvgamma(shape=shapehi0, rate=1)
        Rvarm1[v] ~ dinvgamma(shape=shapelo0, rate=Rratem1[v])
        ##
        Shapelo[v] ~ dinvgamma(shape=minshape, scale=maxshape)
        Shapehi[v] ~ dinvgamma(shape=minshape, scale=maxshape)
        ## Rrate1[v] ~ dinvgamma(shape=shapehi0, rate=1)
        ## Rvar1[v] ~ dinvgamma(shape=shapelo0, rate=Rrate1[v])
    }
    }
    for(k in 1:nclusters){
        for(v in 1:nvar){
            Rmean[v, k] ~ dnorm(mean=Rmean1[v], var=Rvarm1[v])
            Rrate[v, k] ~ dinvgamma(shape=Shapehi[v], rate=Rvar1[v])
            Rvar[v, k] ~ dinvgamma(shape=Shapelo[v], rate=Rrate[v, k])
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
minshape <- 1
maxshape <- 1
shapehi0 <- 1
shapelo0 <- 1
##
alpha0 <- 2^(minalpha:maxalpha)
walpha0 <- 2^((-nalpha2):nalpha2)
shape0 <- 2^(minshape:maxshape)
##wshape0 <- 2^((-nshape2):nshape2)
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
stargets <- sapply(confnimble$getSamplers(), function(xx)xx$target)   
confnimble$removeSamplers(c('Alpha', 'Shapelo','Shapehi'))
##
##    
for(no in intersect(c('Alpha','Shapehi[1]','Shapehi[2]','Shapelo[1]','Shapelo[2]'), stargets)){confnimble$addSampler(target=no,type='slice')}
print(confnimble)
##
## confnimble$printSamplers(executionOrder=TRUE)
    ## confnimble$getSamplerExecutionOrder()
## confnimble$setSamplerExecutionOrder(rev(confnimble$getSamplerExecutionOrder()))
    ## sampledv <- unique(sub('^(.+)\\[.*\\]', '\\1', stargets)
##
sorder <- c('K','Rmean','Rrate','Rvar','W','Alpha')    
neworder <- foreach(var=sorder, .combine=c)%do%{grep(paste0('^',var,'(\\[.+\\])*$'),stargets)}
confnimble$setSamplerExecutionOrder(neworder)
print(confnimble)
##
mcsampler <- buildMCMC(confnimble)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)
    ##
    gc()
    ##
set.seed(thisseed)
Cfinitemixnimble$setInits(initsFunction())
todelete <- Cmcsampler$run(niter=1, thin=1, thin2=1, nburnin=0, time=FALSE, reset=TRUE, resetMV=TRUE)
todelete <- Cmcsampler$run(niter=10000, thin=10000, thin2=1, nburnin=0, time=FALSE, reset=FALSE, resetMV=FALSE)
todelete <- Cmcsampler$run(niter=100*100, thin=100, thin2=1, nburnin=0, time=TRUE, reset=FALSE, resetMV=FALSE)
    ##
    samplertimes <- Cmcsampler$getTimes()
names(samplertimes) <- sapply(confnimble$getSamplers(),function(x)x$target)
sprefixes <- unique(sub('^([^[]+)(\\[.*\\])', '\\1', names(samplertimes)))
cat(paste0('\nSampler times:\n'))
print(sort(sapply(sprefixes, function(x)sum(samplertimes[grepl(x,names(samplertimes))])),decreasing=T))
##
    t(as.matrix(Cmcsampler$mvSamples))
    }
thistime <- Sys.time()-thistime
print(thistime)
mcsamples <- t(matrix(mcsamples, nrow=nrow(mcsamples)/ncores, dimnames=list(rownames(mcsamples)[1:(nrow(mcsamples)/ncores)],NULL)))

stopCluster(cluster)
stopImplicitCluster()
registerDoSEQ()
gc()


## Sampler times:
##         K      Rvar     Rmean     Rrate         W     Alpha 
## 38.344443 23.997377 19.620460  2.459773  1.144812  0.448341 
## 38.383626 24.228465 20.465167  2.471785  1.142326  0.454622 
## 38.058255 23.675885 19.431176  2.508051  1.140872  0.465032 
## 38.064739 23.793607 19.535739  2.443887  1.126356  0.454095 
## 37.74106 23.71560 20.11892  2.41589  1.13146  0.44499 
## 36.838983 22.763314 18.839487  2.385294  1.098117  0.441332 
## 41.62392 27.05462 22.44042  2.50587  1.20309  0.49493 
## 37.112570 23.507219 19.079654  2.272759  1.091747  0.433269 
## 36.743545 23.706603 19.843510  2.159650  1.077770  0.429274 
## 33.814611 21.131176 17.960721  2.023431  1.000664  0.405883 
## > > Time difference of 4.0236 mins


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
pdfname <- paste0('_tnprev-2D_Mn_Sn_A',minalpha,'_',maxalpha,'_S',minshape,'_',maxshape,'_N',npoints)
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
sdssl <- minshape + 0*nimble::rinvgamma(prc*2, shape=minshape, scale=maxshape)
sdssh <- maxshape + 0*nimble::rinvgamma(prc*2, shape=minshape, scale=maxshape)
sdsr <- 1#nimble::rinvgamma(prc*2, shape=baseshape1, rate= nimble::rinvgamma(prc*2, shape=baseshape1, rate=1))
##
set.seed(987)
means <- array(rnorm(2*prc*nclusters, mean=0, sd=1), dim=c(prc,2,nclusters))
sds <- array(sqrt(nimble::rinvgamma(prc*2*nclusters, shape=sdssl, rate=
                                                           nimble::rinvgamma(prc*2*nclusters, shape=sdssh, rate=sdsr))), dim=c(prc,2,nclusters))
##
pdff(pdfname,apaper=3)
##
tplot(y=log2(mcsamples[,extract('Alpha')]), main=paste0('lb-alpha ',mean(mcsamples[,extract('Alpha')])), xlab=NA, ylab=NA)
tplot(y=apply(mcsamples,1,function(rr){length(unique(rr[extract('K')]))}), main=paste0('K  ',mean(apply(mcsamples,1,function(rr){length(unique(rr[extract('K')]))}))), xlab=NA, ylab=NA, ylim=c(0,NA))
tplot(y=log2(mcsamples[,extract('Shapelo')]), main=paste0('lb-shape-lo ',mean(mcsamples[,extract('Shapelo')])), xlab=NA, ylab=NA)
if(length(extract('Shapehi'))){
    tplot(y=log2(mcsamples[,extract('Shapehi')]), main=paste0('lb-shape-hi ',mean(mcsamples[,extract('Shapehi')])), xlab=NA, ylab=NA)
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




## nsam <- 1e4
## testn <- 10
## tests <- sapply(1:nsam, function(xx){
##     mad(rt(n=testn, df=3))
## })
## testt <- mad(rt(n=1e6, df=3))
## 100*mad(tests-testt)/testt
## thist(tests,plot=T)
## abline(v=testt, col=2)

