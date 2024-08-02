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
    
set.seed(123)
nn <- 32

## 6.54423 is the largest value produced by rnorm that I have seen
## 1e-3 -> 2^-10
## 1e-4 -> 2^-14
## 1e-5 -> 2^-17
## 1e-6 -> 2^-20
dx <- 1e-3
ddp <- seq(-32,-8,length.out=128)
dd <- 2^ddp
tran <- function(x,dd){qlogis(x*(1-2*dd)+dd)}
dtran <- function(x,dd){(1-2*dd)/dlogis(plogis(tran(x,dd)))}
## tran <- function(x,dd){qlogis(x*(1-2*dd)+dd)}
## dtran <- function(x,dd){(1-2*dd)/dlogis(plogis(tran(x,dd)))}
tx <- tran(0,dd)
dtx <- dtran(0,dd)
ygrid <- dnorm(tx)*dtx*dx/pnorm(tx)
## dtx <- dtran(dd,dd)-dtran(0,dd)
## ygrid <- dtx/pnorm(tx)dx
tplot(x=ddp,y=1/ygrid)

ddp <- seq(1,15,length.out=128)
tplot(x=ddp,y=qnorm(-ddp*log(10),log.p=T))

ddp <- seq(1,3,length.out=128)
tplot(x=ddp,y=dnorm(-10^ddp)==0)


set.seed(222)
#### Integer
#### with norm transformation IIb
nint <- 16
nmin <- 9
nmax <- nmin-1+nint
scale <- 1#-qnorm(0.5/nint)
lx <- nmin - 0.5*(nmax-nmin)/(nint-1)
sx <- nint*(nmax-nmin)/(nint-1)
tran <- function(x){qnorm((x-lx)/sx)/scale}
invtran <- function(y){round((sx*pnorm(y)+lx-nmin)*(nmax-nmin)/(nint-1))*(nmax-nmin)/(nint-1)+nmin}
nsamples <- 400*8
nsubsamples <- 400
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(1)/scale
shape1s <- c(1) # large scales
shape2s <- c(1) # small scales
scales <- (1/8/scale)^(-2)
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scaleprec <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scaleprec))),nsamples)
##
graphics.off()
pdff('samples_integer_normIIb')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(nmin,nmax,length.out=nint)
extr <- c(1,length(xgrid))
mgrid <- (xgrid[-extr[2]]+xgrid[-extr[1]])/2
mextr <- c(1,length(mgrid))
txgrid <- tran(xgrid)
tmgrid <- tran(mgrid)
##tmgrid <- (txgrid[-extr[2]]+txgrid[-extr[1]])/2
dx <- 2/(nint-1)
ysum <- 0
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- c(
            pnorm(tmgrid[mextr[1]], m[i,acluster], s[i,acluster]),
            pnorm(tmgrid[-mextr[1]], m[i,acluster], s[i,acluster]) - pnorm(tmgrid[-mextr[2]], m[i,acluster], s[i,acluster]),
            pnorm(tmgrid[mextr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
        )
        q[i,acluster]*dens}))
    ## y2 <- rowSums(sapply(1:nclusters,function(acluster){
    ##     dens <- dnorm(xgrid, m[i,acluster], s[i,acluster])*2*dx
    ##     dens[extr[1]] <- pnorm(xgrid[extr[1]], m[i,acluster], s[i,acluster])
    ##     dens[extr[2]] <- pnorm(xgrid[extr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
    ##     q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i<nsubsamples|i==nsamples){
    if(i==nsamples){y <- ysum/nsamples}
    #y2 <- y2*max(y[-extr])/max(y2[-extr])
    ## y[extr] <- y[extr] * max(y[-extr])
    tplot(x=xgrid, y=y, #y=list(y,y2),
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=c(if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}},4), ly=1,lwd=0.5)
    tplot(x=xgrid, y=y,type='p',cex=0.2,
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,add=T,
          col=if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}}, ly=1,lwd=0.5)
    tplot(x=xgrid[extr], y=y[extr],
          type='p',cex=0.075,col=3,add=T)
    ## tplot(x=xgrid[extr], y=y[extr], type='p',cex=0.1,add=T)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    if(i==nsamples){
        abline(h=c(1/nint),lwd=0.5,col=alpha2hex(0.5,c(2)),lty=1)
    }
    ## abline(v=,lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
}
dev.off()


set.seed(222)
#### Integer
#### with norm transformation IV
nint <- 16
nmin <- 9
nmax <- nmin-1+nint
scale <- 1#-qnorm(0.5/nint)
lx <- nmin - 0.5*(nmax-nmin)/(nint-1)
sx <- nint*(nmax-nmin)/(nint-1)
tran <- function(x){qnorm((x-lx)/sx)/scale}
invtran <- function(y){round((sx*pnorm(y)+lx-nmin)*(nmax-nmin)/(nint-1))*(nmax-nmin)/(nint-1)+nmin}
nsamples <- 400*8
nsubsamples <- 400
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(1)/scale
shape1s <- c(1) # large scales
shape2s <- c(1) # small scales
scales <- (1/8/scale)^(-2)
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scaleprec <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scaleprec))),nsamples)
##
graphics.off()
pdff('samples_integer_normIIb')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(nmin,nmax,length.out=nint)
extr <- c(1,length(xgrid))
mgrid <- (xgrid[-extr[2]]+xgrid[-extr[1]])/2
mextr <- c(1,length(mgrid))
txgrid <- tran(xgrid)
tmgrid <- tran(mgrid)
##tmgrid <- (txgrid[-extr[2]]+txgrid[-extr[1]])/2
dx <- 2/(nint-1)
ysum <- 0
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- c(
            pnorm(tmgrid[mextr[1]], m[i,acluster], s[i,acluster]),
            pnorm(tmgrid[-mextr[1]], m[i,acluster], s[i,acluster]) - pnorm(tmgrid[-mextr[2]], m[i,acluster], s[i,acluster]),
            pnorm(tmgrid[mextr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
        )
        q[i,acluster]*dens}))
    ## y2 <- rowSums(sapply(1:nclusters,function(acluster){
    ##     dens <- dnorm(xgrid, m[i,acluster], s[i,acluster])*2*dx
    ##     dens[extr[1]] <- pnorm(xgrid[extr[1]], m[i,acluster], s[i,acluster])
    ##     dens[extr[2]] <- pnorm(xgrid[extr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
    ##     q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i<nsubsamples|i==nsamples){
    if(i==nsamples){y <- ysum/nsamples}
    #y2 <- y2*max(y[-extr])/max(y2[-extr])
    ## y[extr] <- y[extr] * max(y[-extr])
    tplot(x=xgrid, y=y, #y=list(y,y2),
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=c(if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}},4), ly=1,lwd=0.5)
    tplot(x=xgrid, y=y,type='p',cex=0.2,
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,add=T,
          col=if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}}, ly=1,lwd=0.5)
    tplot(x=xgrid[extr], y=y[extr],
          type='p',cex=0.075,col=3,add=T)
    ## tplot(x=xgrid[extr], y=y[extr], type='p',cex=0.1,add=T)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    if(i==nsamples){
        abline(h=c(1/nint),lwd=0.5,col=alpha2hex(0.5,c(2)),lty=1)
    }
    ## abline(v=,lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
}
dev.off()


## set.seed(222)
sd2iqr <- 0.5/qnorm(0.75)
#### Strictly positive
#### with log transformation
scale <- 0.14*sd2iqr#-qnorm(0.5/nint)
location <- 4.3
tran <- function(x){(log(x)-location)/scale}
invtran <- function(y){exp(y*scale+location)}
jac <- function(x){1/(scale*x)}
nsamples <- 400*8
nsubsamples <- 400
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(2)
shape1s <- c(1/2) # large scales
shape2s <- c(1/2) # small scales
scales <- (1/4)^(-2)
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scaleprec <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scaleprec))),nsamples)
##
graphics.off()
pdff('samples_positive_log')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(invtran(-3),invtran(3),length.out=256)
extr <- c(1,length(xgrid))
txgrid <- tran(xgrid)
##tmgrid <- (txgrid[-extr[2]]+txgrid[-extr[1]])/2
dx <- 2/(nint-1)
ysum <- 0
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- c(
            dnorm(txgrid, m[i,acluster], s[i,acluster])*jac(xgrid)
        )
        q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i<nsubsamples|i==nsamples){
    if(i==nsamples){y <- ysum/nsamples}
    #y2 <- y2*max(y[-extr])/max(y2[-extr])
    ## y[extr] <- y[extr] * max(y[-extr])
    tplot(x=xgrid, y=y, #y=list(y,y2),
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=c(if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}},4), ly=1,lwd=0.5)
    tplot(x=xgrid[extr], y=y[extr],
          type='p',cex=0.075,col=3,add=T)
    ## tplot(x=xgrid[extr], y=y[extr], type='p',cex=0.1,add=T)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    abline(v=invtran(c(-1,1)),lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    abline(v=range(test),lwd=0.5,col=alpha2hex(0.5,2),lty=2)
    }
}
dev.off()

##############################################################################




#### Doubly-bounded case
#### Adding tails to last point
dx <- 2^-11
nsamples <- 400
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(2)
shape1s <- c(2) # large scales
shape2s <- c(1) # small scales
scales <- c(8)
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scale <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scale))),nsamples)
##
graphics.off()
pdff('samples_doublybounded_notransf_tails')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(-1, 1, length.out=256)
extr <- c(1,length(xgrid))
ysum <- 0
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- dnorm(xgrid, m[i,acluster], s[i,acluster])
        dens[extr[1]] <- dens[extr[1]]*0 + pnorm(xgrid[extr[1]], m[i,acluster], s[i,acluster])
        dens[extr[2]] <- dens[extr[2]]*0 + pnorm(xgrid[extr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
        q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i==nsamples){y <- ysum/nsamples}
    y[extr] <- y[extr] * max(y[-extr])
    tplot(x=xgrid[-extr], y=y[-extr],
          ylim=c(0,NA),xlim=c(-1-dx,1+dx),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=if(i<nsamples){1}else{if(any(is.infinite(ysum))){2}else{3}}, ly=1,lwd=0.5)
    tplot(x=xgrid[extr], y=y[extr],
          type='p',cex=0.1,add=T)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    abline(v=c(-1,1),lwd=0.5,col=alpha2hex(0.5,7),lty=2)
}
dev.off()


#### Doubly-bounded case
#### with transformation
dd <- 2^-11
tran <- function(x){qnorm(x*(1-2*dd)+dd)}
jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
dx <- 2^-11
##
nsamples <- 400
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(3)
shape1s <- c(1) # large scales
shape2s <- c(1/2) # small scales
scales <- c(8)
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scale <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scale))),nsamples)
##
graphics.off()
pdff('samples_doublybounded_transf')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(0, 1, length.out=256)
txgrid <- tran(xgrid)
extr <- c(1,length(xgrid))
ysum <- 0
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- c(
            pnorm(txgrid[extr[1]], m[i,acluster], s[i,acluster]),
            dnorm(txgrid[-extr], m[i,acluster], s[i,acluster])*jac(xgrid[-extr]),
            pnorm(txgrid[extr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
        )
        q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i==nsamples){y <- ysum/nsamples}
    y[extr] <- y[extr] * max(y[-extr])
    tplot(x=xgrid, y=y,
          ylim=c(0,NA),xlim=c(0-dx,1+dx),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=if(i<nsamples){1}else{if(any(is.infinite(ysum))){2}else{3}}, ly=1,lwd=0.5)
    tplot(x=xgrid[extr], y=y[extr],
          type='p',cex=0.1,add=T)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    abline(v=c(0,1),lwd=0.5,col=alpha2hex(0.5,7),lty=2)
}
dev.off()


set.seed(123)
#### Doubly-bounded case
#### with norm transformation @@@
pdist <- pnorm
qdist <- qnorm
ddist <- dnorm
dd <- pdist(qnorm(2^-6)) # this is the amount of probability left in a tail
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
xmin <- 0
xmax <- 1
xscale <- (xmax-xmin)/(1-2*dd)
xlocation <- xmin - dd*xscale
tran <- function(x){qdist((x-xlocation)/xscale)}
invtran <- function(y){pdist(y)*xscale+xlocation}
jac <- function(y){1/ddist(y)/xscale}
dx <- 1e-3
##
fract <- 400
nsamples <- fract*4
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(1)
shape1s <- c(1) # large scales
shape2s <- c(1) # small scales
scales <- c(1/4)^-2
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scaleprec <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scaleprec))),nsamples)
##
xgrid <- c(xmin,invtran(qdist(dd/8)),seq(xmin, xmax, length.out=256),invtran(qdist(1-dd/8)),xmax)
txgrid <- tran(xgrid)
extr <- c(1,length(xgrid))
extr2 <- c(2,length(xgrid)-1)
ysum <- 0
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
graphics.off()
pdff('samples_doublybounded_norm_extr')
par(mfrow=c(20,20),mar = c(0,0,0,0))
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- c(
            pnorm(txgrid[extr[1]], m[i,acluster], s[i,acluster]),
            dnorm(txgrid[-extr], m[i,acluster], s[i,acluster])*jac(txgrid[-extr]),
            pnorm(txgrid[extr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
        )
        q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i<fract | i==nsamples){
    if(i==nsamples){y <- ysum/nsamples}
    y[extr] <- y[extr] * max(y[-c(extr,extr2)],1/(xmax-xmin))
    tplot(x=xgrid[-c(extr,extr2)], y=y[-c(extr,extr2)],
          ylim=c(0,max(y[-c(extr,extr2)],1/(xmax-xmin))),xlim=range(xgrid[-extr2]),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=if(i<fract){1}else{if(any(is.infinite(ysum))){2}else{3}}, ly=1,lwd=0.5)
    ## tplot(x=xgrid[extr2], y=y[extr2],
    ##       type='p',col=4,cex=0.15,add=T,pch=3)
    tplot(x=xgrid[extr+c(2,-2)], y=y[extr+c(2,-2)],
          type='p',col=4,cex=0.075,add=T,pch=1)
    tplot(x=xgrid[extr], y=y[extr],
          type='p',cex=0.1,col=1,add=T,pch=3)
    abline(h=c(0),lwd=0.5,col=alpha2hex(0.5,c(7,2)),lty=c(1,2))
    if(i==nsamples){
        abline(h=c(1/(xmax-xmin)),lwd=0.5,col=alpha2hex(0.5,c(2)),lty=1)
    }
    abline(v=c(xmin,xmax),lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
}
dev.off()


set.seed(123)
#### Doubly-bounded case
#### with identity transformation @@@
pdist <- function(y){y/2+0.5}
qdist <- function(x){x*2-1}
ddist <- function(y){0.5}
dd <- pdist(qnorm(2^-6)) # this is the amount of probability left in a tail
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
xmin <- 0
xmax <- 150
xscale <- (xmax-xmin)/(1-2*dd)
xlocation <- xmin - dd*xscale
tran <- function(x){qdist((x-xlocation)/xscale)}
invtran <- function(y){pdist(y)*xscale+xlocation}
jac <- function(y){1/ddist(y)/xscale}
dx <- 1e-3
##
fract <- 400
nsamples <- fract*4
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(3)
shape1s <- c(1/2) # large scales
shape2s <- c(1/2) # small scales
scales <- c(1/2)^-2
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scaleprec <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scaleprec))),nsamples)
##
xgrid <- c(xmin,invtran(qdist(dd/8)),seq(xmin, xmax, length.out=256),invtran(qdist(1-dd/8)),xmax)
txgrid <- tran(xgrid)
extr <- c(1,length(xgrid))
extr2 <- c(2,length(xgrid)-1)
ysum <- 0
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
graphics.off()
pdff('samples_doublybounded_id_extr')
par(mfrow=c(20,20),mar = c(0,0,0,0))
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- c(
            pnorm(txgrid[extr[1]], m[i,acluster], s[i,acluster]),
            dnorm(txgrid[-extr], m[i,acluster], s[i,acluster])*jac(txgrid[-extr]),
            pnorm(txgrid[extr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
        )
        q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i<fract | i==nsamples){
    if(i==nsamples){y <- ysum/nsamples}
    y[extr] <- y[extr] * max(y[-c(extr,extr2)],1/(xmax-xmin))
    tplot(x=xgrid[-c(extr,extr2)], y=y[-c(extr,extr2)],
          ylim=c(0,max(y[-c(extr,extr2)],1/(xmax-xmin))),xlim=range(xgrid[-extr2]),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=if(i<fract){1}else{if(any(is.infinite(ysum))){2}else{3}}, ly=1,lwd=0.5)
    ## tplot(x=xgrid[extr2], y=y[extr2],
    ##       type='p',col=4,cex=0.15,add=T,pch=3)
    tplot(x=xgrid[extr+c(2,-2)], y=y[extr+c(2,-2)],
          type='p',col=4,cex=0.075,add=T,pch=1)
    tplot(x=xgrid[extr], y=y[extr],
          type='p',cex=0.1,col=1,add=T,pch=3)
    abline(h=c(0),lwd=0.5,col=alpha2hex(0.5,c(7,2)),lty=c(1,2))
    if(i==nsamples){
        abline(h=c(1/(xmax-xmin)),lwd=0.5,col=alpha2hex(0.5,c(2)),lty=1)
    }
    abline(v=c(xmin,xmax),lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
}
dev.off()

set.seed(123)
#### Doubly-bounded case
#### with logis transformation @@@
pdist <- plogis
qdist <- qlogis
ddist <- dlogis
dd <- pdist(qnorm(2^-6)) # this is the amount of probability left in a tail
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
xmin <- 0
xmax <- 1
xscale <- (xmax-xmin)/(1-2*dd)
xlocation <- xmin - dd*xscale
tran <- function(x){qdist((x-xlocation)/xscale)}
invtran <- function(y){pdist(y)*xscale+xlocation}
jac <- function(y){1/ddist(y)/xscale}
dx <- 1e-3
##
fract <- 400
nsamples <- fract*4
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(1.5)
shape1s <- c(1/2) # large scales
shape2s <- c(1) # small scales
scales <- c(1/4)^-2
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scaleprec <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scaleprec))),nsamples)
##
xgrid <- c(xmin,invtran(qdist(dd/8)),seq(xmin, xmax, length.out=256),invtran(qdist(1-dd/8)),xmax)
txgrid <- tran(xgrid)
extr <- c(1,length(xgrid))
extr2 <- c(2,length(xgrid)-1)
ysum <- 0
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
graphics.off()
pdff('samples_doublybounded_logis_extr')
par(mfrow=c(20,20),mar = c(0,0,0,0))
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- c(
            pnorm(txgrid[extr[1]], m[i,acluster], s[i,acluster]),
            dnorm(txgrid[-extr], m[i,acluster], s[i,acluster])*jac(txgrid[-extr]),
            pnorm(txgrid[extr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
        )
        q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i<fract | i==nsamples){
    if(i==nsamples){y <- ysum/nsamples}
    y[extr] <- y[extr] * max(y[-c(extr,extr2)],1/(xmax-xmin))
    tplot(x=xgrid[-c(extr,extr2)], y=y[-c(extr,extr2)],
          ylim=c(0,max(y[-c(extr,extr2)],1/(xmax-xmin))),xlim=range(xgrid[-extr2]),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=if(i<fract){1}else{if(any(is.infinite(ysum))){2}else{3}}, ly=1,lwd=0.5)
    ## tplot(x=xgrid[extr2], y=y[extr2],
    ##       type='p',col=4,cex=0.15,add=T,pch=3)
    tplot(x=xgrid[extr+c(2,-2)], y=y[extr+c(2,-2)],
          type='p',col=4,cex=0.075,add=T,pch=1)
    tplot(x=xgrid[extr], y=y[extr],
          type='p',cex=0.1,col=1,add=T,pch=3)
    abline(h=c(0),lwd=0.5,col=alpha2hex(0.5,c(7,2)),lty=c(1,2))
    if(i==nsamples){
        abline(h=c(1/(xmax-xmin)),lwd=0.5,col=alpha2hex(0.5,c(2)),lty=1)
    }
    abline(v=c(xmin,xmax),lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
}
dev.off()



nclusters <- 1
nsamples <- 1e5
shape1s <- c(1/2) # large scales
shape2s <- c(1/2) # small scales
scales <- 1^(-2) #
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scale <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- log10(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scale)))/2
thist(s,plot=T)


#### Doubly-bounded case
#### with logis transformation
dd <- 2^-10
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
tran <- function(x){qlogis(x*(1-2*dd)+dd)}
jac <- function(x){1/dlogis(tran(x))*(1-2*dd)}
dx <- 1e-3
##
fract <- 400
nsamples <- fract*4
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(1.5)
shape1s <- c(1.25) # large scales
shape2s <- c(1) # small scales
scales <-  1/2^(-2)
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scale <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scale))),nsamples)
##
xgrid <- c(0,-dd*7/8,seq(0, 1, length.out=256),1+dd*7/8,1)
txgrid <- tran(xgrid)
extr <- c(1,length(xgrid))
ysum <- 0
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
graphics.off()
pdff('samples_doublybounded_logis')
par(mfrow=c(20,20),mar = c(0,0,0,0))
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- c(
            pnorm(txgrid[extr[1]], m[i,acluster], s[i,acluster]),
            dnorm(txgrid[-extr], m[i,acluster], s[i,acluster])*jac(xgrid[-extr]),
            pnorm(txgrid[extr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
        )
        q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i<fract | i==nsamples){
    if(i==nsamples){y <- ysum/nsamples}
    y[extr] <- y[extr] * max(y[-extr])
    tplot(x=xgrid[-extr], y=y[-extr],
          ylim=c(0,max(y,1)),xlim=c(0-dd*7/8,1+dd*7/8),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=if(i<fract){1}else{if(any(is.infinite(ysum))){2}else{3}}, ly=1,lwd=0.5)
    tplot(x=xgrid[extr+c(1,-1)], y=y[extr+c(1,-1)],
          type='p',col=4,cex=0.15,add=T,pch=3)
    tplot(x=xgrid[extr+c(2,-2)], y=y[extr+c(2,-2)],
          type='p',cex=0.15,add=T,pch=2)
    tplot(x=xgrid[extr], y=y[extr],
          type='p',cex=0.075,col=2,add=T)
    abline(h=c(0),lwd=0.5,col=alpha2hex(0.5,c(7,2)),lty=c(1,2))
    if(i==nsamples){
        abline(h=c(1),lwd=0.5,col=alpha2hex(0.5,c(2)),lty=1)
    }
    abline(v=c(0,1),lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
}
dev.off()





#### Integer
#### without transformation
dx <- 2^-11
nint <- 8+1
nsamples <- 400
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(1)
shape1s <- c(2) # large scales
shape2s <- c(1/2) # small scales
scales <- c(8)
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scale <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scale))),nsamples)
##
graphics.off()
pdff('samples_integer')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(-1,1,length.out=nint)
dx <- 2/(nint-1)
extr <- c(1,length(xgrid))
ysum <- 0
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- pnorm(xgrid+dx, m[i,acluster], s[i,acluster]) - pnorm(xgrid-dx, m[i,acluster], s[i,acluster])
        dens[extr[1]] <- pnorm(xgrid[extr[1]]+dx, m[i,acluster], s[i,acluster])
        dens[extr[2]] <- pnorm(xgrid[extr[2]]-dx, m[i,acluster], s[i,acluster], lower.tail=F)
        q[i,acluster]*dens}))
    ## y2 <- rowSums(sapply(1:nclusters,function(acluster){
    ##     dens <- dnorm(xgrid, m[i,acluster], s[i,acluster])*2*dx
    ##     dens[extr[1]] <- pnorm(xgrid[extr[1]], m[i,acluster], s[i,acluster])
    ##     dens[extr[2]] <- pnorm(xgrid[extr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
    ##     q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i==nsamples){y <- ysum/nsamples}
    #y2 <- y2*max(y[-extr])/max(y2[-extr])
    ## y[extr] <- y[extr] * max(y[-extr])
    tplot(x=xgrid, y=y, #y=list(y,y2),
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=c(if(i<nsamples){1}else{if(any(is.infinite(ysum))){2}else{3}},4), ly=1,lwd=0.5)
    tplot(x=xgrid, y=y,type='p',cex=0.2,
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,add=T,
          col=if(i<nsamples){1}else{if(any(is.infinite(ysum))){2}else{3}}, ly=1,lwd=0.5)
    ## tplot(x=xgrid[extr], y=y[extr], type='p',cex=0.1,add=T)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    ## abline(v=,lwd=0.5,col=alpha2hex(0.5,7),lty=2)
}
dev.off()




#### Integer
#### no transformation
nint <- 16
dd <- qnorm(0.5/nint)/2 + 0.5
tran <- function(x){(x*(1-2*dd)+dd)*2-1}
## dd <- 2^-11
##
nsamples <- 400*8
nsubsamples <- 400
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(1)
shape1s <- c(2) # large scales
shape2s <- c(1/2) # small scales
scales <- 1/4^(-2)
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scale <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scale))),nsamples)
##
graphics.off()
pdff('samples_integer_notrans')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(0,1,length.out=nint)
extr <- c(1,length(xgrid))
mgrid <- (xgrid[-extr[2]]+xgrid[-extr[1]])/2
mextr <- c(1,length(mgrid))
txgrid <- tran(xgrid)
tmgrid <- tran(mgrid)
##tmgrid <- (txgrid[-extr[2]]+txgrid[-extr[1]])/2
dx <- 2/(nint-1)
ysum <- 0
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- c(
            pnorm(tmgrid[mextr[1]], m[i,acluster], s[i,acluster]),
            pnorm(tmgrid[-mextr[1]], m[i,acluster], s[i,acluster]) - pnorm(tmgrid[-mextr[2]], m[i,acluster], s[i,acluster]),
            pnorm(tmgrid[mextr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
        )
        q[i,acluster]*dens}))
    ## y2 <- rowSums(sapply(1:nclusters,function(acluster){
    ##     dens <- dnorm(xgrid, m[i,acluster], s[i,acluster])*2*dx
    ##     dens[extr[1]] <- pnorm(xgrid[extr[1]], m[i,acluster], s[i,acluster])
    ##     dens[extr[2]] <- pnorm(xgrid[extr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
    ##     q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i<nsubsamples|i==nsamples){
    if(i==nsamples){y <- ysum/nsamples}
    #y2 <- y2*max(y[-extr])/max(y2[-extr])
    ## y[extr] <- y[extr] * max(y[-extr])
    tplot(x=xgrid, y=y, #y=list(y,y2),
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=c(if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}},4), ly=1,lwd=0.5)
    tplot(x=xgrid, y=y,type='p',cex=0.2,
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,add=T,
          col=if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}}, ly=1,lwd=0.5)
    tplot(x=xgrid[extr], y=y[extr],
          type='p',cex=0.075,col=3,add=T)
    ## tplot(x=xgrid[extr], y=y[extr], type='p',cex=0.1,add=T)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    if(i==nsamples){
        abline(h=c(1/nint),lwd=0.5,col=alpha2hex(0.5,c(2)),lty=1)
    }
    ## abline(v=,lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
}
dev.off()



#### Integer
#### with norm transformation
nint <- 16
dd <- pnorm(qnorm(0.5/nint))
tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## dd <- 2^-11
##
nsamples <- 400*8
nsubsamples <- 400
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(1)
shape1s <- c(2) # large scales
shape2s <- c(1) # small scales
scales <- 1/4^(-2)
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scale <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scale))),nsamples)
##
graphics.off()
pdff('samples_integer_norm')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(0,1,length.out=nint)
extr <- c(1,length(xgrid))
mgrid <- (xgrid[-extr[2]]+xgrid[-extr[1]])/2
mextr <- c(1,length(mgrid))
txgrid <- tran(xgrid)
tmgrid <- tran(mgrid)
##tmgrid <- (txgrid[-extr[2]]+txgrid[-extr[1]])/2
dx <- 2/(nint-1)
ysum <- 0
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- c(
            pnorm(tmgrid[mextr[1]], m[i,acluster], s[i,acluster]),
            pnorm(tmgrid[-mextr[1]], m[i,acluster], s[i,acluster]) - pnorm(tmgrid[-mextr[2]], m[i,acluster], s[i,acluster]),
            pnorm(tmgrid[mextr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
        )
        q[i,acluster]*dens}))
    ## y2 <- rowSums(sapply(1:nclusters,function(acluster){
    ##     dens <- dnorm(xgrid, m[i,acluster], s[i,acluster])*2*dx
    ##     dens[extr[1]] <- pnorm(xgrid[extr[1]], m[i,acluster], s[i,acluster])
    ##     dens[extr[2]] <- pnorm(xgrid[extr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
    ##     q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i<nsubsamples|i==nsamples){
    if(i==nsamples){y <- ysum/nsamples}
    #y2 <- y2*max(y[-extr])/max(y2[-extr])
    ## y[extr] <- y[extr] * max(y[-extr])
    tplot(x=xgrid, y=y, #y=list(y,y2),
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=c(if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}},4), ly=1,lwd=0.5)
    tplot(x=xgrid, y=y,type='p',cex=0.2,
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,add=T,
          col=if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}}, ly=1,lwd=0.5)
    tplot(x=xgrid[extr], y=y[extr],
          type='p',cex=0.075,col=3,add=T)
    ## tplot(x=xgrid[extr], y=y[extr], type='p',cex=0.1,add=T)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    if(i==nsamples){
        abline(h=c(1/nint),lwd=0.5,col=alpha2hex(0.5,c(2)),lty=1)
    }
    ## abline(v=,lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
}
dev.off()

set.seed(111)
#### Integer
#### with norm transformation II
nint <- 75
dd <- pnorm(qnorm(0.5/nint))
const <- qnorm(1-dd)
tran <- function(x){qnorm(x*(1-2*dd)+dd)/const}
## dd <- 2^-11
##
nsamples <- 400*8
nsubsamples <- 400
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(1)/const
shape1s <- c(2) # large scales
shape2s <- c(1) # small scales
scales <- (1/4/const)^(-2)
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scale <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scale))),nsamples)
##
graphics.off()
pdff('samples_integer_normII')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(0,1,length.out=nint)
extr <- c(1,length(xgrid))
mgrid <- (xgrid[-extr[2]]+xgrid[-extr[1]])/2
mextr <- c(1,length(mgrid))
txgrid <- tran(xgrid)
tmgrid <- tran(mgrid)
##tmgrid <- (txgrid[-extr[2]]+txgrid[-extr[1]])/2
dx <- 2/(nint-1)
ysum <- 0
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- c(
            pnorm(tmgrid[mextr[1]], m[i,acluster], s[i,acluster]),
            pnorm(tmgrid[-mextr[1]], m[i,acluster], s[i,acluster]) - pnorm(tmgrid[-mextr[2]], m[i,acluster], s[i,acluster]),
            pnorm(tmgrid[mextr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
        )
        q[i,acluster]*dens}))
    ## y2 <- rowSums(sapply(1:nclusters,function(acluster){
    ##     dens <- dnorm(xgrid, m[i,acluster], s[i,acluster])*2*dx
    ##     dens[extr[1]] <- pnorm(xgrid[extr[1]], m[i,acluster], s[i,acluster])
    ##     dens[extr[2]] <- pnorm(xgrid[extr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
    ##     q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i<nsubsamples|i==nsamples){
    if(i==nsamples){y <- ysum/nsamples}
    #y2 <- y2*max(y[-extr])/max(y2[-extr])
    ## y[extr] <- y[extr] * max(y[-extr])
    tplot(x=xgrid, y=y, #y=list(y,y2),
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=c(if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}},4), ly=1,lwd=0.5)
    tplot(x=xgrid, y=y,type='p',cex=0.2,
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,add=T,
          col=if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}}, ly=1,lwd=0.5)
    tplot(x=xgrid[extr], y=y[extr],
          type='p',cex=0.075,col=3,add=T)
    ## tplot(x=xgrid[extr], y=y[extr], type='p',cex=0.1,add=T)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    if(i==nsamples){
        abline(h=c(1/nint),lwd=0.5,col=alpha2hex(0.5,c(2)),lty=1)
    }
    ## abline(v=,lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
}
dev.off()



#### Integer
#### with norm transformation III
nint <- 16
dd <- pnorm(-1)
const <- -qnorm(0.5/nint)
tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## dd <- 2^-11
##
nsamples <- 400*8
nsubsamples <- 400
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(1)/const
shape1s <- c(1) # large scales
shape2s <- c(1) # small scales
scales <- (1/2/const)^(-2)
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scale <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scale))),nsamples)
##
graphics.off()
pdff('samples_integer_normIII')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(0,1,length.out=nint)
extr <- c(1,length(xgrid))
mgrid <- (xgrid[-extr[2]]+xgrid[-extr[1]])/2
mextr <- c(1,length(mgrid))
txgrid <- tran(xgrid)
tmgrid <- tran(mgrid)
##tmgrid <- (txgrid[-extr[2]]+txgrid[-extr[1]])/2
dx <- 2/(nint-1)
ysum <- 0
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- c(
            pnorm(tmgrid[mextr[1]], m[i,acluster], s[i,acluster]),
            pnorm(tmgrid[-mextr[1]], m[i,acluster], s[i,acluster]) - pnorm(tmgrid[-mextr[2]], m[i,acluster], s[i,acluster]),
            pnorm(tmgrid[mextr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
        )
        q[i,acluster]*dens}))
    ## y2 <- rowSums(sapply(1:nclusters,function(acluster){
    ##     dens <- dnorm(xgrid, m[i,acluster], s[i,acluster])*2*dx
    ##     dens[extr[1]] <- pnorm(xgrid[extr[1]], m[i,acluster], s[i,acluster])
    ##     dens[extr[2]] <- pnorm(xgrid[extr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
    ##     q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i<nsubsamples|i==nsamples){
    if(i==nsamples){y <- ysum/nsamples}
    #y2 <- y2*max(y[-extr])/max(y2[-extr])
    ## y[extr] <- y[extr] * max(y[-extr])
    tplot(x=xgrid, y=y, #y=list(y,y2),
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=c(if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}},4), ly=1,lwd=0.5)
    tplot(x=xgrid, y=y,type='p',cex=0.2,
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,add=T,
          col=if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}}, ly=1,lwd=0.5)
    tplot(x=xgrid[extr], y=y[extr],
          type='p',cex=0.075,col=3,add=T)
    ## tplot(x=xgrid[extr], y=y[extr], type='p',cex=0.1,add=T)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    if(i==nsamples){
        abline(h=c(1/nint),lwd=0.5,col=alpha2hex(0.5,c(2)),lty=1)
    }
    ## abline(v=,lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
}
dev.off()





#### Integer
#### with logis transformation
nint <- 16
dd <- plogis(qnorm(0.5/nint))
tran <- function(x){qlogis(x*(1-2*dd)+dd)}
##
nint <- 75
nsamples <- 400*8
nsubsamples <- 400
nclusters <- 64
alphas <- c(1,2,0.5)
means <- c(0)
sds <- c(1)
shape1s <- c(1) # large scales
shape2s <- c(1/2) # small scales
scales <- 1/4^(-2)
##
alpha <- sample(rep(alphas,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sds,2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,means,sd),nsamples)
shape1 <- sample(rep(shape1s,2),nsamples*nclusters,replace=T)
shape2 <- sample(rep(shape2s,2),nsamples*nclusters,replace=T)
scale <- sample(rep(scales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shape1,rate=nimble::rinvgamma(nsamples*nclusters,shape=shape2,scale=scale))),nsamples)
##
graphics.off()
pdff('samples_integer_logis')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(0,1,length.out=nint)
extr <- c(1,length(xgrid))
mgrid <- (xgrid[-extr[2]]+xgrid[-extr[1]])/2
mextr <- c(1,length(mgrid))
txgrid <- tran(xgrid)
tmgrid <- tran(mgrid)
dx <- 2/(nint-1)
ysum <- 0
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        dens <- c(
            pnorm(tmgrid[mextr[1]], m[i,acluster], s[i,acluster]),
            pnorm(tmgrid[-mextr[1]], m[i,acluster], s[i,acluster]) - pnorm(tmgrid[-mextr[2]], m[i,acluster], s[i,acluster]),
            pnorm(tmgrid[mextr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
        )
        q[i,acluster]*dens}))
    ## y2 <- rowSums(sapply(1:nclusters,function(acluster){
    ##     dens <- dnorm(xgrid, m[i,acluster], s[i,acluster])*2*dx
    ##     dens[extr[1]] <- pnorm(xgrid[extr[1]], m[i,acluster], s[i,acluster])
    ##     dens[extr[2]] <- pnorm(xgrid[extr[2]], m[i,acluster], s[i,acluster], lower.tail=F)
    ##     q[i,acluster]*dens}))
    ysum <- ysum+y
    if(i<nsubsamples|i==nsamples){
    if(i==nsamples){y <- ysum/nsamples}
    #y2 <- y2*max(y[-extr])/max(y2[-extr])
    ## y[extr] <- y[extr] * max(y[-extr])
    tplot(x=xgrid, y=y, #y=list(y,y2),
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=c(if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}},4), ly=1,lwd=0.5)
    tplot(x=xgrid, y=y,type='p',cex=0.2,
          ylim=c(0,NA),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,add=T,
          col=if(i<nsubsamples){1}else{if(any(is.infinite(ysum))){2}else{3}}, ly=1,lwd=0.5)
    tplot(x=xgrid[extr], y=y[extr],
          type='p',cex=0.075,col=3,add=T)
    abline(h=c(0),lwd=0.5,col=alpha2hex(0.5,c(7,2)),lty=c(1,2))
    ## tplot(x=xgrid[extr], y=y[extr], type='p',cex=0.1,add=T)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    if(i==nsamples){
        abline(h=c(1/nint),lwd=0.5,col=alpha2hex(0.5,c(2)),lty=1)
    }
    ## abline(v=,lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
}
dev.off()




##
## fb <- function(x,d){log(x*(1-2*d)+d)-log(1-d-x*(1-2*d))}
dd <- 2^-11
tran <- function(x){qlogis(x*(1-2*dd)+dd)}

dgrid <- seq(-12,-10.5,length.out=512)
ddiff <- 2*pnorm(fb(0,2^dgrid))-pnorm(fb(1/2^14,2^dgrid))
tplot(dgrid,ddiff)
## 2^-11

q2 <- 0.83
q1 <- 0.5
q3 <- 0.92

dd <- 2^-11
nsam <- 400
k <- 64
alch <- c(1)
sdch <- c(1.5)
sh1ch <- c(1.5)
sh2ch <- c(1.5)
scch <- c(8)
##
muu <- fb0(q2,dd)
sii <- fb0(q3,dd)-fb0(q1,dd)
fb <- function(x){(fb0(x,dd)-muu)/sii}
##
al <- sample(rep(alch,2),nsam,replace=T)
q <- extraDistr::rdirichlet(n=nsam,alpha=matrix(al/k,nsam,k))
sd <- sample(rep(sdch,2),nsam*k,replace=T)
m <- matrix(rnorm(nsam*k,0,sd),nsam)
sh1 <- sample(rep(sh1ch,2),nsam*k,replace=T)
sh2 <- sample(rep(sh2ch,2),nsam*k,replace=T)
sc <- sample(rep(scch,2),nsam*k,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsam*k,shape=sh1,rate=nimble::rinvgamma(nsam*k,shape=sh2,scale=sc))),nsam)
##
graphics.off()
pdff('samples_doublybounded')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(0, 1, length.out=128)
ysum <- 0
for(i in 1:nsam){
    y <- rowSums(sapply(1:k,function(cc)q[i,cc]*dnorm(fb(xgrid), m[i,cc], s[i,cc])))/sii/dlogis(xgrid*(1-2*dd)+dd)*(1-2*dd)
    ysum <- ysum+y
    tplot(x=xgrid,y=if(i<nsam){y}else{ysum/nsam},ylim=c(0,NA),ylabels=NA,xlabels=NA,
          mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,col=if(i<nsam){1}else{2},
          ly=1,lwd=0.5,xlim=c(-2*dd,1+2*dd),
          xticks=NA,yticks=NA)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    abline(v=c(dd,1-dd),lwd=0.5,col=alpha2hex(0.5,7),lty=1)
}
dev.off()





## fb <- function(x,d){log(x*(1-2*d)+d)-log(1-d-x*(1-2*d))}
fb <- function(x,d){qlogis(x*(1-2*d)+d)}
dfb <- function(x,d){}

dgrid <- seq(-12,-10.5,length.out=512)
ddiff <- 2*pnorm(fb(0,2^dgrid))-pnorm(fb(1/2^14,2^dgrid))
tplot(dgrid,ddiff)
## 2^-11

dd <- 2^-11
nsam <- 400
k <- 64
alch <- c(1/2,1,2)
sdch <- c(1.5)
sh1ch <- c(1.5)
sh2ch <- c(1.5)
scch <- c(8)
##
al <- sample(rep(alch,2),nsam,replace=T)
q <- extraDistr::rdirichlet(n=nsam,alpha=matrix(al/k,nsam,k))
sd <- sample(rep(sdch,2),nsam*k,replace=T)
m <- matrix(rnorm(nsam*k,0,sd),nsam)
sh1 <- sample(rep(sh1ch,2),nsam*k,replace=T)
sh2 <- sample(rep(sh2ch,2),nsam*k,replace=T)
sc <- sample(rep(scch,2),nsam*k,replace=T)
    s <- matrix(sqrt(nimble::rinvgamma(nsam*k,shape=sh1,rate=nimble::rinvgamma(nsam*k,shape=sh2,scale=sc))),nsam)
##
graphics.off()
pdff('test3.pdf')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(0, 1, length.out=128)
ysum <- 0
for(i in 1:nsam){
    y <- rowSums(sapply(1:k,function(cc)q[i,cc]*dnorm(fb(xgrid,dd), m[i,cc], s[i,cc])))/dlogis(fb(xgrid,dd))
    ysum <- ysum+y
    tplot(x=xgrid,y=if(i<nsam){y}else{ysum/nsam},ylim=c(0,NA),ylabels=NA,xlabels=NA,
          mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,col=if(i<nsam){1}else{2},
          ly=1,lwd=0.5,xlim=c(-2*dd,1+2*dd),
          xticks=NA,yticks=NA)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    abline(v=c(0,1),lwd=0.5,col=alpha2hex(0.5,7),lty=1)
}
dev.off()



## fb <- function(x,d){log(x*(1-2*d)+d)-log(1-d-x*(1-2*d))}
fb <- function(x,d){log(x*(1-d)+d)}

dgrid <- seq(-12,-10.5,length.out=512)
ddiff <- 2*pnorm(fb(0,2^dgrid))-pnorm(fb(1/2^14,2^dgrid))
tplot(dgrid,ddiff)
## 2^-11

dd <- 2^-10.75
nsam <- 400
k <- 64
alch <- c(1/2,1,2)
sdch <- c(2)
sh1ch <- c(2)
sh2ch <- c(1)
scch <- c(4)
##
al <- sample(rep(alch,2),nsam,replace=T)
q <- extraDistr::rdirichlet(n=nsam,alpha=matrix(al/k,nsam,k))
sd <- sample(rep(sdch,2),nsam*k,replace=T)
m <- matrix(rnorm(nsam*k,0,sd),nsam)
sh1 <- sample(rep(sh1ch,2),nsam*k,replace=T)
sh2 <- sample(rep(sh2ch,2),nsam*k,replace=T)
sc <- sample(rep(scch,2),nsam*k,replace=T)
    s <- matrix(sqrt(nimble::rinvgamma(nsam*k,shape=sh1,rate=nimble::rinvgamma(nsam*k,shape=sh2,scale=sc))),nsam)
##
graphics.off()
pdff('test3log')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(0, 3, length.out=128)
ysum <- 0
for(i in 1:nsam){
    y <- rowSums(sapply(1:k,function(cc)q[i,cc]*dnorm(fb(xgrid,dd), m[i,cc], s[i,cc])))*dd/(xgrid*(1-dd)+dd)
    ysum <- ysum+y
    tplot(x=xgrid,y=if(i<nsam){y}else{ysum/nsam},ylim=c(0,NA),ylabels=NA,xlabels=NA,
          mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,col=if(i<nsam){1}else{2},
          ly=1,lwd=0.5,xlim=c(-dd,3),
          xticks=NA,yticks=NA)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    abline(v=c(0,1),lwd=0.5,col=alpha2hex(0.5,7),lty=1)
}
dev.off()





## fb <- function(x,d){log(x*(1-2*d)+d)-log(1-d-x*(1-2*d))}
fb <- function(x,d){log(x*(1-d)+d)}

dgrid <- seq(-12,-10.5,length.out=512)
ddiff <- 2*pnorm(fb(0,2^dgrid))-pnorm(fb(1/2^14,2^dgrid))
tplot(dgrid,ddiff)
## 2^-11

dd <- 2^-10.75
nsam <- 400
k <- 64
alch <- c(1/2,1,2)
sdch <- c(2)
sh1ch <- c(0.5)
sh2ch <- c(0.5)
scch <- c(1)
##
al <- sample(rep(alch,2),nsam,replace=T)
q <- extraDistr::rdirichlet(n=nsam,alpha=matrix(al/k,nsam,k))
sd <- sample(rep(sdch,2),nsam*k,replace=T)
m <- matrix(rnorm(nsam*k,0,sd),nsam)
sh1 <- sample(rep(sh1ch,2),nsam*k,replace=T)
sh2 <- sample(rep(sh2ch,2),nsam*k,replace=T)
sc <- sample(rep(scch,2),nsam*k,replace=T)
    s <- matrix(sqrt(nimble::rinvgamma(nsam*k,shape=sh1,rate=nimble::rinvgamma(nsam*k,shape=sh2,scale=sc))),nsam)
##
graphics.off()
pdff('test3u')
par(mfrow=c(20,20),mar = c(0,0,0,0))
xgrid <- seq(-3.5, 3.5, length.out=128)
ysum <- 0
for(i in 1:nsam){
    y <- rowSums(sapply(1:k,function(cc)q[i,cc]*dnorm(xgrid, m[i,cc], s[i,cc])))
    ysum <- ysum+y
    tplot(x=xgrid,y=if(i<nsam){y}else{ysum/nsam},ylim=c(0,NA),ylabels=NA,xlabels=NA,
          mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,col=if(i<nsam){1}else{2},
          ly=1,lwd=0.5,xlim=c(-3.5,3.5),
          xticks=NA,yticks=NA)
    abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
    abline(v=c(-1,1),lwd=0.5,col=alpha2hex(0.5,7),lty=1)
}
dev.off()













































## fb <- function(x){pnorm(x)*(nn+2)/nn-1/nn}
fb <- function(x){qnorm((nn*x+1)/(nn+2))}
fl <- function(x){log(x)}
##
nsam <- 400
##
## combinations <- t(expand.grid(al,sd,sh1,sh2,sc))
k <- 64
##

alch <- c(1,0.5,2)
sdch <- c(3)
sh1ch <- c(1,1/2)
sh2ch <- c(1,1/2)
scch <- c(1)
    ##
    al <- sample(rep(alch,2),nsam,replace=T)
q <- extraDistr::rdirichlet(n=nsam,alpha=matrix(al/k,nsam,k))
    sd <- sample(rep(sdch,2),nsam*k,replace=T)
m <- matrix(rnorm(nsam*k,0,sd),nsam)
    sh1 <- sample(rep(sh1ch,2),nsam*k,replace=T)
    sh2 <- sample(rep(sh2ch,2),nsam*k,replace=T)
    sc <- sample(rep(scch,2),nsam*k,replace=T)
##
    s <- matrix(sqrt(nimble::rinvgamma(nsam*k,shape=sh1,rate=nimble::rinvgamma(nsam*k,shape=sh2,scale=sc))),nsam)
##
graphics.off()
    title <- 'varioushpdiscrete'
    pdff(paste0('priors_',title))
    ##
    par(mfrow=c(20,20),mar = c(0,0,0,0))
    xgrid <- seq(-3, 3, length.out=128)
    for(i in 1:nsam){
        y <- rowSums(sapply(1:k,function(cc)q[i,cc]*dnorm(xgrid, m[i,cc], s[i,cc])))
        tplot(x=xgrid,y=y,ylim=c(0,NA),ylabels=NA,xlabels=NA,
              mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,
              ly=1,lwd=0.5,xlim=c(NA,NA),
              xticks=NA,yticks=NA)
            abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
        abline(v=c(-1,1),lwd=0.5,col=alpha2hex(0.5,7),lty=1)
        }
##
par(mfrow=c(20,20),mar = c(0,0,0,0))
    xgrid <- seq(0, 1, length.out=128)
    for(i in 1:nsam){
        y <- rowSums(sapply(1:k,function(cc)q[i,cc]*dnorm(fb(xgrid), m[i,cc], s[i,cc])))*dnorm((nn*xgrid+1)/(nn+2))*nn/(nn+2)
        tplot(x=xgrid,y=y,ylim=c(0,NA),ylabels=NA,xlabels=NA,
              mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,
              ly=1,lwd=0.5,xlim=c(-2/nn, 1+2/nn),
              xticks=NA,yticks=NA)
            abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
        abline(v=c(0,1),lwd=0.5,col=alpha2hex(0.5,7),lty=1)
        }
    ##
    par(mfrow=c(20,20),mar = c(0,0,0,0))
    xgrid <- seq(0, 3, length.out=128)
    for(i in 1:nsam){
        y <- rowSums(sapply(1:k,function(cc)q[i,cc]*dnorm(fl(xgrid), m[i,cc], s[i,cc])))/xgrid
        tplot(x=xgrid,y=y,ylim=c(0,NA),ylabels=NA,xlabels=NA,
              mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,
              ly=1,lwd=0.5,xlim=c(-0.1,NA),
              xticks=NA,yticks=NA)
            abline(h=0,lwd=0.5,col=alpha2hex(0.5,7),lty=1)
        abline(v=c(0,1),lwd=0.5,col=alpha2hex(0.5,7),lty=1)
        }
dev.off()


    NULL
    }






dt <- fread('~/repositories/ADBayes/_scripts2/data_transformed_shuffled.csv')

outp <- t(sapply(dt2,function(x){c(uni=length(unique(x)),occ=length(x)/length(unique(x)),D=min(diff(sort(unique(x)))),rD=min(diff(sort(unique(x))))/diff(range(x)),min=min(x),max=max(x),tquant(x,(1:3)/4))}))
signif(outp[order(outp[,1]),],2)


v

##                  uni   occ       D      rD     min      max     25%      50%      75%
## Apoe4_             2 340.0 1.0e+00 1.0e+00  0.0000   1.0000  0.0000   1.0000   1.0000
## Subgroup_num_      2 340.0 1.0e+00 1.0e+00  0.0000   1.0000  0.0000   0.0000   1.0000
## Gender_num_        2 340.0 1.0e+00 1.0e+00  0.0000   1.0000  0.0000   0.0000   1.0000
## GDTOTAL_gds        7  97.0 1.0e+00 1.7e-01  0.0000   6.0000  1.0000   1.0000   2.0000
## AVDEL30MIN_neuro  16  42.0 1.0e+00 6.7e-02  0.0000  15.0000  0.0000   3.0000   6.0000
## AVDELTOT_neuro    16  42.0 1.0e+00 6.7e-02  0.0000  15.0000  8.0000  11.0000  13.0000
## CATANIMSC_neuro   29  23.0 1.0e+00 3.3e-02  5.0000  35.0000 13.0000  17.0000  20.0000
## ANARTERR_neuro    47  14.0 1.0e+00 2.0e-02  0.0000  50.0000  6.0000  10.0000  18.0000
## RAVLT_immediate   52  13.0 1.0e+00 1.8e-02 11.0000  67.0000 26.0000  32.0000  39.0000
## TRAASCOR_neuro    84   8.1 1.0e+00 7.3e-03 13.0000 150.0000 30.0000  37.0000  48.0000
## TRABSCOR_neuro   190   3.6 1.0e+00 3.7e-03 33.0000 300.0000 73.0000 100.0000 140.0000
## AGE              270   2.5 1.0e-01 2.7e-03 55.0000  91.0000 69.0000  74.0000  79.0000
## LRHHC_n_long     680   1.0 1.6e-09 3.9e-07  0.0023   0.0065  0.0037   0.0042   0.0047

Apoe4_
Subgroup_num_    
Gender_num_      
GDTOTAL_gds      
AVDEL30MIN_neuro 
AVDELTOT_neuro   
CATANIMSC_neuro  
ANARTERR_neuro   
RAVLT_immediate  
TRAASCOR_neuro   
TRABSCOR_neuro   
AGE              
LRHHC_n_long
