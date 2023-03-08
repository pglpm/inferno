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

#### Currently used ####

#### continuous variate
sdovermad <- 1/qnorm(0.75)
sdovermad2 <- 0.5/qnorm(0.75)
## dt <- fread('~/repositories/ledley-jaynes_machine/scripts/ingrid_data_nogds6.csv')
## varinfo <- data.matrix(read.csv('~/repositories/ledley-jaynes_machine/scripts/varinfo.csv',row.names=1))
set.seed(123)
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
xlocation <- 0
xscale <- 1
xmin <- -5
xmax <- 5
tran <- function(x){(x-xlocation)/xscale}
invtran <- function(y){y*xscale+xlocation}
jac <- function(y){1/xscale}
##
## hyperparameters
rowcol <- c(20,20)
nsamples <- 1e6
nclusters <- 64
alpha0 <- 2^((-3):3)
rmean0 <- 0
zeta <- 1
rvar0 <- (zeta)^2
rshapein0 <- 1 # large scales
rshapeout0 <- 1 # small scales
hwidth <- 2 # number of powers of 2 to consider in either direction
rvarscales <- (zeta * 2^((-hwidth):hwidth))^2
##
xsamples <- rnorm(nsamples,
                  mean=rnorm(nsamples,mean=rmean0,sd=sqrt(rvar0)),
                  sd=sqrt(
                      extraDistr::rbetapr(nsamples,shape1=rshapein0,shape2=rshapeout0,
                                          scale=sample(rvarscales,nsamples,replace=T))
                  )
                  )
IQR(xsamples)/2
mad(xsamples, constant=1)
IQR(xsamples)*sdovermad2/2
rm(xsamples)
## thist(xsamples[xsamples<6&xsamples>-6],plot=T)
## abline(v=c(-1,1))
graphics.off()
pdff('priorsamples_real')
nsamples2 <- min(2^14,nsamples)
alpha <- sample(rep(alpha0,2),nsamples2,replace=T)
q <- extraDistr::rdirichlet(n=nsamples2,alpha=matrix(alpha/nclusters,nsamples2,nclusters))
sd <- sample(rep(sqrt(rvar0),2),nsamples2*nclusters,replace=T)
m <- matrix(rnorm(nsamples2*nclusters,rmean0,sd),nsamples2)
shapein <- sample(rep(rshapein0,2),nsamples2*nclusters,replace=T)
shapeout <- sample(rep(rshapeout0,2),nsamples2*nclusters,replace=T)
scalevar <- sample(rep(rvarscales,2),nsamples2*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples2*nclusters,shape=shapeout,rate=nimble::rinvgamma(nsamples2*nclusters,shape=shapein,rate=scalevar))),nsamples2)
##
txgrid <- seq(xmin, xmax, length.out=256)
xgrid <- invtran(txgrid)
ysum <- 0
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
par(mfrow=rowcol,mar = c(0,0,0,0))
for(i in 1:nsamples2){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        q[i,acluster] *
            dnorm(txgrid, m[i,acluster], s[i,acluster]) * jac(txgrid)
    }))
    ysum <- ysum+y
    if(i < prod(rowcol) | i==nsamples2){
    if(i == nsamples2){y <- ysum/nsamples2}
    ## if(!is.null(data)){
    ##     his <- thist(data)
    ##     ymax <- max(y,his$density)
    ## }else{ymax <- NULL}
    ymax <- NULL
    tplot(x=xgrid, y=y,
          ylim=c(0,max(y,ymax)),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=(if(i < prod(rowcol)){1}else{if(any(is.infinite(ysum))){2}else{3}}), ly=1,lwd=(if(i < prod(rowcol)){0.5}else{1}))
    ## tplot(x=xgrid[extr2], y=y[extr2],
    ##       type='p',col=4,cex=0.15,add=T,pch=3)
    ## if(!is.null(data)){
    ##     tplot(x=his$mids,y=his$density,type='l',lwd=0.5,add=T,alpha=0.25,col=4)
    ## }
    abline(h=c(0),lwd=0.5,col=alpha2hex2(0.5,c(7,2)),lty=c(1,2))
##    if(i==nsamples2){
    if(TRUE){
        abline(v=c(-1,1),lwd=0.5,col=alpha2hex(7,0.5),lty=1)
        ## abline(v=c(xmin,xmax),lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
    }
}
dev.off()

#### continuous variate - log transformation
sdovermad2 <- 0.5/qnorm(0.75)
## dt <- fread('~/repositories/ledley-jaynes_machine/scripts/ingrid_data_nogds6.csv')
## varinfo <- data.matrix(read.csv('~/repositories/ledley-jaynes_machine/scripts/varinfo.csv',row.names=1))
graphics.off()
pdff('priorsamples_real_log')
set.seed(123)
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
xlocation <- 0
xq1 <- 1
xq3 <- 3
xscale <- (log(xq3)-log(xq1))*sdovermad2/2
xlocation <- mean(log(c(xq1,xq3)))
xmin <- 2^(-10)
xmax <- 4
tran <- function(x){(log(x)-xlocation)/xscale}
invtran <- function(y){exp(xscale*y+xlocation)}
jac <- function(y){exp(-(xscale*y+xlocation))/xscale}
##
## hyperparameters
rowcol <- c(20,20)
nsamples <- 1e6
nclusters <- 64
alpha0 <- 2^((-3):3)
rmean0 <- 0
rvar0 <- (0+2*sdovermad2)^2
rshapein0 <- 1 # large scales
rshapeout0 <- 1 # small scales
hwidth <- 2 # number of powers of 2 to consider in either direction
rvarscales <- ((0+1*sdovermad2) * 2^((-hwidth):hwidth))^2
##
xsamples <- rnorm(nsamples,
                  mean=rnorm(nsamples,mean=rmean0,sd=sqrt(rvar0)),
                  sd=sqrt(
                      extraDistr::rbetapr(nsamples,shape1=rshapein0,shape2=rshapeout0,
                                          scale=sample(rvarscales,nsamples,replace=T))
                  )
                  )
IQR(invtran(xsamples))/2
mad(invtran(xsamples), constant=1)
IQR(invtran(xsamples))*sdovermad2/2
rm(xsamples)
## thist(xsamples[xsamples<6&xsamples>-6],plot=T)
## abline(v=c(-1,1))
nsamples2 <- min(2^14,nsamples)
alpha <- sample(rep(alpha0,2),nsamples2,replace=T)
q <- extraDistr::rdirichlet(n=nsamples2,alpha=matrix(alpha/nclusters,nsamples2,nclusters))
sd <- sample(rep(sqrt(rvar0),2),nsamples2*nclusters,replace=T)
m <- matrix(rnorm(nsamples2*nclusters,rmean0,sd),nsamples2)
shapein <- sample(rep(rshapein0,2),nsamples2*nclusters,replace=T)
shapeout <- sample(rep(rshapeout0,2),nsamples2*nclusters,replace=T)
scalevar <- sample(rep(rvarscales,2),nsamples2*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples2*nclusters,shape=shapeout,rate=nimble::rinvgamma(nsamples2*nclusters,shape=shapein,rate=scalevar))),nsamples2)
##
xgrid <- seq(xmin, xmax, length.out=256)
txgrid <- tran(xgrid)
ysum <- 0
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
par(mfrow=rowcol,mar = c(0,0,0,0))
for(i in 1:nsamples2){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        q[i,acluster] *
            dnorm(txgrid, m[i,acluster], s[i,acluster]) * jac(txgrid)
    }))
    ysum <- ysum+y
    if(i < prod(rowcol) | i==nsamples2){
    if(i == nsamples2){y <- ysum/nsamples2}
    ## if(!is.null(data)){
    ##     his <- thist(data)
    ##     ymax <- max(y,his$density)
    ## }else{ymax <- NULL}
    ymax <- NULL
    tplot(x=xgrid, y=y,
          ylim=c(0,max(y,ymax)),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=(if(i < prod(rowcol)){1}else{if(any(is.infinite(ysum))){2}else{3}}), ly=1,lwd=(if(i < prod(rowcol)){0.5}else{1}))
    ## tplot(x=xgrid[extr2], y=y[extr2],
    ##       type='p',col=4,cex=0.15,add=T,pch=3)
    ## if(!is.null(data)){
    ##     tplot(x=his$mids,y=his$density,type='l',lwd=0.5,add=T,alpha=0.25,col=4)
    ## }
    abline(h=c(0),lwd=0.5,col=alpha2hex2(0.5,c(7,2)),lty=c(1,2))
##    if(i==nsamples2){
    if(TRUE){
        abline(v=c(xq1,xq3),lwd=0.5,col=alpha2hex(7,0.5),lty=1)
        ## abline(v=c(xmin,xmax),lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
    }
}
dev.off()




#### Plot of mixture of Gaussians
nclusters <- 64
alpha0 <- 2^((-3):3)
rmean0 <- 0
rvar0 <- (1)^2
rshapein0 <- 1 # large scales
rshapeout0 <- 1 # small scales
hwidth <- 2 # number of powers of 2 to consider in either direction
rvarscales <- (1 * 2^((-hwidth):hwidth))^2
##
nsamples <- 2^24
xsamples <- rnorm(nsamples,
                  mean=rnorm(nsamples,mean=rmean0,sd=sqrt(rvar0)),
                  sd=sqrt(
                      extraDistr::rbetapr(nsamples,shape1=rshapein0,shape2=rshapeout0,
                                          scale=sample(rvarscales,nsamples,replace=T))
                  )
                  )

thismad <- mad(xsamples,constant=1)
thismad
xr <- ceiling(max(abs(xsamples)))
xgrid <- seq(-xr,xr,by=0.1)
his <- thist(xsamples,n=xgrid)
his$density <- (his$density + rev(his$density))/2
pdff('Gaussmix')
tplot(x=his$mids, y=list(his$density,dnorm(his$mids,sd=thismad/qnorm(3/4)),dcauchy(his$mids,scale=thismad),dlogis(his$mids,scale=thismad/qlogis(3/4))),
      xlim=c(-1,1)*6,
      xticks=seq(-xr,xr,by=1), xlabels=seq(-xr,xr,by=1),
                                        #sapply(seq(-3,3,by=1),function(i)as.expression(bquote(.(i*2)*bar(sigma))))
      xlab=expression(italic(x)/bar(sigma)), ylab='density',
      lwd=c(3,2,2,5),lty=c(1,2,4,3), alpha=c(0,rep(0.25,3)),
      mar=c(NA,5,2,1))
abline(v=c(-1,1)*thismad,col=alpha2hex(7,0.25),lwd=2)
dev.off()


#### Calculate and save function Q
nint <- 512
seqnint <- (1:(nint-1))/nint # seq(1/nint,(nint-1)/nint,length.out=nint-1)
nsamplesx <- 2^24
testgr <- c(NULL,
            tquant(rnorm(nsamplesx,
                  mean=rnorm(nsamplesx,mean=rmean0,sd=sqrt(rvar0)),
                  sd=sqrt(
                      extraDistr::rbetapr(nsamplesx,shape1=rshapein0,shape2=rshapeout0,
                                          scale=sample(rvarscales,nsamplesx,replace=T))
                  )
                  ),
                  seqnint
                  ),
            NULL)
testgr <- (testgr-rev(testgr))/2
approxq <- approxfun(x=seqnint, y=testgr, yleft=-Inf, yright=+Inf)
saveRDS(approxq,'Qfunction512.rsd')

xss <- foreach(nint=rev(c(5,
                10,32,100,
                256)))%do%{
                (1:(nint-1))/nint}

nint <- 256
xgrid <- seq(1/nint,(nint-1)/nint,length.out=nint-1)
pdff('Qfunction2')
tplot(x=xgrid,y=list(approxq(xgrid),qnorm(xgrid,sd=thismad/qnorm(3/4)),qcauchy(xgrid,scale=thismad),qlogis(xgrid,scale=thismad/qlogis(3/4))),
      lwd=c(3,2,2,5),lty=c(1,2,4,3), alpha=c(0,rep(0.25,3)),
      ylim=range(approxq(xgrid)), 
      ## xticks=c(0,0.25,0.5,0.75,1),xlabels=c(0,expression(italic(m)/4),expression(italic(m)/2),expression(3*italic(m)/4),expression(italic(m))),
      xlab=expression(italic(x)), ylab=expression(italic(Q)(italic(x))),
      mar=c(NA,5,2,1))
dev.off()




#### doubly-censored variate
sdovermad2 <- 0.5/qnorm(0.75)
## dt <- fread('~/repositories/ledley-jaynes_machine/scripts/ingrid_data_nogds6.csv')
## varinfo <- data.matrix(read.csv('~/repositories/ledley-jaynes_machine/scripts/varinfo.csv',row.names=1))
graphics.off()
pdff('priorsamples_2censored')
set.seed(123)
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
xmin <- -2
xmax <- 2
##
## hyperparameters
rowcol <- c(20,20)
nsamples <- 1e4
nclusters <- 64
alpha0 <- 2^((-3):3)
dmean0 <- 0
dvar0 <- (1)^2
dshapein0 <- 1 # large scales
dshapeout0 <- 1 # small scales
hwidth <- 2 # number of powers of 2 to consider in either direction
dvarscales <- (1 * 2^((-hwidth):hwidth))^2
##
alpha <- sample(rep(alpha0,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sqrt(dvar0),2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,dmean0,sd),nsamples)
shapein <- sample(rep(dshapein0,2),nsamples*nclusters,replace=T)
shapeout <- sample(rep(dshapeout0,2),nsamples*nclusters,replace=T)
scalevar <- sample(rep(dvarscales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shapeout,rate=nimble::rinvgamma(nsamples*nclusters,shape=shapein,rate=scalevar))),nsamples)
##
xgrid <- seq(xmin, xmax, length.out=256)
bgrid <- xgrid[c(1,length(xgrid))]
ysum <- ybsum <- 0
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
par(mfrow=rowcol,mar = c(0,0,0,0))
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        q[i,acluster] *
            dnorm(xgrid, m[i,acluster], s[i,acluster])
    }))
    yb <- rowSums(sapply(1:nclusters,function(acluster){
        q[i,acluster] *
            pnorm(bgrid*c(1,-1), c(1,-1)*m[i,acluster], s[i,acluster])
            ## c(pnorm(bgrid[1], c(1,-1)*m[i,acluster], s[i,acluster]),
            ##   pnorm(bgrid[2], m[i,acluster], s[i,acluster], lower.tail=F))
    }))
    ysum <- ysum+y
    ybsum <- ybsum+yb
    if(i < prod(rowcol) | i==nsamples){
        if(i == nsamples){
            y <- ysum/nsamples
        yb <- ybsum/nsamples}
    ## if(!is.null(data)){
    ##     his <- thist(data)
    ##     ymax <- max(y,his$density)
    ## }else{ymax <- NULL}
    ymax <- NULL
    tplot(x=xgrid, y=y,
          ylim=c(0,max(y,ymax)),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=(if(i < prod(rowcol)){1}else{if(any(is.infinite(ysum))){2}else{3}}),
          ly=1,lwd=(if(i < prod(rowcol)){0.5}else{1}))
    tplot(x=cbind(c(xmin,xmin),c(xmax,xmax)), y=rbind(0,yb)*max(y),
          type='l', lty=1, lwd=0.1,
          col=(if(i < prod(rowcol)){1}else{if(any(is.infinite(ysum))){2}else{3}}),
          add=T)
    tplot(x=c(xmin,xmax), y=yb*max(y),
          type='p', cex=0.2,
          col=(if(i < prod(rowcol)){1}else{if(any(is.infinite(ysum))){2}else{3}}),
          add=T)
    ## tplot(x=xgrid[extr2], y=y[extr2],
    ##       type='p',col=4,cex=0.15,add=T,pch=3)
    ## if(!is.null(data)){
    ##     tplot(x=his$mids,y=his$density,type='l',lwd=0.5,add=T,alpha=0.25,col=4)
    ## }
    abline(h=c(0,max(y)),lwd=0.5,col=alpha2hex2(0.5,c(7,7)),lty=c(1,2))
        abline(v=c(-1,1),lwd=0.5,col=alpha2hex(7,0.5),lty=1)
    ## if(i==nsamples){
    ##     abline(v=invtran(c(-1,1)),lwd=0.5,col=alpha2hex(2,0.5),lty=1)
    ##     ## abline(v=c(xmin,xmax),lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    ## }
    }
}
dev.off()


#### integer variate
approxq <- readRDS('Qfunction512.rsd')
sdovermad2 <- 0.5/qnorm(0.75)
## dt <- fread('~/repositories/ledley-jaynes_machine/scripts/ingrid_data_nogds6.csv')
## varinfo <- data.matrix(read.csv('~/repositories/ledley-jaynes_machine/scripts/varinfo.csv',row.names=1))
graphics.off()
pdff(paste0('priorsamples_integer_q'))
for(nint in c(5,
                10,32,100,
              256)){
set.seed(123)
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
nmin <- 1
nmax <- nint
dd <- (nint-2)/(2*nint^2)#
tsqrt <- function(x){sign(x)*sqrt(abs(x))}
tran <- function(x){qnorm((x-nmin)/(nmax-nmin)*(1-2*dd)+dd,mean=0,sd=transd)}
##
## hyperparameters
rowcol <- c(20,20)
nsamples <- 1e4
nclusters <- 64
alpha0 <- 2^((-3):3)
imean0 <- 0
ivar0 <- (1)^2#(7/8)^2
ishapein0 <- 1 # large scales
ishapeout0 <- 1 # small scales
hwidth <- 2 # number of powers of 2 to consider in either direction
ivarscales <- (1 * 2^((-hwidth):hwidth))^2
##
alpha <- sample(rep(alpha0,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sqrt(ivar0),2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,imean0,sd),nsamples)
shapein <- sample(rep(ishapein0,2),nsamples*nclusters,replace=T)
shapeout <- sample(rep(ishapeout0,2),nsamples*nclusters,replace=T)
scalevar <- sample(rep(ivarscales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shapeout,rate=nimble::rinvgamma(nsamples*nclusters,shape=shapein,rate=scalevar))),nsamples)
##
xgrid <- seq(nmin,nmax,length.out=nint)
transd <- 2*sdovermad2 #*log(nint)
mgrid <- approxq(seq(0,1,length.out=nint+1))
##mgrid <- qnorm(seq(0,1,length.out=nint+1), mean=0, sd=transd)
#mgrid <- c(-Inf, tsqrt((1:(nint-1))-nint/2)*sdovermad2/2, +Inf)
## kkk <- 1
## mgrid <- qnorm(c(0, seq(1/(kkk*nint), (kkk*nint-1)/(kkk*nint), length.out=nint-1), 1))
ysum <- 0
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
par(mfrow=rowcol,mar = c(0,0,0,0))
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        q[i,acluster] *
    diff(pnorm(mgrid, m[i,acluster], s[i,acluster]))
    }))
    ysum <- ysum+y
    if(i < prod(rowcol) | i==nsamples){
    if(i == nsamples){y <- ysum/nsamples}
    ## if(!is.null(data)){
    ##     his <- thist(data)
    ##     ymax <- max(y,his$density)
    ## }else{ymax <- NULL}
    ymax <- NULL
    tplot(x=xgrid, y=y,
          ylim=c(0,max(y,ymax)),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=(if(i < prod(rowcol)){1}else{if(any(is.infinite(ysum))){2}else{3}}), ly=1,lwd=(if(i < prod(rowcol)){0.5}else{1}))
    ## tplot(x=xgrid[extr2], y=y[extr2],
    ##       type='p',col=4,cex=0.15,add=T,pch=3)
    ## if(!is.null(data)){
    ##     tplot(x=his$mids,y=his$density,type='l',lwd=0.5,add=T,alpha=0.25,col=4)
    ## }
    abline(h=c(0),lwd=0.5,col=alpha2hex2(0.5,c(7,2)),lty=c(1,2))
    if(i==nsamples){text(nint/2,max(y)*0.9,nint,cex=0.5)}
    ## if(i==nsamples){
    ##     abline(v=invtran(c(-1,1)),lwd=0.5,col=alpha2hex(2,0.5),lty=1)
    ##     ## abline(v=c(xmin,xmax),lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    ## }
    }
}
}
dev.off()


#### integer variate - just one value of n
approxq <- readRDS('Qfunction512.rsd')
sdovermad2 <- 0.5/qnorm(0.75)
## dt <- fread('~/repositories/ledley-jaynes_machine/scripts/ingrid_data_nogds6.csv')
## varinfo <- data.matrix(read.csv('~/repositories/ledley-jaynes_machine/scripts/varinfo.csv',row.names=1))
graphics.off()
nint <- 32
set.seed(123)
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
nmin <- 1
nmax <- nint
dd <- (nint-2)/(2*nint^2)#
tsqrt <- function(x){sign(x)*sqrt(abs(x))}
tran <- function(x){qnorm((x-nmin)/(nmax-nmin)*(1-2*dd)+dd,mean=0,sd=transd)}
##
## hyperparameters
rowcol <- c(20,20)
nsamples <- 1e4
nclusters <- 64
alpha0 <- 2^((-3):3)
imean0 <- 0
ivar0 <- (1)^2#(7/8)^2
ishapein0 <- 1 # large scales
ishapeout0 <- 1 # small scales
hwidth <- 2 # number of powers of 2 to consider in either direction
ivarscales <- (1 * 2^((-hwidth):hwidth))^2
##
alpha <- sample(rep(alpha0,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
sd <- sample(rep(sqrt(ivar0),2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,imean0,sd),nsamples)
shapein <- sample(rep(ishapein0,2),nsamples*nclusters,replace=T)
shapeout <- sample(rep(ishapeout0,2),nsamples*nclusters,replace=T)
scalevar <- sample(rep(ivarscales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shapeout,rate=nimble::rinvgamma(nsamples*nclusters,shape=shapein,rate=scalevar))),nsamples)
##
xgrid <- seq(nmin,nmax,length.out=nint)
transd <- 2*sdovermad2 #*log(nint)
mgrid <- approxq(seq(0,1,length.out=nint+1))
##mgrid <- qnorm(seq(0,1,length.out=nint+1), mean=0, sd=transd)
#mgrid <- c(-Inf, tsqrt((1:(nint-1))-nint/2)*sdovermad2/2, +Inf)
## kkk <- 1
## mgrid <- qnorm(c(0, seq(1/(kkk*nint), (kkk*nint-1)/(kkk*nint), length.out=nint-1), 1))
ysum <- 0
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
pdff(paste0('priorsamples_integer_',nint))
par(mfrow=rowcol,mar = c(0,0,0,0))
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        q[i,acluster] *
    diff(pnorm(mgrid, m[i,acluster], s[i,acluster]))
    }))
    ysum <- ysum+y
    if(i < prod(rowcol) | i==nsamples){
    if(i == nsamples){y <- ysum/nsamples}
    ## if(!is.null(data)){
    ##     his <- thist(data)
    ##     ymax <- max(y,his$density)
    ## }else{ymax <- NULL}
    ymax <- NULL
    tplot(x=xgrid, y=y,
          ylim=c(0,max(y,ymax)),xlim=range(xgrid),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          col=(if(i < prod(rowcol)){1}else{if(any(is.infinite(ysum))){2}else{3}}), ly=1,lwd=(if(i < prod(rowcol)){0.5}else{1}))
    ## tplot(x=xgrid[extr2], y=y[extr2],
    ##       type='p',col=4,cex=0.15,add=T,pch=3)
    ## if(!is.null(data)){
    ##     tplot(x=his$mids,y=his$density,type='l',lwd=0.5,add=T,alpha=0.25,col=4)
    ## }
    abline(h=c(0),lwd=0.5,col=alpha2hex2(0.5,c(7,2)),lty=c(1,2))
    ## if(i==nsamples){
    ##     abline(v=invtran(c(-1,1)),lwd=0.5,col=alpha2hex(2,0.5),lty=1)
    ##     ## abline(v=c(xmin,xmax),lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    ## }
    }
}
dev.off()









#### continuous variate
sdovermad <- 1/qnorm(0.75)
sdovermad2 <- 0.5/qnorm(0.75)
## dt <- fread('~/repositories/ledley-jaynes_machine/scripts/ingrid_data_nogds6.csv')
## varinfo <- data.matrix(read.csv('~/repositories/ledley-jaynes_machine/scripts/varinfo.csv',row.names=1))
set.seed(123)
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
xlocation <- 0
xscale <- 1
xmin <- -5
xmax <- 5
##
## hyperparameters
rowcol <- c(20,20)
nsamples <- 1e6
nclusters <- 64
alpha0 <- 2^((-3):3)
rmean0 <- 0
zeta <- 1
rvar0 <- (zeta)^2
rshapein0 <- 1 # large scales
hwidth <- 2 # number of powers of 2 to consider in either direction
rvarscales <- (zeta * 2^((-hwidth):hwidth))^2
####
alpha <- sample(rep(alpha0,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
##
msh <- 2
msd <- sqrt(nimble::rinvgamma(nsamples, shape=msh, rate=1))
mm <- rnorm(nsamples, mean=0, sd=1)
m <- matrix(rnorm(nsamples*nclusters, mean=mm, sd=msd),nrow=nsamples)
##
## shapein <- 1 # sample(rep(rshapein0,2),nsamples*nclusters,replace=T)
swidth <- 2
scentre <- 1
rshapeout0 <- (scentre * 2^((-swidth):swidth)) # small scales
shapeout <- sample(rep(rshapeout0,2), nsamples, replace=T)
rshape <- 1
ratevar <- nimble::rinvgamma(nsamples, shape=rshape, rate=1)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters,shape=shapeout,rate=ratevar)),nrow=nsamples)
##
maxxsamples <- 2^10
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
xmax <- 0
exclu <- 2
nxsamples <- 2^10
graphics.off()
##pdff(paste0('priorsamples_testgammas_msh',msh,'_sc',scentre,'_rs',rshape))
pdff(paste0('priorsamples_testgammas'), apaper=3)
par(mfrow=rowcol,mar = c(0,0,0,0))
for(i in 1:prod(rowcol)){
    lab <- sample(1:nclusters, maxxsamples, prob=q[i,], replace=T)
    xgrid <- tquant(rnorm(maxxsamples, mean=m[i,lab], sd=s[i,lab]), c(exclu/2,100-exclu/2)/100)
    if(i < prod(rowcol)){
    xgrid <- seq(min(-1,xgrid[1]), max(1,xgrid[2]), length.out=256)
        y <- rowSums(sapply(1:nclusters,function(acluster){
            q[i,acluster] *
                dnorm(xgrid, m[i,acluster], s[i,acluster])
        }))
    }else{
        xgrid <- seq(-5,5,length.out=128)
        y <- sapply(xgrid, function(xx){
            mean(dnorm(xx, m[,1], s[,1]))
            })
        ## y <- foreach(j=1:nsamples, .combine='+', .inorder=F)%dopar%{
        ##     dnorm(xgrid, m[j,1], s[j,1])
        ## }/nsamples
    }
        tplot(x=xgrid, y=y,
              xlim=range(c(xgrid,-1,1)),
          ylim=c(0,NA),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          lwd=0.5,
          col=(if(i < prod(rowcol)){1}else{3}))
    ## tplot(x=xgrid[extr2], y=y[extr2],
    ##       type='p',col=4,cex=0.15,add=T,pch=3)
    ## if(!is.null(data)){
    ##     tplot(x=his$mids,y=his$density,type='l',lwd=0.5,add=T,alpha=0.25,col=4)
    ## }
    abline(h=c(0),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
    abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2,7)),lty=c(1,2,1))
}
for(i in 1:prod(rowcol)){
    labs <- sample(1:nclusters, nxsamples, prob=q[i,], replace=T)
    x <- rnorm(nxsamples, mean=m[i,labs], sd=s[i,labs])
    y <- rnorm(nxsamples, mean=m[i+round(nsamples/2),labs], sd=s[i+round(nsamples/2),labs])
            tplot(x=x,y=y,type='p',pch='.',alpha=0.75,
                  xlim=range(tquant(x, c(exclu/2,100-exclu/2)/100),-1,1),
                  ylim=range(tquant(y, c(exclu/2,100-exclu/2)/100),-1,1),
                  xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
                  xticks=NA,yticks=NA,
                  mar=c(1,1,1,1)*0.5)
            abline(h=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=par('usr')[1:2],col='black',lwd=0.5)
            abline(h=par('usr')[3:4],col='black',lwd=0.5)
}
dev.off()

#### continuous variate - betaprime, no hyperpriors
sdovermad <- 1/qnorm(0.75)
sdovermad2 <- 0.5/qnorm(0.75)
## dt <- fread('~/repositories/ledley-jaynes_machine/scripts/ingrid_data_nogds6.csv')
## varinfo <- data.matrix(read.csv('~/repositories/ledley-jaynes_machine/scripts/varinfo.csv',row.names=1))
set.seed(123)
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
xlocation <- 0
xscale <- 1
xmin <- -5
xmax <- 5
##
## hyperparameters
rowcol <- c(20,20)
nsamples <- 1e6
nclusters <- 64
rmean0 <- 0
zeta <- 1
rvar0 <- (zeta)^2
rshapein0 <- 1 # large scales
hwidth <- 2 # number of powers of 2 to consider in either direction
rvarscales <- (zeta * 2^((-hwidth):hwidth))^2
#### Weights
alphas <- sample(rep(2^((-3):3), 2), nsamples, replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alphas/nclusters,nsamples,nclusters))
#### Means
meansHshShape <- 2
## meansHsd <- sample(rep((1 * 2^((-2):2)), 2), nsamples, replace=T)
meansHsd <- 1 # sqrt(nimble::rinvgamma(nsamples, shape=meansHshShape, rate=nimble::rinvgamma(nsamples, shape=meansHshShape, rate=1)))
meansHmean <- 0 # rnorm(nsamples, mean=0, sd=1)
means <- matrix(rnorm(nsamples*nclusters, mean=meansHmean, sd=meansHsd), nrow=nsamples)
#### Sds
sdsHshapes <- 1 # sample(rep((1 * 2^((-2):2)), 2), nsamples, replace=T)
sdsHrates <- 1 # sample(rep((1 * 2^((-2):2)), 2), nsamples, replace=T)
sds <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters, shape=sdsHshapes, rate=nimble::rinvgamma(nsamples*nclusters, shape=sdsHshapes, rate=sdsHrates))), nrow=nsamples)
##
maxxsamples <- 2^10
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
xmax <- 0
exclu <- 2
nxsamples <- 2^10
graphics.off()
##pdff(paste0('priorsamples_testgammas_msh',msh,'_sc',scentre,'_rs',rshape))
pdff(paste0('priorsamples_testbetaprimes_nohyper'), apaper=3)
par(mfrow=rowcol,mar = c(0,0,0,0))
for(i in 1:prod(rowcol)){
    lab <- sample(1:nclusters, maxxsamples, prob=q[i,], replace=T)
    xgrid <- tquant(rnorm(maxxsamples, mean=means[i,lab], sd=sds[i,lab]), c(exclu/2,100-exclu/2)/100)
    if(i < prod(rowcol)){
    xgrid <- seq(min(-1,xgrid[1]), max(1,xgrid[2]), length.out=256)
        y <- rowSums(sapply(1:nclusters,function(acluster){
            q[i,acluster] *
                dnorm(xgrid, means[i,acluster], sds[i,acluster])
        }))
    }else{
        xgrid <- seq(-5,5,length.out=128)
        y <- sapply(xgrid, function(xx){
            mean(dnorm(xx, means[,1], sds[,1]))
            })
        ## y <- foreach(j=1:nsamples, .combine='+', .inorder=F)%dopar%{
        ##     dnorm(xgrid, means[j,1], sds[j,1])
        ## }/nsamples
    }
        tplot(x=xgrid, y=y,
              xlim=range(c(xgrid,-1,1)),
          ylim=c(0,NA),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          lwd=0.5,
          col=(if(i < prod(rowcol)){1}else{3}))
    ## tplot(x=xgrid[extr2], y=y[extr2],
    ##       type='p',col=4,cex=0.15,add=T,pch=3)
    ## if(!is.null(data)){
    ##     tplot(x=his$mids,y=his$density,type='l',lwd=0.5,add=T,alpha=0.25,col=4)
    ## }
    abline(h=c(0),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
    abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2,7)),lty=c(1,2,1))
}
for(i in 1:prod(rowcol)){
    labs <- sample(1:nclusters, nxsamples, prob=q[i,], replace=T)
    x <- rnorm(nxsamples, mean=means[i,labs], sd=sds[i,labs])
    y <- rnorm(nxsamples, mean=means[i+round(nsamples/2),labs], sd=sds[i+round(nsamples/2),labs])
            tplot(x=x,y=y,type='p',pch='.',alpha=0.75,
                  xlim=range(tquant(x, c(exclu/2,100-exclu/2)/100),-1,1),
                  ylim=range(tquant(y, c(exclu/2,100-exclu/2)/100),-1,1),
                  xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
                  xticks=NA,yticks=NA,
                  mar=c(1,1,1,1)*0.5)
            abline(h=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=par('usr')[1:2],col='black',lwd=0.5)
            abline(h=par('usr')[3:4],col='black',lwd=0.5)
}
dev.off()

#### continuous variate - betaprime, discrete scale hyperpriors
sdovermad <- 1/qnorm(0.75)
sdovermad2 <- 0.5/qnorm(0.75)
## dt <- fread('~/repositories/ledley-jaynes_machine/scripts/ingrid_data_nogds6.csv')
## varinfo <- data.matrix(read.csv('~/repositories/ledley-jaynes_machine/scripts/varinfo.csv',row.names=1))
set.seed(123)
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
xlocation <- 0
xscale <- 1
xmin <- -5
xmax <- 5
##
## hyperparameters
rowcol <- c(20,20)
nsamples <- 1e6
nclusters <- 64
rmean0 <- 0
zeta <- 1
rvar0 <- (zeta)^2
rshapein0 <- 1 # large scales
hwidth <- 2 # number of powers of 2 to consider in either direction
rvarscales <- (zeta * 2^((-hwidth):hwidth))^2
#### Weights
alphas <- sample(rep(2^((-3):3), 2), nsamples, replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alphas/nclusters,nsamples,nclusters))
#### Means
meansHshShape <- 2
meansHsd <- sample(rep((1 * 2^((-2):2)), 2), nsamples, replace=T) # sqrt(nimble::rinvgamma(nsamples, shape=meansHshShape, rate=nimble::rinvgamma(nsamples, shape=meansHshShape, rate=1)))
meansHmean <- 0 # rnorm(nsamples, mean=0, sd=1)
means <- matrix(rnorm(nsamples*nclusters, mean=meansHmean, sd=meansHsd),
                nrow=nsamples)
#### Sds
sdsHshapes <- 1 # sample(rep((1 * 2^((-2):2)), 2), nsamples, replace=T)
sdsHrates <- sample(rep((1 * 2^((-2):2))^2, 2), nsamples, replace=T)
sds <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters, shape=sdsHshapes,
                                     rate=nimble::rinvgamma(nsamples*nclusters, shape=sdsHshapes, rate=sdsHrates))),
              nrow=nsamples)
##
maxxsamples <- 2^10
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
xmax <- 0
exclu <- 2
nxsamples <- 2^10
graphics.off()
##pdff(paste0('priorsamples_testgammas_msh',msh,'_sc',scentre,'_rs',rshape))
pdff(paste0('priorsamples_testbetaprimes_discrscalehyper'), apaper=3)
par(mfrow=rowcol,mar = c(0,0,0,0))
for(i in 1:prod(rowcol)){
    lab <- sample(1:nclusters, maxxsamples, prob=q[i,], replace=T)
    xgrid <- tquant(rnorm(maxxsamples, mean=means[i,lab], sd=sds[i,lab]), c(exclu/2,100-exclu/2)/100)
    if(i < prod(rowcol)){
    xgrid <- seq(min(-1,xgrid[1]), max(1,xgrid[2]), length.out=256)
        y <- rowSums(sapply(1:nclusters,function(acluster){
            q[i,acluster] *
                dnorm(xgrid, means[i,acluster], sds[i,acluster])
        }))
    }else{
        xgrid <- seq(-5,5,length.out=128)
        y <- sapply(xgrid, function(xx){
            mean(dnorm(xx, means[,1], sds[,1]))
            })
        ## y <- foreach(j=1:nsamples, .combine='+', .inorder=F)%dopar%{
        ##     dnorm(xgrid, means[j,1], sds[j,1])
        ## }/nsamples
    }
        tplot(x=xgrid, y=y,
              xlim=range(c(xgrid,-1,1)),
          ylim=c(0,NA),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          lwd=0.5,
          col=(if(i < prod(rowcol)){1}else{3}))
    ## tplot(x=xgrid[extr2], y=y[extr2],
    ##       type='p',col=4,cex=0.15,add=T,pch=3)
    ## if(!is.null(data)){
    ##     tplot(x=his$mids,y=his$density,type='l',lwd=0.5,add=T,alpha=0.25,col=4)
    ## }
    abline(h=c(0),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
    abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2,7)),lty=c(1,2,1))
}
for(i in 1:prod(rowcol)){
    labs <- sample(1:nclusters, nxsamples, prob=q[i,], replace=T)
    x <- rnorm(nxsamples, mean=means[i,labs], sd=sds[i,labs])
    y <- rnorm(nxsamples, mean=means[i+round(nsamples/2),labs], sd=sds[i+round(nsamples/2),labs])
            tplot(x=x,y=y,type='p',pch='.',alpha=0.75,
                  xlim=range(tquant(x, c(exclu/2,100-exclu/2)/100),-1,1),
                  ylim=range(tquant(y, c(exclu/2,100-exclu/2)/100),-1,1),
                  xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
                  xticks=NA,yticks=NA,
                  mar=c(1,1,1,1)*0.5)
            abline(h=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=par('usr')[1:2],col='black',lwd=0.5)
            abline(h=par('usr')[3:4],col='black',lwd=0.5)
}
dev.off()

#### continuous variate - betaprime, discrete shape hyperpriors
sdovermad <- 1/qnorm(0.75)
sdovermad2 <- 0.5/qnorm(0.75)
## dt <- fread('~/repositories/ledley-jaynes_machine/scripts/ingrid_data_nogds6.csv')
## varinfo <- data.matrix(read.csv('~/repositories/ledley-jaynes_machine/scripts/varinfo.csv',row.names=1))
set.seed(123)
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
xlocation <- 0
xscale <- 1
xmin <- -5
xmax <- 5
##
## hyperparameters
rowcol <- c(20,20)
nsamples <- 1e6
nclusters <- 64
rmean0 <- 0
zeta <- 1
rvar0 <- (zeta)^2
rshapein0 <- 1 # large scales
hwidth <- 2 # number of powers of 2 to consider in either direction
rvarscales <- (zeta * 2^((-hwidth):hwidth))^2
#### Weights
alphas <- sample(rep(2^((-3):3), 2), nsamples, replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alphas/nclusters,nsamples,nclusters))
#### Means
meansHshShape <- 2
meansHsd <- sample(rep((1 * 2^((-2):2)), 2), nsamples, replace=T) # sqrt(nimble::rinvgamma(nsamples, shape=meansHshShape, rate=nimble::rinvgamma(nsamples, shape=meansHshShape, rate=1)))
meansHmean <- 0 # rnorm(nsamples, mean=0, sd=1)
means <- matrix(rnorm(nsamples*nclusters, mean=meansHmean, sd=meansHsd),
                nrow=nsamples)
#### Sds
sdsHshapes <- sample(rep((1 * 2^((-2):2)), 2), nsamples, replace=T)
sdsHrates <- 1 # sample(rep((1 * 2^((-2):2))^2, 2), nsamples, replace=T)
sds <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters, shape=sdsHshapes,
                                     rate=nimble::rinvgamma(nsamples*nclusters, shape=sdsHshapes, rate=sdsHrates))),
              nrow=nsamples)
##
maxxsamples <- 2^10
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
xmax <- 0
exclu <- 2
nxsamples <- 2^10
graphics.off()
##pdff(paste0('priorsamples_testgammas_msh',msh,'_sc',scentre,'_rs',rshape))
pdff(paste0('priorsamples_testbetaprimes_discrshapehyper'), apaper=3)
par(mfrow=rowcol,mar = c(0,0,0,0))
for(i in 1:prod(rowcol)){
    lab <- sample(1:nclusters, maxxsamples, prob=q[i,], replace=T)
    xgrid <- tquant(rnorm(maxxsamples, mean=means[i,lab], sd=sds[i,lab]), c(exclu/2,100-exclu/2)/100)
    if(i < prod(rowcol)){
    xgrid <- seq(min(-1,xgrid[1]), max(1,xgrid[2]), length.out=256)
        y <- rowSums(sapply(1:nclusters,function(acluster){
            q[i,acluster] *
                dnorm(xgrid, means[i,acluster], sds[i,acluster])
        }))
    }else{
        xgrid <- seq(-5,5,length.out=128)
        y <- sapply(xgrid, function(xx){
            mean(dnorm(xx, means[,1], sds[,1]))
            })
        ## y <- foreach(j=1:nsamples, .combine='+', .inorder=F)%dopar%{
        ##     dnorm(xgrid, means[j,1], sds[j,1])
        ## }/nsamples
    }
        tplot(x=xgrid, y=y,
              xlim=range(c(xgrid,-1,1)),
          ylim=c(0,NA),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          lwd=0.5,
          col=(if(i < prod(rowcol)){1}else{3}))
    ## tplot(x=xgrid[extr2], y=y[extr2],
    ##       type='p',col=4,cex=0.15,add=T,pch=3)
    ## if(!is.null(data)){
    ##     tplot(x=his$mids,y=his$density,type='l',lwd=0.5,add=T,alpha=0.25,col=4)
    ## }
    abline(h=c(0),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
    abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2,7)),lty=c(1,2,1))
}
for(i in 1:prod(rowcol)){
    labs <- sample(1:nclusters, nxsamples, prob=q[i,], replace=T)
    x <- rnorm(nxsamples, mean=means[i,labs], sd=sds[i,labs])
    y <- rnorm(nxsamples, mean=means[i+round(nsamples/2),labs], sd=sds[i+round(nsamples/2),labs])
            tplot(x=x,y=y,type='p',pch='.',alpha=0.75,
                  xlim=range(tquant(x, c(exclu/2,100-exclu/2)/100),-1,1),
                  ylim=range(tquant(y, c(exclu/2,100-exclu/2)/100),-1,1),
                  xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
                  xticks=NA,yticks=NA,
                  mar=c(1,1,1,1)*0.5)
            abline(h=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=par('usr')[1:2],col='black',lwd=0.5)
            abline(h=par('usr')[3:4],col='black',lwd=0.5)
}
dev.off()

#### continuous variate - betaprime
sdovermad <- 1/qnorm(0.75)
sdovermad2 <- 0.5/qnorm(0.75)
## dt <- fread('~/repositories/ledley-jaynes_machine/scripts/ingrid_data_nogds6.csv')
## varinfo <- data.matrix(read.csv('~/repositories/ledley-jaynes_machine/scripts/varinfo.csv',row.names=1))
set.seed(123)
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
xlocation <- 0
xscale <- 1
xmin <- -5
xmax <- 5
##
## hyperparameters
rowcol <- c(20,20)
nsamples <- 1e6
nclusters <- 64
rmean0 <- 0
zeta <- 1
rvar0 <- (zeta)^2
rshapein0 <- 1 # large scales
hwidth <- 2 # number of powers of 2 to consider in either direction
rvarscales <- (zeta * 2^((-hwidth):hwidth))^2
#### Weights
alphas <- sample(rep(2^((-3):3), 2), nsamples, replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alphas/nclusters,nsamples,nclusters))
#### Means
meansHshShape <- 2
## meansHsd <- sample(rep((1 * 2^((-2):2)), 2), nsamples, replace=T)
meansHsd <- sqrt(nimble::rinvgamma(nsamples, shape=meansHshShape,
                                   rate=nimble::rinvgamma(nsamples, shape=meansHshShape, rate=1)))
meansHmean <- rnorm(nsamples, mean=0, sd=1)
means <- matrix(rnorm(nsamples*nclusters, mean=meansHmean, sd=meansHsd),
                nrow=nsamples)
#### Sds
sdsHshapes <- sample(rep((1 * 2^((-2):2)), 2), nsamples, replace=T)
sdsHrates <- sample(rep((1 * 2^((-2):2)), 2), nsamples, replace=T)
sds <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters, shape=sdsHshapes,
                                     rate=nimble::rinvgamma(nsamples*nclusters, shape=sdsHshapes, rate=sdsHrates))),
              nrow=nsamples)
##
maxxsamples <- 2^10
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
xmax <- 0
exclu <- 2
nxsamples <- 2^10
graphics.off()
##pdff(paste0('priorsamples_testgammas_msh',msh,'_sc',scentre,'_rs',rshape))
pdff(paste0('priorsamples_testbetaprimes'), apaper=3)
par(mfrow=rowcol,mar = c(0,0,0,0))
for(i in 1:prod(rowcol)){
    lab <- sample(1:nclusters, maxxsamples, prob=q[i,], replace=T)
    xgrid <- tquant(rnorm(maxxsamples, mean=means[i,lab], sd=sds[i,lab]), c(exclu/2,100-exclu/2)/100)
    if(i < prod(rowcol)){
    xgrid <- seq(min(-1,xgrid[1]), max(1,xgrid[2]), length.out=256)
        y <- rowSums(sapply(1:nclusters,function(acluster){
            q[i,acluster] *
                dnorm(xgrid, means[i,acluster], sds[i,acluster])
        }))
    }else{
        xgrid <- seq(-5,5,length.out=128)
        y <- sapply(xgrid, function(xx){
            mean(dnorm(xx, means[,1], sds[,1]))
            })
        ## y <- foreach(j=1:nsamples, .combine='+', .inorder=F)%dopar%{
        ##     dnorm(xgrid, means[j,1], sds[j,1])
        ## }/nsamples
    }
        tplot(x=xgrid, y=y,
              xlim=range(c(xgrid,-1,1)),
          ylim=c(0,NA),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          lwd=0.5,
          col=(if(i < prod(rowcol)){1}else{3}))
    ## tplot(x=xgrid[extr2], y=y[extr2],
    ##       type='p',col=4,cex=0.15,add=T,pch=3)
    ## if(!is.null(data)){
    ##     tplot(x=his$mids,y=his$density,type='l',lwd=0.5,add=T,alpha=0.25,col=4)
    ## }
    abline(h=c(0),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
    abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2,7)),lty=c(1,2,1))
}
for(i in 1:prod(rowcol)){
    labs <- sample(1:nclusters, nxsamples, prob=q[i,], replace=T)
    x <- rnorm(nxsamples, mean=means[i,labs], sd=sds[i,labs])
    y <- rnorm(nxsamples, mean=means[i+round(nsamples/2),labs], sd=sds[i+round(nsamples/2),labs])
            tplot(x=x,y=y,type='p',pch='.',alpha=0.75,
                  xlim=range(tquant(x, c(exclu/2,100-exclu/2)/100),-1,1),
                  ylim=range(tquant(y, c(exclu/2,100-exclu/2)/100),-1,1),
                  xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
                  xticks=NA,yticks=NA,
                  mar=c(1,1,1,1)*0.5)
            abline(h=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=par('usr')[1:2],col='black',lwd=0.5)
            abline(h=par('usr')[3:4],col='black',lwd=0.5)
}
dev.off()

#### continuous variate - double betaprime
sdovermad <- 1/qnorm(0.75)
sdovermad2 <- 0.5/qnorm(0.75)
## dt <- fread('~/repositories/ledley-jaynes_machine/scripts/ingrid_data_nogds6.csv')
## varinfo <- data.matrix(read.csv('~/repositories/ledley-jaynes_machine/scripts/varinfo.csv',row.names=1))
set.seed(123)
## tran <- function(x){qnorm(x*(1-2*dd)+dd)}
## jac <- function(x){1/dnorm(x*(1-2*dd)+dd)*(1-2*dd)}
xlocation <- 0
xscale <- 1
xmin <- -5
xmax <- 5
##
## hyperparameters
rowcol <- c(20,20)
nsamples <- 1e6
nclusters <- 64
rmean0 <- 0
zeta <- 1
rvar0 <- (zeta)^2
rshapein0 <- 1 # large scales
hwidth <- 2 # number of powers of 2 to consider in either direction
rvarscales <- (zeta * 2^((-hwidth):hwidth))^2
#### Weights
alphas <- sample(rep(2^((-3):3), 2), nsamples, replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alphas/nclusters,nsamples,nclusters))
#### Means
meansHshShape <- 2
## meansHsd <- sample(rep((1 * 2^((-2):2)), 2), nsamples, replace=T)
meansHsd <- sqrt(nimble::rinvgamma(nsamples, shape=meansHshShape,
                                   rate=nimble::rinvgamma(nsamples, shape=meansHshShape, rate=1)))
meansHmean <- rnorm(nsamples, mean=0, sd=1)
means <- matrix(rnorm(nsamples*nclusters, mean=meansHmean, sd=meansHsd),
                nrow=nsamples)
#### Sds
sdsHshapes <- sample(rep((1 * 2^((-2):2)), 2), nsamples, replace=T)
sdsHratesShape <- 2
sdsHrates <- nimble::rinvgamma(nsamples, shape=sdsHratesShape,
                                     rate=nimble::rinvgamma(nsamples, shape=sdsHratesShape, rate=1))
sds <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters, shape=sdsHshapes,
                                     rate=nimble::rinvgamma(nsamples*nclusters, shape=sdsHshapes, rate=sdsHrates))),
              nrow=nsamples)
##
maxxsamples <- 2^10
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
xmax <- 0
exclu <- 2
nxsamples <- 2^10
graphics.off()
##pdff(paste0('priorsamples_testgammas_msh',msh,'_sc',scentre,'_rs',rshape))
pdff(paste0('priorsamples_testbetaprimesdouble'), apaper=3)
par(mfrow=rowcol,mar = c(0,0,0,0))
for(i in 1:prod(rowcol)){
    lab <- sample(1:nclusters, maxxsamples, prob=q[i,], replace=T)
    xgrid <- tquant(rnorm(maxxsamples, mean=means[i,lab], sd=sds[i,lab]), c(exclu/2,100-exclu/2)/100)
    if(i < prod(rowcol)){
    xgrid <- seq(min(-1,xgrid[1]), max(1,xgrid[2]), length.out=256)
        y <- rowSums(sapply(1:nclusters,function(acluster){
            q[i,acluster] *
                dnorm(xgrid, means[i,acluster], sds[i,acluster])
        }))
    }else{
        xgrid <- seq(-5,5,length.out=128)
        y <- sapply(xgrid, function(xx){
            mean(dnorm(xx, means[,1], sds[,1]))
            })
        ## y <- foreach(j=1:nsamples, .combine='+', .inorder=F)%dopar%{
        ##     dnorm(xgrid, means[j,1], sds[j,1])
        ## }/nsamples
    }
        tplot(x=xgrid, y=y,
              xlim=range(c(xgrid,-1,1)),
          ylim=c(0,NA),
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5,
          lwd=0.5,
          col=(if(i < prod(rowcol)){1}else{3}))
    ## tplot(x=xgrid[extr2], y=y[extr2],
    ##       type='p',col=4,cex=0.15,add=T,pch=3)
    ## if(!is.null(data)){
    ##     tplot(x=his$mids,y=his$density,type='l',lwd=0.5,add=T,alpha=0.25,col=4)
    ## }
    abline(h=c(0),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
    abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2,7)),lty=c(1,2,1))
}
for(i in 1:prod(rowcol)){
    labs <- sample(1:nclusters, nxsamples, prob=q[i,], replace=T)
    x <- rnorm(nxsamples, mean=means[i,labs], sd=sds[i,labs])
    y <- rnorm(nxsamples, mean=means[i+round(nsamples/2),labs], sd=sds[i+round(nsamples/2),labs])
            tplot(x=x,y=y,type='p',pch='.',alpha=0.75,
                  xlim=range(tquant(x, c(exclu/2,100-exclu/2)/100),-1,1),
                  ylim=range(tquant(y, c(exclu/2,100-exclu/2)/100),-1,1),
                  xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
                  xticks=NA,yticks=NA,
                  mar=c(1,1,1,1)*0.5)
            abline(h=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=par('usr')[1:2],col='black',lwd=0.5)
            abline(h=par('usr')[3:4],col='black',lwd=0.5)
}
dev.off()


## goerueretal2010:
## W(s|b,(bw)^-1) = s^(b/2-1) exp[-s/2 * wb]
## G(s|b,a) = s^(b/2-1) exp[-s/2 * b/a]


pdff('justtest_gamma', apaper=3)
nsamples <- 2^17
xlim <- c(-5,5)
##
m <- rnorm(nsamples,mean=0,sd=1)
v <- nimble::rinvgamma(nsamples,1,1)
x <- rnorm(nsamples,mean=m,sd=sqrt(v))
y <- rnorm(nsamples,mean=m,sd=sqrt(v))
tplot(x,y,type='p', pch='.', alpha=0,
      xlim=range(c(tquant(x,c(2.5,97.5)/100),xlim)),
      ylim=range(c(tquant(y,c(2.5,97.5)/100),xlim))
      )
##
v <- nimble::rinvgamma(nsamples,1,1)
x <- rnorm(nsamples,mean=0,sd=sqrt(v))
y <- rnorm(nsamples,mean=0,sd=sqrt(v))
tplot(x,y,type='p', pch='.', alpha=0,
      xlim=range(c(tquant(x,c(2.5,97.5)/100),xlim)),
      ylim=range(c(tquant(y,c(2.5,97.5)/100),xlim))
      )
##
vx <- nimble::rinvgamma(nsamples,1,1)
vy <- nimble::rinvgamma(nsamples,1,1)
x <- rnorm(nsamples,mean=0,sd=sqrt(vx))
y <- rnorm(nsamples,mean=0,sd=sqrt(vy))
tplot(x,y,type='p', pch='.', alpha=0,
      xlim=range(c(tquant(x,c(2.5,97.5)/100),xlim)),
      ylim=range(c(tquant(y,c(2.5,97.5)/100),xlim))
      )
dev.off()
