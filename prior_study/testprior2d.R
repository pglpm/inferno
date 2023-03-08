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

#### Currently used ####

#### continuous variate per cluster
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
rowcol <- c(20,20)
nsamples <- 1e6
## hyperparameters
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
##
nsamples <- 400
alpha <- sample(rep(alpha0,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
##
sd <- sample(rep(sqrt(rvar0),2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,rmean0,sd),nsamples)
##
shapein <- rshapein0
shapeout <- rshapeout0
scalevar <- sample(rep(rvarscales,2),nsamples*nclusters,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters, shape=shapeout,
                                   rate=nimble::rinvgamma(nsamples*nclusters, shape=shapein, rate=scalevar))),
            nsamples)
##
xgrid <- seq(xmin, xmax, length.out=256)
ysum <- 0
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
graphics.off()
pdff('priorsamples_real_percluster')
par(mfrow=rowcol,mar = c(0,0,0,0))
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        q[i,acluster] *
            dnorm(xgrid, m[i,acluster], s[i,acluster])
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
##    if(i==nsamples2){
    if(TRUE){
        abline(v=c(-1,1),lwd=0.5,col=alpha2hex(7,0.5),lty=1)
        ## abline(v=c(xmin,xmax),lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
    }
}
dev.off()



#### continuous variate common across clusters
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
rowcol <- c(20,20)
nsamples <- 1e6
## hyperparameters
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
##
nsamples <- 400
alpha <- sample(rep(alpha0,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
##
sd <- sample(rep(sqrt(rvar0),2),nsamples*nclusters,replace=T)
m <- matrix(rnorm(nsamples*nclusters,rmean0,sd),nsamples)
##
shapein <- rshapein0
shapeout <- rshapeout0
scalevar <- sample(rep(rvarscales,2),nsamples,replace=T)
s <- matrix(sqrt(nimble::rinvgamma(nsamples*nclusters, shape=shapeout,
                                   rate=nimble::rinvgamma(nsamples, shape=shapein, rate=scalevar))),
            nsamples)
##
xgrid <- seq(xmin, xmax, length.out=256)
ysum <- 0
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
graphics.off()
pdff('priorsamples_real_commoncluster')
par(mfrow=rowcol,mar = c(0,0,0,0))
for(i in 1:nsamples){
    y <- rowSums(sapply(1:nclusters,function(acluster){
        q[i,acluster] *
            dnorm(xgrid, m[i,acluster], s[i,acluster])
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
##    if(i==nsamples2){
    if(TRUE){
        abline(v=c(-1,1),lwd=0.5,col=alpha2hex(7,0.5),lty=1)
        ## abline(v=c(xmin,xmax),lwd=0.5,col=alpha2hex(0.5,7),lty=2)
    }
    }
}
dev.off()




#### 2D

#### continuous variate per cluster
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
rowcol <- c(10,10)
nsamples <- 1e6
## hyperparameters
nclusters <- 64
alpha0 <- 2^((-3):3)
rmean0 <- 0
zeta <- 1
rvar0 <- (2*zeta)^2
rshapein0 <- 1 # large scales
rshapeout0 <- 1 # small scales
hwidth <- 2 # number of powers of 2 to consider in either direction
rvarscales <- (zeta * 2^((-hwidth):hwidth))^2
##
##
nsamples <- 1000
alpha <- sample(rep(alpha0,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
##
sd <- sample(rep(sqrt(rvar0),2),nsamples*nclusters*2,replace=T)
m <- array(rnorm(nsamples*nclusters*2,rmean0,sd),dim=c(nsamples,nclusters,2))
##
shapein <- rshapein0
shapeout <- rshapeout0
scalevar <- sample(rep(rvarscales,2),nsamples*nclusters*2,replace=T)
s <- array(sqrt(nimble::rinvgamma(nsamples*nclusters*2, shape=shapeout,
                                   rate=nimble::rinvgamma(nsamples*nclusters*2, shape=shapein, rate=scalevar))),
            dim=c(nsamples,nclusters,2))
##
xgrid <- seq(xmin, xmax, length.out=256)
ysum <- 0
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
graphics.off()
pdff('priorsamples_real2D_percluster')
par(mfrow=rowcol,mar = c(0,0,0,0))
ycol <- xcol <- NA*numeric(nsamples)
for(i in 1:nsamples){
if(i < nsamples){nxsamples <- 256}else{nxsamples <- 1}
    labs <- sample(1:nclusters, nxsamples, prob=q[i,], replace=T)
    x <- rnorm(nxsamples, m[i,labs,1], s[i,labs,1])
    y <- rnorm(nxsamples, m[i,labs,2], s[i,labs,2])
    xcol[i] <- x[1]
    ycol[i] <- y[1]
if(i < prod(rowcol) || i == nsamples){
    if(i == nsamples){
        x <- xcol[sample(1:length(xcol),256,replace=T)]
        y <- ycol[sample(1:length(ycol),256,replace=T)]
    }
    tplot(x=x,y=y,type='p',pch='.',alpha=0.5,
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5)
    abline(h=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.5,c(7,2)),lty=c(1,2))
    abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.5,c(7,2)),lty=c(1,2))
    }
}
dev.off()


#### continuous variate common among clusters
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
rowcol <- c(10,10)
nsamples <- 1e6
## hyperparameters
nclusters <- 64
alpha0 <- 2^((-3):3)
rmean0 <- 0
zeta <- 1
rvar0 <- (2*zeta)^2
rshapein0 <- 1 # large scales
rshapeout0 <- 1 # small scales
hwidth <- 2 # number of powers of 2 to consider in either direction
rvarscales <- (zeta * 2^((-hwidth):hwidth))^2
##
##
nsamples <- 1000
alpha <- sample(rep(alpha0,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
##
sd <- sample(rep(sqrt(rvar0),2),nsamples*nclusters*2,replace=T)
m <- array(rnorm(nsamples*nclusters*2,rmean0,sd),dim=c(nsamples,2,nclusters))
##
shapein <- rshapein0
shapeout <- rshapeout0
scalevar <- sample(rep(rvarscales,2),nsamples*2,replace=T)
s <- array(sqrt(nimble::rinvgamma(nsamples*nclusters*2, shape=shapeout,
                                   rate=nimble::rinvgamma(nsamples*2, shape=shapein, rate=scalevar))),
            dim=c(nsamples,2,nclusters))
##
xgrid <- seq(xmin, xmax, length.out=256)
ysum <- 0
## tplot(x=xgrid,y=dnorm(txgrid)*jac(xgrid))
graphics.off()
pdff('priorsamples_real2D_commoncluster')
par(mfrow=rowcol,mar = c(0,0,0,0))
ycol <- xcol <- NA*numeric(nsamples)
for(i in 1:nsamples){
if(i < nsamples){nxsamples <- 256}else{nxsamples <- 1}
    labs <- sample(1:nclusters, nxsamples, prob=q[i,], replace=T)
    x <- rnorm(nxsamples, m[i,1,labs], s[i,1,labs])
    y <- rnorm(nxsamples, m[i,2,labs], s[i,2,labs])
    xcol[i] <- x[1]
    ycol[i] <- y[1]
if(i < prod(rowcol) || i == nsamples){
    if(i == nsamples){
        x <- xcol[sample(1:length(xcol),256,replace=T)]
        y <- ycol[sample(1:length(ycol),256,replace=T)]
    }
    tplot(x=x,y=y,type='p',pch='.',alpha=0.5,
          xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
          xticks=NA,yticks=NA,
          mar=c(1,1,1,1)*0.5)
    abline(h=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.5,c(7,2)),lty=c(1,2))
    abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.5,c(7,2)),lty=c(1,2))
    }
}
dev.off()



#### continuous variate per cluster
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
rowcol <- c(16,16)
## hyperparameters
nclusters <- 64
alpha0 <- 2^((-3):3)
rmean0 <- 0
zeta <- 1
hwidthmu <- 1 # number of powers of 2 to consider in either direction
rvar0 <- (zeta * 2^((-hwidthmu):hwidthmu))^2
rshapein0 <- 1 # large scales
rshapeout0 <- 1 # small scales
hwidthsi <- 0 # number of powers of 2 to consider in either direction
rvarscales <- (zeta * 2^((-hwidthsi):hwidthsi))^2
##
##
nsamples <- 2^15
maxxsamples <- 2^10
excl <- 5
alpha <- sample(rep(alpha0,2),nsamples,replace=T)
q <- extraDistr::rdirichlet(n=nsamples,alpha=matrix(alpha/nclusters,nsamples,nclusters))
##
sd <- sample(rep(sqrt(rvar0),2),nsamples*2,replace=T)
m <- array(rnorm(nsamples*2*nclusters, mean=rmean0, sd=sd),dim=c(nsamples,2,nclusters))
##
shapein <- rshapein0
shapeout <- rshapeout0
scalevar <- sample(rep(rvarscales,2),nsamples*2,replace=T)
commonfactor <- 1
s <- array(sqrt(nimble::rinvgamma(nsamples*2*nclusters, shape=shapeout,
                                   rate=nimble::rinvgamma(nsamples*2*commonfactor, shape=shapein, rate=scalevar))),
            dim=c(nsamples,2,nclusters))
##
graphics.off()
pdff(paste0('priorsamples_real2D_com',commonfactor,'_mu',hwidthmu,'_si',hwidthsi,'_Q',excl), apaper=3)
pdf1 <- dev.cur()
pdff(paste0('priorsamples_real1D_com',commonfactor,'_mu',hwidthmu,'_si',hwidthsi,'_Q',excl), apaper=3)
pdf2 <- dev.cur()
for(exclu in c(0,excl)){
    y <- 0
    dev.set(pdf1)
    par(mfrow=rowcol,mar = c(0,0,0,0))
    dev.set(pdf2)
    par(mfrow=rowcol,mar = c(0,0,0,0))
    ycol <- xcol <- NA*numeric(nsamples)
    for(i in 1:nsamples){
        ##set.seed(321)
        if(i < prod(rowcol)){nxsamples <- maxxsamples}else{nxsamples <- 2}
        labs <- sample(1:nclusters, nxsamples, prob=q[i,], replace=T)
        x <- rnorm(nxsamples, mean=m[i,1,labs], sd=s[i,1,labs])
        y <- rnorm(nxsamples, mean=m[i,2,labs], sd=s[i,2,labs])
        xcol[i] <- sample(rep(x,2),1)
        ycol[i] <- sample(rep(y,2),1)
        if(i < prod(rowcol) || i == nsamples){
            if(i == nsamples){
                x <- xcol[sample(1:length(xcol),maxxsamples,replace=T)]
                y <- ycol[sample(1:length(ycol),maxxsamples,replace=T)]
            }
            dev.set(pdf1)
            tplot(x=x,y=y,type='p',pch='.',alpha=0.75,
                  xlim=tquant(x, c(exclu/2,100-exclu/2)/100),
                  ylim=tquant(y, c(exclu/2,100-exclu/2)/100),
                  xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
                  xticks=NA,yticks=NA,
                  mar=c(1,1,1,1)*0.5)
            abline(h=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=par('usr')[1:2],col='black',lwd=0.5)
            abline(h=par('usr')[3:4],col='black',lwd=0.5)
            ##
            ##set.seed(321)
            boundx <- tquant(x, c(exclu/2,100-exclu/2)/100)
            xgrid <- seq(boundx[1],boundx[2], length.out=256)
            if(i < nsamples){
                ygrid <- rowSums(sapply(1:nclusters,function(acluster){
                q[i,acluster] * dnorm(xgrid, m[i,1,acluster], s[i,1,acluster])
                }))
            }else{
                ygrid <- foreach(j=1:nsamples, .combine='+')%dopar%{rowSums(sapply(1:nclusters,function(acluster){
                q[j,acluster] * dnorm(xgrid, m[j,1,acluster], s[j,1,acluster])
                }))}
                }
            ## histx <- thist(x[x >= boundx[1] & x <= boundx[2]], n=128)
            dev.set(pdf2)
            tplot(x=xgrid,y=ygrid,alpha=0,lwd=1,
                  col=(if(i < prod(rowcol)){1}else{3}),
                  ## xlim=tquant(x, c(exclu/2,100-exclu/2)/100),
                  ylim=c(0,NA),
                  xlabels=NA,ylabels=NA, xlab=NA,ylab=NA,
                  xticks=NA,yticks=NA,
                  mar=c(1,1,1,1)*0.5)
            abline(h=c(0),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            abline(v=c(-1,0,1),lwd=0.5,col=alpha2hex2(0.25,c(7,2)),lty=c(1,2))
            ##abline(v=par('usr')[1:2],col='black',lwd=0.5)
            ##abline(h=par('usr')[3:4],col='black',lwd=0.5)
        }
    }
}
dev.off(pdf1)
dev.off(pdf2)













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
nsamplesx <- 2^24
testgr <- c(NULL,
            tquant(rnorm(nsamplesx,
                  mean=rnorm(nsamplesx,mean=rmean0,sd=sqrt(rvar0)),
                  sd=sqrt(
                      extraDistr::rbetapr(nsamplesx,shape1=rshapein0,shape2=rshapeout0,
                                          scale=sample(rvarscales,nsamplesx,replace=T))
                  )
                  ),
                  seq(1/nint,(nint-1)/nint,length.out=nint-1)),
            NULL)
testgr <- (testgr-rev(testgr))/2
approxq <- approxfun(x=seq(1/nint,(nint-1)/nint,length.out=nint-1),y=testgr,yleft=-Inf,yright=+Inf)
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


