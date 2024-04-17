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
