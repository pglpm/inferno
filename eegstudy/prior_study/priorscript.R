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
    ncores <- 4}
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
fb <- function(x){pnorm(x)*(nn+2)/nn-1/nn}
fl <- function(x){exp(x)}
##
nsam <- 400
##
al <- c(1/2,1,2)
sd <- c(1/2,1,2)
sh1 <- c(1/2,1,2)
sh2 <- c(1/2,1,2)
sc <- c(1,2,4,8)
combinations <- t(expand.grid(al,sd,sh1,sh2,sc))
k <- 64
##
foreach(combo=combinations,.combine=c)%dopar%{
    al <- combo[1]
    sd <- combo[2]
    sh1 <- combo[3]
    sh2 <- combo[4]
    sc <- combo[5]
    ##
    q <- extraDistr::rdirichlet(n=nsam,alpha=rep(al/k,k))
    m <- matrix(rnorm(nsam*k,0,sd),k)
    s <- matrix(sqrt(nimble::rinvgamma(nsam*k,shape=sh1,rate=nimble::rinvgamma(nsam*k,shape=sh2,scale=sc))),k)
    graphics.off()
    title <- paste0('k',k,'_a',al,'_s',sd,'_a',sh1,'_b',sh2,'_t',sc)
    ##
    pdff(paste0('priors_',title))
    par(mfrow=c(20,20),mar = c(0,0,0,0))
    for(i in 1:nsam){
        cc <- extraDistr::rcat(1e6,q[i,])
        points <- rnorm(1e6,m[cc,i],s[cc,i])
        ##
        h <- thist(fb(points),n=100)
        tplot(x=h$mids,y=h$density,ylim=c(0,NA),ylabels=NA,xlabels=NA,
              mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,
              ly=1,xlim=c(-2/nn,1+2/nn),lwd=0.5,
              xticks=NA,yticks=NA)
            abline(h=0,lwd=0.5,col=7,lty=1)
        ##
        h <- thist(fl(points),n=100)
        tplot(x=h$mids,y=h$density,ylim=c(0,NA),ylabels=NA,xlabels=NA,
              mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,
              ly=1,xlim=c(0,3),lwd=0.5,
              xticks=NA,yticks=NA)
    abline(h=0,lwd=0.5,col=7,lty=1)
        ##
        h <- thist((points),n=100)
        tplot(x=h$mids,y=h$density,ylim=c(0,NA),ylabels=NA,xlabels=NA,
              mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,
              ly=1,xlim=c(-3,3),lwd=0.5,
              xticks=NA,yticks=NA)
    abline(h=0,lwd=0.5,col=7,lty=1)
        NULL
    }
    dev.off()
}
