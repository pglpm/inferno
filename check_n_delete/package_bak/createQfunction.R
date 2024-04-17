#### Calculate and save transformation function for ordinal variates
createQfunction <- function(nint=8192, nsamples=2^24L, mean=0, sd=1, shapelo=0.5, shapehi=0.5, rate=1, file=paste0('Qfunction',nint), plot=F){
    ##
    seqnint <- (1:(nint-1))/nint
    xsamples <- rnorm(nsamples,
                      mean=rnorm(nsamples,mean=mean,sd=sd),
                      sd=sqrt(nimble::rinvgamma(nsamples, shape=shapelo, rate=nimble::rinvgamma(nsamples, shape=shapehi, rate=rate))
                          ##extraDistr::rbetapr(nsamples, shape1=shapehi, shape2=shapelo, scale=scale)
                      )
                      )
    ##
    thismad <- mad(xsamples,constant=1)
    oquants <- c(NULL,
                 quantile(x=xsamples, probs=seqnint, na.rm=T, type=6),
                 NULL)
    rm(xsamples)
    oquants <- (oquants-rev(oquants))/2
    approxq <- approxfun(x=seqnint, y=oquants, yleft=-Inf, yright=+Inf)
    if(is.character(plot)){
        ##
        nint <- 256
        xgrid <- seq(1/nint,(nint-1)/nint,length.out=nint-1)
        pdff(plot)
        tplot(x=xgrid,y=list(approxq(xgrid),qnorm(xgrid,sd=thismad/qnorm(3/4)),qcauchy(xgrid,scale=thismad)#,qlogis(xgrid,scale=1/qlogis(3/4))
                             ),
              lwd=c(3,2,2,5),lty=c(1,2,4,3), alpha=c(0,rep(0.25,3)),
              ylim=range(approxq(xgrid)), 
              ## xticks=c(0,0.25,0.5,0.75,1),xlabels=c(0,expression(italic(m)/4),expression(italic(m)/2),expression(3*italic(m)/4),expression(italic(m))),
              xlab=expression(italic(x)), ylab=expression(italic(Q)(italic(x))),
              mar=c(NA,5,1,1))
        legend('topleft', c(expression(italic(Q)),'Gauss','Cauchy'), lwd=c(3,2,2,5), lty=c(1,2,4,3), col=c(1,2,3,4), bty='n')
        dev.off()
    }
    if(is.character(file)){
        saveRDS(approxq, paste0(file,'.rds'))
    }
    approxq
}

nsamples <- 2^24L
mean <- 0
sd <- 1
shapelo <- shapehi <- 0.5
rate <- 1
##
means <- rnorm(nsamples,mean=mean,sd=sd)
sds <- sqrt(nimble::rinvgamma(nsamples, shape=shapelo, rate=nimble::rinvgamma(nsamples, shape=shapehi, rate=rate)))

Qf <- readRDS('Qfunction8192.rds')

nint <- 128
seqnint <- (1:(nint-1))/nint
xvals <- Qf(seqnint)
range(xvals)

dsamples <- foreach(x=xvals, .combine=c)%dopar%{mean(dnorm(x, mean=means, sd=sds))}
dsamples <- (dsamples+rev(dsamples))/2
##
##
approxq <- approxfun(x=xvals, y=dsamples, yleft=0, yright=0)
file <- paste0('DQfunction',nint)
saveRDS(approxq, paste0(file,'.rds'))
##

rg <- 2
tplot(list(xvals,xvals/max(xvals)*rg),dsamples,xlim=c(-rg,rg))


dq128 <- readRDS('DQfunction128.rds')
dq256 <- readRDS('DQfunction256.rds')
dq512 <- readRDS('DQfunction512.rds')
dq1024 <- readRDS('DQfunction1024.rds')
dq2048 <- readRDS('DQfunction2048.rds')
rg <- 0.1
xgrid <- seq(-rg,rg,length.out=1024)
tplot(xgrid,list(dq128(xgrid),dq256(xgrid),dq512(xgrid),dq1024(xgrid),dq2048(xgrid)),lty=1,alpha=0.5)

rg <- 1000
xgrid <- seq(0,rg,length.out=1024*1000)
for(fu in list(dq128,dq256,dq512,dq1024,dq2048)){
    print(any(diff(fu(xgrid))>0))
    print(min(fu(xgrid)))
}

rg <- 0.01
xgrid <- 0.05+seq(-rg,rg,length.out=1024)
tplot(xgrid,Qf(xgrid),lty=1,alpha=0.5)





nint2 <- 1024
seqnint2 <- (1:(nint2-1))/nint2
xvals2 <- Qf(seqnint2)
dsamples2 <- diff(seqnint2)/diff(xvals2)
dsamples2 <- (dsamples2+rev(dsamples2))/2
xvals2 <- (xvals2[-1]+xvals2[-length(xvals2)])/2
##
rg <- 1
tplot(xvals2,list(dq128(xvals2),dsamples2),xlim=c(-rg,rg), lty=1, alpha=0.5)


tplot(list(xvals,xvals/max(xvals)*10),xsamples,xlim=c(-10,10))



approxq <- approxfun(x=xvals, y=xsamples, yleft=0, yright=0)
file <- paste0('DQfunction',nint)
saveRDS(approxq, paste0(file,'.rds'))




fu <- dq2048
qq2 <- 0.8
qq2 <- Qf(qq2)
ygrid <- seq(0,1,length.out=1024)
xgrid <- Qf(ygrid)-qq2
tplot(ygrid,dnorm(xgrid,mean=0,sd=0.5)/fu(xgrid+qq2),ylim=c(0,NA))

fu <- dnorm
qq2 <- 0.8
qq2 <- qnorm(qq2)
ygrid <- seq(0,1,length.out=1024)
xgrid <- qnorm(ygrid)-qq2
tplot(ygrid,dnorm(xgrid,mean=0,sd=0.5)/fu(xgrid+qq2),ylim=c(0,NA))



smoothf <- function(x){
    (2*x+c(x[-1],0)+c(0,x[-length(x)]))/4
}

rg <- 1000
xgrid <- seq(-rg,rg,length.out=1024*1000)
fu <- dq2048(xgrid)
print(any(diff(fu)>0))
for(i in 1:100000){
    fu <- smoothf(fu)
    if(!any(diff(fu[xgrid>=0])>0)){
        print(i)
        break
    }
}


rg <- 0.1
tplot(xgrid,list(fu,dq2048(xgrid),dq128(xgrid)),xlim=c(-rg,rg))
