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
mcsamples <- readRDS('_newMC4096ter/_jointmcsamples-_newMC4096ter-4096.rds')
varinfo <- data.matrix(read.csv('_newMC4096ter/variateparameters.csv',row.names=1,header=T))
mainvar <- readRDS('cogvars.rds')
source(paste0('functions_mcmc.R'))
varNames <- rownames(varinfo)
othervars <- setdiff(varNames, mainvar)
alldata <- fread('data_ep.csv', sep=',')
## mcsamples <- mcsamples[1:4,]
## testv <- c(sapply(0:2,function(i)rownames(varinfo)[varinfo[,'type']==i][1:2]))

gc()
## Xnames <- sapply(0:2,function(i)rownames(varinfo[varinfo[,'type']==i,])[1])
Xnames <- setdiff(varNames, mainvar)
## Xnames <- c( "DMN_T_CPL", "DMN_T_CC")
## Xnames <- setdiff(varNames, c(mainvar,"DMN_T_CPL", "DMN_T_CC"))
Ynames <- 'PRM_cor_delayed'
nmcsubsamples <- nrow(mcsamples)
set.seed(321)
subsamples <- sample(1:nrow(mcsamples), nmcsubsamples, replace=F)


Xvalues <- samplesXmc(pointspermcsample=1,mcsamples=mcsamples,Xnames=Xnames,variateparameters=varinfo,seed=123,inorder=F)[,1,]
nXsamples <- nrow(Xvalues)


##
locations <- varinfo[,'location']
scales <- varinfo[,'scale']
names(scales) <- names(locations) <- rownames(varinfo)
##
rY <- Ynames[varinfo[Ynames,'type']==0]
iY <- Ynames[varinfo[Ynames,'type']==3]
cY <- Ynames[varinfo[Ynames,'type']==1]
bY <- Ynames[varinfo[Ynames,'type']==2]
ordYnames <- c(rY,iY,cY,bY)
nrY <- length(rY)
niY <- length(iY)
ncY <- length(cY)
nbY <- length(bY)
##
rX <- Xnames[varinfo[Xnames,'type']==0]
iX <- Xnames[varinfo[Xnames,'type']==3]
cX <- Xnames[varinfo[Xnames,'type']==1]
bX <- Xnames[varinfo[Xnames,'type']==2]
ordXnames <- c(rX,iX,cX,bX)
nrX <- length(rX)
niX <- length(iX)
ncX <- length(cX)
nbX <- length(bX)
##
Qi <- grep('q',colnames(mcsamples))
nclusters <- length(Qi)
sclusters <- seq_len(nclusters)
if(nrY>0){
    totake <- varinfo[rY,'index']
    imeanRY <- array(sapply(paste0('meanR\\[',totake,','),
                            grep,colnames(mcsamples)),
                     dim=c(nclusters,nrY), dimnames=list(NULL,rY))
    ivarRY <- array(sapply(paste0('varR\\[',totake,','),
                           grep,colnames(mcsamples)),
                    dim=c(nclusters,nrY), dimnames=list(NULL,rY))
}
if(niY>0){
    totake <- varinfo[iY,'index']
    iprobIY <- array(sapply(paste0('probI\\[',totake,','),
                            grep,colnames(mcsamples)),
                     dim=c(nclusters,niY), dimnames=list(NULL,iY))
    isizeIY <- array(sapply(paste0('sizeI\\[',totake,','),
                            grep,colnames(mcsamples)),
                     dim=c(nclusters,niY), dimnames=list(NULL,iY))
}
if(ncY>0){
    totake <- varinfo[cY,'index']
    iprobCY <- sapply(paste0('probC\\[',totake,','),
                      grep,colnames(mcsamples))
    ncategories <- length(iprobCY)/ncY/nclusters
    scategories <- seq_len(ncategories)
    iprobCY <- aperm(array(iprobCY,
                           dim=c(nclusters,ncategories,ncY),
                           dimnames=list(NULL,NULL,nY)),
                     c(2,1,3))
}
if(nbY>0){
    totake <- varinfo[bY,'index']
    iprobBY <- array(sapply(paste0('probB\\[',totake,','),
                            grep,colnames(mcsamples)),
                     dim=c(nclusters,nbY), dimnames=list(NULL,bY))
}
##
if(nrX>0){
    totake <- varinfo[rX,'index']
    imeanRX <- array(t(sapply(paste0('meanR\\[',totake,','),
                              grep,colnames(mcsamples))),
                     dim=c(nrX,nclusters), dimnames=list(rX,NULL))
    ivarRX <- array(t(sapply(paste0('varR\\[',totake,','),
                             grep,colnames(mcsamples))),
                    dim=c(nrX,nclusters), dimnames=list(rX,NULL))
}
if(niX>0){
    totake <- varinfo[iX,'index']
    iprobIX <- array(t(sapply(paste0('probI\\[',totake,','),
                              grep,colnames(mcsamples))),
                     dim=c(niX,nclusters), dimnames=list(iX,NULL))
    isizeIX <- array(t(sapply(paste0('sizeI\\[',totake,','),
                              grep,colnames(mcsamples))),
                     dim=c(niX,nclusters), dimnames=list(iX,NULL))
}
if(ncX>0){
    totake <- varinfo[cX,'index']
    iprobCX <- sapply(paste0('probC\\[',totake,','),
                      grep,colnames(mcsamples))
    ncategories <- length(iprobCX)/ncX/nclusters
    scategories <- seq_len(ncategories)
    iprobCX <- aperm(array(iprobCX,
                           dim=c(nclusters,ncategories,ncX),
                           dimnames=list(NULL,NULL,cX)),
                     c(3,1,2))
}
if(nbX>0){
    totake <- varinfo[bX,'index']
    iprobBX <- array(t(sapply(paste0('probB\\[',totake,','),
                              grep,colnames(mcsamples))),
                     dim=c(nbX,nclusters), dimnames=list(bX,NULL))
}
##

## Go through the selected MC samples
##graphics.off()
##pdff('testsdreduction')
##
## Draw X samples
if(nrX>0){
    drawrX <- (t(Xvalues[,rX,drop=F])-locations[rX])/scales[rX]
}else{drawrX <- NULL}
##        
drawiX <- NULL
##
if(ncX>0){
    drawcX <- (t(Xvalues[,cX,drop=F])-locations[cX])/scales[cX]
}else{drawcX <- NULL}
##
if(nbX>0){
    drawbX <- (t(Xvalues[,bX,drop=F])-locations[bX])/scales[bX]
}else{drawbX <- NULL}

## ## Calculate SD of Y conditional on the X samples
##
## First the conditional mixture weights
Ymin <- (0-locations[Ynames])/scales[Ynames]
Ymax <- (100-locations[Ynames])/scales[Ynames]
fn <- function(...){Map('+',...)}
allmoments <- foreach(amcsample=t(mcsamples), .combine=fn, .inorder=F)%dopar%{
    amcsample <- c(amcsample)
    aq <- amcsample[Qi]
    ##
    if(nrY>0){
        mu0 <- replace(imeanRY,1:length(imeanRY),amcsample[imeanRY])
        si0 <- sqrt(replace(ivarRY,1:length(ivarRY),amcsample[ivarRY]))
        a0 <- (Ymin-mu0)/si0
        b0 <- (Ymax-mu0)/si0
        da0 <- dnorm(a0)
        db0 <- dnorm(b0)
        z0 <- pnorm(b0) - pnorm(a0)
        ##
        meanRY <- mu0 - si0*(db0-da0)/z0
        varRY <- si0*si0*(1 - (b0*db0-a0*da0)/z0 - ((db0-da0)/z0)^2)
    }
    if(ncY>0){
        probCY <- replace(iprobCY,1:length(iprobCY),amcsample[iprobCY])
    }
    if(nbY>0){
        probBY <- replace(iprobBY,1:length(iprobBY),amcsample[iprobBY])
    }
    ##
    if(nrX>0){
        meanRX <- replace(imeanRX,1:length(imeanRX),amcsample[imeanRX])
        sdRX <- sqrt(replace(ivarRX,1:length(ivarRX),amcsample[ivarRX]))
    }
    if(ncX>0){        
        probCX <- replace(iprobCX,1:length(iprobCX),amcsample[iprobCX])
    }
    if(nbX>0){
        probBX <- replace(iprobBX,1:length(iprobBX),amcsample[iprobBX])
    }
    ##
    pX <- t(log(aq) + t(#rows=clusters
                      vapply(sclusters, function(acluster){#cols=clusters
                          ## real covariates (rows)
                          (if(nrX>0){
                               colSums(dnorm(x=drawrX, mean=meanRX[,acluster], sd=sdRX[,acluster], log=T), na.rm=T)
                           }else{0}) +
                              ## integer covariates (rows)
                              (if(niX>0){
                                   colSums(dbinom(x=drawiX, prob=probIX[,acluster], size=sizeIX[,acluster], log=T), na.rm=T)
                               }else{0}) +
                              ## category variates (cols)
                              (if(ncX>0){
                                   colSums(matrix(extraDistr::dcat(x=drawcX, prob=probCX[,acluster,], log=T), nrow=ncX), na.rm=T)
                               }else{0}) +
                              ## binary covariates (rows)
                              (if(nbX>0){
                                   colSums(matrix(extraDistr::dbern(x=drawbX, prob=probBX[,acluster], log=T), nrow=nbX), na.rm=T)
                               }else{0})
                      }, numeric(nXsamples))
                      ))
    ##
    pX <- exp(pX - apply(pX,1,max,na.rm=T)) # max(lpX,na.rm=T)
    ## pX <- cbind(aq,t(pX)) # rows=clusters, cols=Xsamples
    ## cpX <- colSums(pX)
    ## cpX[cpX==0] <- 1
    pX <- cbind(aq,t(pX/rowSums(pX))) # rows=clusters, cols=Xsamples
    pX[is.na(pX)] <- 0
    ##
    ## Now the means and SDs of Y
    if(nrY>0){
        rM <- t(sapply(rY, function(acov){
            out <- meanRY[,acov]
            out2 <- colSums(pX*out)
            rbind(out2,
                  colSums(pX * varRY[,acov]) +
                  colSums(pX * outer(out,out2,'-')^2)
                  )
            ## out3 <- colSums(pX*out)/cpX
            ## rbind(out3 * scales[acov]+locations[acov],
            ##       colSums(pX * (varRY[,acov] + out*out))/cpX * scales[acov]*scales[acov] +2*scales[acov]*locations[acov]*out3 + locations[acov]*locations[acov]
            ##       )
        })) # Y, moments, samples
    }else{rM <- NULL}
    ##
    if(niY>0){
        iM <- t(sapply(iY, function(acov){
            out <- probIY[,acov] * sizeIY[,acov]
            rbind(colSums(pX*out),
                  colSums(pX*out*
                          (out - probIY[,acov]+1))
                  )
        })) # Y, moments, samples
    }else{iM <- NULL}
    ##
    if(ncY>0){
        cM <- t(sapply(cY, function(acov){
            out <- probCY[,,acov]*scategories
            rbind(colSums(pX*colSums(out)),
                  colSums(pX*colSums(out*scategories))
                  )
        })) # Y, moments, samples
    }else{cM <- NULL}
    ##
    if(nbY>0){
        bM <- t(sapply(bY, function(acov){
            out <- colSums(pX*probBY[,acov])
            rbind(out,
                  out
                  )
        })) # Y, moments, samples
    }else{bM <- NULL}
##    if(exists('out')){rm(out)}
    ##
    out <- array(rbind(rM,iM,cM,bM),
                  dim=c(length(ordYnames), 2, nXsamples+1),
                  dimnames=list(ordYnames,c('mean','variance'),NULL)
                  )[order(match(ordYnames,Ynames)),,,drop=F]
    ##    test[is.na(test)] <- 0
    if(!any(is.na(out))){list(out,1)}else{list(0,0)}
}
## END sampling



Ysds <- sqrt(allmoments[[1]][1,2,]/allmoments[[2]])*scales[Ynames]

uncYsds <- Ysds[1]
condYsds <- Ysds[-1]
bestX <- order(condYsds)[1:3]
mehX <- order(abs(condYsds-uncYsds))[1:3]
showX <- c(bestX,mehX)

histo <- thist(condYsds,plot=F)
graphics.off()
pdff(paste0('sdreductionHisto-',Ynames,'--using-allgraph'))
tplot(x=histo$breaks,y=histo$density,
      xlab='standard deviation for prediction',
      ylab='probability',
      main=paste0('median improvement in SD: ',signif((median(condYsds)/uncYsds-1)*100,3),'%')
      )
abline(v=uncYsds,col=darkgrey,lwd=3,lty=2)
text(uncYsds,mean(par('usr')[3:4]),paste0('without knowledge\nof graph variates'),col=darkgrey,pos=4)
abline(v=median(condYsds),col=palette()[4],lwd=3,lty=1)
xti <- par('xaxp')
xti <- seq(xti[1],xti[2],length.out=xti[3]+1)
axis(side=3,at=xti,labels=paste0(signif((xti/uncYsds-1)*100,2),'%'),tick=F)
dev.off()

## rgY <- extendrange(range(alldata[[Ynames]]),by=1/4)
rgYplot <- range(alldata[[Ynames]])
rgY <- range(varinfo[Ynames,c('thmax','thmin')])
Ygrid <- matrix(seq(rgY[1],rgY[2],length.out=512),ncol=1,dimnames=list(NULL,Ynames))
##

uncdistrY <- samplesFmc(Y=Ygrid,X=NULL,mcsamples=mcsamples,variateparameters=varinfo,inorder=F)

distrYbest <- samplesFmc(Y=Ygrid,X=Xvalues[bestX[1],,drop=F],mcsamples=mcsamples,variateparameters=varinfo,inorder=F)

distrYmeh <- samplesFmc(Y=Ygrid,X=Xvalues[mehX[2],,drop=F],mcsamples=mcsamples,variateparameters=varinfo,inorder=F)

apmean <- sum(rowMeans(uncdistrY,na.rm=T)*c(Ygrid))/sum(rowMeans(uncdistrY,na.rm=T))
apsd <- sqrt(sum(rowMeans(uncdistrY,na.rm=T)*c(Ygrid-apmean)^2)/sum(rowMeans(uncdistrY,na.rm=T)))

besmean <- sum(rowMeans(distrYbest,na.rm=T)*c(Ygrid))/sum(rowMeans(distrYbest,na.rm=T))
bessd <- sqrt(sum(rowMeans(distrYbest,na.rm=T)*c(Ygrid-apmean)^2)/sum(rowMeans(distrYbest,na.rm=T)))




graphics.off()
## pdff(paste0('sdreduction2-',Ynames,'--using-DMN_T'))
pdff(paste0('sdreductionExamples-',Ynames,'--using-allgraph'))
## pdff(paste0('sdreduction2-',Ynames,'--using-allgraph_minus_DMN_T'))
tplot(x=Ygrid,y=list(rowMeans(uncdistrY,na.rm=T),
                     rowMeans(distrYbest,na.rm=T),
                     rowMeans(distrYmeh,na.rm=T)),
      lty=c(1,2,4), lwd=c(6,3,3), col=c(7,1,2),
      xlab=Ynames, ylab='probability',ylim=c(0,NA)
      )
legend(x='topleft',legend=c('no predictors','good graph-predictor value','bad graph-predictor value'),
       lty=c(1,2,4),lwd=c(6,3,3),col=c(7,1,2),
       bty='n')
redhist <- thist(condYsds[!is.na(condYsds)],n=-abs(uncYsds-min(condYsds,na.rm=T))/20)
rgh <- range(c(min(condYsds,na.rm=T), quantile(condYsds,0.95,na.rm=T), uncYsds))
tplot(x=redhist$breaks,y=redhist$density,xlim=rgh,
      xlab=paste0('standard deviation of probability for ',Ynames),ylab='probability',ylim=c(0,NA))
abline(v=uncYsds,lwd=4,lty=1,col=grey)
abline(v=(meds <- median(condYsds,na.rm=T)),lwd=4,lty=2,col=bluepurple)
text(uncYsds,max(redhist$density)/2,paste0('without\nknowledge\nof graph data:\n',signif(uncYsds,3)),col=darkgrey,pos=2)
legend(x='top',legend=paste0('median: ',signif(meds,3)),text.col=bluepurple,bty='n')
dev.off()





stop()











redhisto <- thist(allredu,n=-5)
tplot(x=redhisto$breaks,y=redhisto$density,xlab='median reduction',
      ylab='probability',
      xtransf=function(x){paste0(x,'%')})
abline(v=0,lty=2,lwd=4,col=redpurple)
abline(v=median(allredu),lty=3,lwd=4,col=green)
dev.off()



for(i in 1:nmcsubsamples){
    origsdy <- sqrt(allmoments[1,2,1,i])
    newsdy <- sqrt(allmoments[1,2,-1,i])
##    sdhisto <- thist(newsdy,n=-abs(origsdy-quant(newsdy,1/8))/4)
##    sdhisto <- thist(newsdy,n=-abs(diff(quant(newsdy,c(0.1,0.9))))/10)
    sdhisto <- thist(newsdy,n=-0.5)
    rg <- range(c(quant(newsdy,c(0.1,0.9)),origsdy))
    redu <- round(median(newsdy/origsdy-1)*100)
    allredu[i] <- redu
    if(i<17){
    tplot(x=sdhisto$mids,y=sdhisto$density,
          xlim=rg,ylim=c(0,NA),
          xlab=paste0('SD of ',Ynames),
          ylab='probability',
          lty=1,lwd=2,alpha=0,col=7,
          main=paste0('median change: ',redu,'%'))
    abline(v=origsdy,lwd=3,lty=2,col=redpurple)
    text(origsdy,max(sdhisto$density)/2,'without\nknowledge\nof graph data',col=redpurple,pos=2)
    }
}
redhisto <- thist(allredu,n=-5)
tplot(x=redhisto$breaks,y=redhisto$density,xlab='median reduction',
      ylab='probability',
      xtransf=function(x){paste0(x,'%')})
abline(v=0,lty=2,lwd=4,col=redpurple)
abline(v=median(allredu),lty=3,lwd=4,col=green)
dev.off()





f <- function(x){pnorm(x)}
nn <- 64
f <- function(x){(x<0)*((1)/(nn*exp(nn*log(abs(1-x))))-1/nn)+
                     (x>1)*((-1)/(nn*exp(nn*log(abs(x))))+1/nn+1)+
                     (x>=0 & x<=1)*x }
##
tplot(x=(xgrid <- seq(-0.5,1.5,length.out=512)),y=f(xgrid))

f <- function(x){pcauchy(x)}

qnorm(seq(0,1,length.out=33))

nn <- 16
qnorm((seq(0,1,length.out=33)+1/nn)*nn/(nn+2))

set.seed(123)
nn <- 32
f <- function(x){pnorm(x)*(nn+2)/nn-1/nn}
##
nsam <- 400
sd <- 1/2
al <- 1
sh1 <- 1
sh2 <- 1
sc <- 8
k <- 64
q <- extraDistr::rdirichlet(n=nsam,alpha=rep(al/k,k))
m <- matrix(rnorm(nsam*k,0,sd),k)
s <- matrix(sqrt(nimble::rinvgamma(nsam*k,shape=sh1,rate=nimble::rinvgamma(nsam*k,shape=sh2,scale=sc))),k)
graphics.off()
pdff('test')
par(mfrow=c(20,20),mar = c(0,0,0,0))
for(i in 1:nsam){
    cc <- extraDistr::rcat(1e6,q[i,])
    points <- rnorm(1e6,m[cc,i],s[cc,i])
    h <- thist(points,n=100)
    tplot(x=h$mids,y=h$density,ylim=c(0,NA),ylabels=NA,xlabels=NA,
          mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,
          ly=1,xlim=c(-2/nn,1+2/nn),lwd=0.5,
          xticks=NA,yticks=NA)
#    tplot(x=f(seq(-2,2,length.out=16)),y=rep(0,16),type='p',pch=16,cex=0.5,col=2,add=T)
}
dev.off()
##
file.copy('test.pdf',paste0('prior_k',k,'_a',al,'_sd',sd,'_sa',sh1,'_sb',sh2,'_sc',sc,'.pdf'))
    
set.seed(123)
nn <- 32
f1 <- function(x){pnorm(x)*(nn+2)/nn-1/nn}
f2 <- function(x){exp(x)}
##
nsam <- 40
sd <- 1/2
al <- 1
sh1 <- 1
sh2 <- 1
sc <- 8
k <- 64
q <- extraDistr::rdirichlet(n=nsam,alpha=rep(al/k,k))
m <- matrix(rnorm(nsam*k,0,sd),k)
s <- matrix(sqrt(nimble::rinvgamma(nsam*k,shape=sh1,rate=nimble::rinvgamma(nsam*k,shape=sh2,scale=sc))),k)
graphics.off()
pdff('test')
par(mfrow=c(20,20),mar = c(0,0,0,0))
for(i in 1:nsam){
    cc <- extraDistr::rcat(1e6,q[i,])
    points <- rnorm(1e6,m[cc,i],s[cc,i])
    h <- thist(points,n=100)
    tplot(x=h$mids,y=h$density,ylim=c(0,NA),ylabels=NA,xlabels=NA,
          mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,
          ly=1,xlim=c(-2/nn,1+2/nn),lwd=0.5,
          xticks=NA,yticks=NA)
    abline(h=0,lwd=0.5,col=7,lty=1)
#    tplot(x=f(seq(-2,2,length.out=16)),y=rep(0,16),type='p',pch=16,cex=0.5,col=2,add=T)
}
dev.off()


##
file.copy('test.pdf',paste0('prior_k',k,'_a',al,'_sd',sd,'_sa',sh1,'_sb',sh2,'_sc',sc,'.pdf'))
    



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
combinations <- expand.grid(al,sd,sh1,sh2,sc)
k <- 64
q <- extraDistr::rdirichlet(n=nsam,alpha=rep(al/k,k))
m <- matrix(rnorm(nsam*k,0,sd),k)
s <- matrix(sqrt(nimble::rinvgamma(nsam*k,shape=sh1,rate=nimble::rinvgamma(nsam*k,shape=sh2,scale=sc))),k)
graphics.off()
title <- paste0('k',k,'_a',al,'_s',sd,'_a',sh1,'_b',sh2,'_t',sc)
##
pdff(paste0('bounded_',title))
pdfb <- dev.cur()
par(mfrow=c(20,20),mar = c(0,0,0,0))
pdff(paste0('log_',title))
pdfl <- dev.cur()
par(mfrow=c(20,20),mar = c(0,0,0,0))
pdff(paste0('id_',title))
pdfi <- dev.cur()
par(mfrow=c(20,20),mar = c(0,0,0,0))
for(i in data.table(t(combinations))){
    cc <- extraDistr::rcat(1e6,q[i,])
    points <- rnorm(1e6,m[cc,i],s[cc,i])
    ##
    dev.set(pdfb)
    h <- thist(fb(points),n=100)
    tplot(x=h$mids,y=h$density,ylim=c(0,NA),ylabels=NA,xlabels=NA,
          mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,
          ly=1,xlim=c(-2/nn,1+2/nn),lwd=0.5,
          xticks=NA,yticks=NA)
    ##
    dev.set(pdfl)
    h <- thist(fl(points),n=100)
    tplot(x=h$mids,y=h$density,ylim=c(0,NA),ylabels=NA,xlabels=NA,
          mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,
          ly=1,xlim=c(0,3),lwd=0.5,
          xticks=NA,yticks=NA)
    ##
    dev.set(pdfi)
    h <- thist((points),n=100)
    tplot(x=h$mids,y=h$density,ylim=c(0,NA),ylabels=NA,xlabels=NA,
          mar=c(1,1,1,1)*0.5,xlab=NA,ylab=NA,
          ly=1,xlim=c(-3,3),lwd=0.5,
          xticks=NA,yticks=NA)
}
dev.off(pdfb)
dev.off(pdfl)
dev.off(pdfi)




dev.off()
##
file.copy('test.pdf',paste0('prior_k',k,'_a',al,'_sd',sd,'_sa',sh1,'_sb',sh2,'_sc',sc,'.pdf'))



graphics.off()
foreach(i=1:8, .combine=c)%dopar%{
    graphics.off()
    pdff(paste0('testcom',i))
    par(mfrow=c(2,2),mar = c(0,0,0,0))
    for(j in 1:4){
        tplot(i+(0:1),j+(0:1))
        }
dev.off()
}






dx <- 1e-4
ygrid <- seq(-16,0,length.out=256)
tplot(x=ygrid,
      y=(pnorm(ygrid,log.p=T) - dnorm(ygrid,log=T) + dlogis(ygrid,log=T))/log(10))
abline(v=-12*log(2))


tplot(x=xgrid,
      y=pnorm(xgrid+dx)/(2*pnorm(xgrid)))
tplot(x=xgrid,
      y=0.5+0.5*(dnorm(xgrid)*dx/pnorm(xgrid)), add=T,col=2,lty=2)


aa <- -(1+exp(-8))/(exp(8)-exp(-8))
aa2 <- -(1+exp(8))/(exp(16)-1)
bb <- (1+exp(8))/(exp(8)-exp(-8))
dd <- exp(log(MM-mm)+log(2^-expo+2+2^expo)-log(2^-expo+2^expo))

expo <- 12
expo2 <- expo/2
mm <- 0
MM <- 100
aa <- (mm*2^expo2-MM*2^-expo2)/(2^expo2-2^-expo2)
bb <- (MM*2^expo2-mm*2^-expo2)/(2^expo2-2^-expo2)
dd <- exp(log(MM-mm)+log1p(2^-expo)-log1p(-2^-expo))
##
aa2 <- (mm*(1+2^expo)-MM*(1+2^-expo))/(2^expo-2^-expo)
bb2 <- (MM*(1+2^expo)-mm*(1+2^-expo))/(2^expo-2^-expo)
dd2 <- exp(log(MM-mm)+log(2^-expo+2+2^expo)-log(2^expo-2^-expo))
##
xgrid <- seq(mm,MM,length.out=256)
tplot(x=xgrid, y=qlogis((xgrid-aa)/dd))
##
ygrid <- seq(-8,8,length.out=256)
tplot(y=ygrid, x=plogis(ygrid)*dd+aa, add=T,col=2,lty=2)

xgrid2 <- plogis(qlogis((xgrid-aa)/dd))*dd+aa


tplot(x=xgrid,
      y=qlogis((xgrid-aa)/(bb-aa)),add=T,col=2,lty=2)
## tplot(x=xgrid,
##       y=qlogis(xgrid, location=(aa+bb)/2, scale=dd),add=T,col=2,lty=2)

tplot(x=xgrid,
      y=qlogis((xgrid-aa2)/(bb2-aa2)),add=T,col=2,lty=2)

tplot(x=xgrid,
      y=qlogis((xgrid-aa)/dd),add=T,col=3,lty=3)



tplot(x=xgrid,
      y=log((xgrid-aa)/(bb-xgrid)))
tplot(x=xgrid,
      y=qlogis((xgrid-aa)/(bb-aa)),add=T,col=2,lty=2)
tplot(x=xgrid,
      y=qlogis((xgrid-aa)/dd),add=T,col=3,lty=3)

mm <- 0
MM <- 1
me <- 1
q1 <- 0.1
q2 <- 1
##
fn <- function(a){
    si <- log(q2-a)-log(q1-a)
    ya <- (log(mm-a) - log(me-a))/si
    pnorm(ya)/dnorm(ya)*si*(mm-a)
}
##
lagrid <- seq(0,-14,length.out=512)
agrid <- -2^lagrid
tplot(x=lagrid,
      y=log10(fn(agrid)))




mm <- 0
MM <- 100
dd <- MM-mm
me <- 50
q1 <- 40
q2 <- 60
tra <- function(y){qlogis((y-mm)/dd*(1-2*a)+a)}
##
fn <- function(a){
    si <- tra(q2)-tra(q1)
    mu <- tra(me)
    pnorm(mm)/dnorm(mm)*dlogis(mmsi*(mm-a)
}
##
lagrid <- seq(0,-14,length.out=512)
agrid <- -2^lagrid
tplot(x=lagrid,
      y=log10(fn(agrid)))
