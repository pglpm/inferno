library('png')
library('data.table')
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
variateparameters <- fread('_newMC4096ter/variateparameters.csv')
cogvars <- readRDS('cogvars.rds')
mainvar <- readRDS('mainvar.rds')
varNames <- variateparameters$variate

gc()
Ynames <- 'PRM_cor_delayed'
Xnames <- setdiff(varNames, mainvar)
## Xnames <- c( "DMN_T_CPL", "DMN_T_CC")
## Xnames <- setdiff(varNames, c(mainvar,"DMN_T_CPL", "DMN_T_CC"))
nXsamples <- 1024L*8L
nmcsubsamples <- 1024L
set.seed(321)
subsamples <- sample(1:nrow(mcsamples), nmcsubsamples, replace=F)
##
allvariates <- rownames(variateparameters)
locations <- variateparameters[,'location']
scales <- variateparameters[,'scale']
    ##
    rY <- Ynames[variateparameters[Ynames,'type']==0]
    iY <- Ynames[variateparameters[Ynames,'type']==3]
    cY <- Ynames[variateparameters[Ynames,'type']==1]
    bY <- Ynames[variateparameters[Ynames,'type']==2]
    ordYnames <- c(rY,iY,cY,bY)
    nrY <- length(rY)
    niY <- length(iY)
    ncY <- length(cY)
    nbY <- length(bY)
    ##
    rX <- Xnames[variateparameters[Xnames,'type']==0]
    iX <- Xnames[variateparameters[Xnames,'type']==3]
    cX <- Xnames[variateparameters[Xnames,'type']==1]
    bX <- Xnames[variateparameters[Xnames,'type']==2]
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
        totake <- variateparameters[rY,'index']
        imeanRY <- array(sapply(paste0('meanR\\[',totake,','),
                                grep,colnames(mcsamples)),
                         dim=c(nclusters,nrY), dimnames=list(NULL,rY))
        ivarRY <- array(sapply(paste0('varR\\[',totake,','),
                                grep,colnames(mcsamples)),
                         dim=c(nclusters,nrY), dimnames=list(NULL,rY))
        }
    if(niY>0){
        totake <- variateparameters[iY,'index']
        iprobIY <- array(sapply(paste0('probI\\[',totake,','),
                                grep,colnames(mcsamples)),
                         dim=c(nclusters,niY), dimnames=list(NULL,iY))
        isizeIY <- array(sapply(paste0('sizeI\\[',totake,','),
                                grep,colnames(mcsamples)),
                         dim=c(nclusters,niY), dimnames=list(NULL,iY))
        }
    if(ncY>0){
        totake <- variateparameters[cY,'index']
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
        totake <- variateparameters[bY,'index']
        iprobBY <- array(sapply(paste0('probB\\[',totake,','),
                                grep,colnames(mcsamples)),
                         dim=c(nclusters,nbY), dimnames=list(NULL,bY))
    }
    ##
    if(nrX>0){
        totake <- variateparameters[rX,'index']
        imeanRX <- array(t(sapply(paste0('meanR\\[',totake,','),
                                grep,colnames(mcsamples))),
                         dim=c(nrX,nclusters), dimnames=list(rX,NULL))
        ivarRX <- array(t(sapply(paste0('varR\\[',totake,','),
                               grep,colnames(mcsamples))),
                        dim=c(nrX,nclusters), dimnames=list(rX,NULL))
    }
    if(niX>0){
        totake <- variateparameters[iX,'index']
        iprobIX <- array(t(sapply(paste0('probI\\[',totake,','),
                                grep,colnames(mcsamples))),
                         dim=c(niX,nclusters), dimnames=list(iX,NULL))
        isizeIX <- array(t(sapply(paste0('sizeI\\[',totake,','),
                                grep,colnames(mcsamples))),
                         dim=c(niX,nclusters), dimnames=list(iX,NULL))
    }
    if(ncX>0){
        totake <- variateparameters[cX,'index']
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
        totake <- variateparameters[bX,'index']
        iprobBX <- array(t(sapply(paste0('probB\\[',totake,','),
                                grep,colnames(mcsamples))),
                         dim=c(nbX,nclusters), dimnames=list(bX,NULL))
    }
##
    ## Go through the selected MC samples
##graphics.off()
##pdff('testsdreduction')
allmoments <- bind_as_dim(foreach(amcsample=t(mcsamples[subsamples,]))%dorng%{
        amcsample <- c(amcsample)
        aq <- amcsample[Qi]
        ##
        if(nrY>0){
        meanRY <- replace(imeanRY,1:length(imeanRY),amcsample[imeanRY])
        varRY <- replace(ivarRY,1:length(ivarRY),amcsample[ivarRY])
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
        ## Draw X samples
        drawclusters <- sample(x=sclusters, size=nXsamples, prob=aq, replace=TRUE)
        if(nrX>0){
             drawrX <- array(rnorm(n=nXsamples*nrX, mean=meanRX[,drawclusters], sd=sdRX[,drawclusters]), dim=c(nrX,nXsamples), dimnames=list(rX,NULL))
         }else{drawrX <- NULL}
        ##        
        drawiX <- NULL
        ##
        if(ncX>0){
             drawcX <- array(extraDistr::rcat(n=nXsamples*ncX, prob=matrix(probCX[,drawclusters,],nrow=ncX*nXsamples)), dim=c(ncX,nXsamples), dimnames=list(cX,NULL))
         }else{drawcX <- NULL}
        ##
        if(nbX>0){
             drawbX <- array(extraDistr::rbern(n=nbX*nXsamples, prob=probBX[,drawclusters]), dim=c(nbX,nXsamples), dimnames=list(bX,NULL))
         }else{drawbX <- NULL}
        ##
        ## ## Calculate SD of Y conditional on the X samples
        ##
        ## First the conditional mixture weights
        pX <- t(exp(#cols=clusters
                log(aq) + t(#rows=clusters
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
                }, numeric(nXsamples)))
            ))
        pX <- cbind(aq,t(pX/rowSums(pX))) # rows=clusters, cols=Xsamples
        ##
        ## Now the means and SDs of Y
        if(nrY>0){
            rM <- t(sapply(rY, function(acov){
                out <- meanRY[,acov]
                rbind(out2 <- colSums(pX*out),
                      colSums(pX * (varRY[,acov] + out*out)) -out2*out2)
            })) # Y, moments, samples
        }else{rM <- NULL}
        ##
        if(niY>0){
            iM <- t(sapply(iY, function(acov){
                out <- probIY[,acov] * sizeIY[,acov]
                rbind(out2 <- colSums(pX*out),
                      colSums(pX*out*
                              (out - probIY[,acov]+1)) -out2*out2)
            })) # Y, moments, samples
        }else{iM <- NULL}
        if(exists('out2')){rm(out2)}
        ##
        if(ncY>0){
            cM <- t(sapply(cY, function(acov){
                out <- probCY[,,acov]*scategories
                rbind(colSums(pX*colSums(out)),
                      colSums(pX*colSums(out*scategories)))
            })) # Y, moments, samples
        }else{cM <- NULL}
        ##
        if(nbY>0){
            bM <- t(sapply(bY, function(acov){
                rbind(out <- colSums(pX*probBY[,acov]),
                      out*(1-out))
            })) # Y, moments, samples
        }else{bM <- NULL}
        if(exists('out')){rm(out)}
        ##
        array(rbind(rM,iM,cM,bM),
              dim=c(length(ordYnames), 2, nXsamples+1),
              dimnames=list(ordYnames,c('mean','var'),NULL)
              )[order(match(ordYnames,Ynames)),,,drop=F] *
            c(scales[Ynames],scales[Ynames]*scales[Ynames]) +
            c(locations[Ynames],locations[Ynames]*0L)
}, -1)
## END sampling

graphics.off()
## pdff(paste0('sdreduction-',Ynames,'--using-DMN_T'))
## pdff(paste0('sdreduction-',Ynames,'--using-allgraph'))
pdff(paste0('sdreduction-',Ynames,'--using-allgraph_minus_DMN_T'))
allredu <- numeric(nmcsubsamples)
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
