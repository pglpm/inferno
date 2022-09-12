## Author: PGL  Porta Mana
## Created: 2022-09-09T18:23:21+0200
## Last-Updated: 2022-09-11T15:42:41+0200
################
## Test script for VB's data analysis
################

## load customized plot functions
source('~/work/pglpm_plotfunctions.R')
##
## Read MCMC seed from command line
mcmcseed = as.integer(commandArgs(trailingOnly=TRUE))[1]
if(is.na(mcmcseed) | (!is.na(mcmcseed) & mcmcseed <=0)){mcmcseed <- 941}
print(paste0('MCMC seed = ',mcmcseed))
##
#### Packages and setup ####
library('data.table')
library('png')
library('foreach')
library('doFuture')
library('doRNG')
registerDoFuture()
print('availableCores:')
print(availableCores())
print('availableCores-multicore:')
print(availableCores('multicore'))
ncores <- 4#availableCores()-1
print(paste0('using ',ncores,' cores'))
if(ncores>1){
    if(.Platform$OS.type=='unix'){
        plan(multicore, workers=ncores)
    }else{
        plan(multisession, workers=ncores)
    }
}else{
    plan(sequential)
}
#
source('functions_mcmc.R') # load functions for post-MCMC calculations

saveinfofile <- 'variateinfo_1.csv'
variateinfo <- fread(saveinfofile, sep=',')
covNames <- variateinfo$variate
covTypes <- variateinfo$type
covMins <- variateinfo$min
covMaxs <- variateinfo$max
names(covTypes) <- names(covMins) <- names(covMaxs) <- covNames
##
maincov <- 'Group'
##
realCovs <- covNames[covTypes=='real']
integerCovs <- covNames[covTypes=='integer']
binaryCovs <- covNames[covTypes=='binary']
covNames <- c(realCovs, integerCovs, binaryCovs)
covTypes <- covTypes[covNames]
covMins <- covMins[covNames]
covMaxs <- covMaxs[covNames]
nrcovs <- length(realCovs)
nicovs <- length(integerCovs)
nbcovs <- length(binaryCovs)
ncovs <- length(covNames)


## Combine parallel MCMC samples
dirname <- '_mcmcsaga-V13-D848-K64-I1024'
##
npar <- 8
ntotal <- 1024*4
parmList <- mcsamples2parmlist(
foreach(i=1:npar, .combine=rbind)%dopar%{
        temp <- readRDS(paste0(dirname, '/_mcsamples-R_mcmctest_',i,'_1-V13-D848-K64-I1024.rds'))
        nskip <- ceiling(nrow(temp)/(ntotal/npar))
        if(nskip<1){print('WARNING! not enough points')}
        temp[seq(nrow(temp),1,by=-nskip),]
        ## if(any(is.na(nrow(temp)+1-rev(seq(1,nrow(temp),by=nskip)[1:(ntotal/npar)])))){print('WARNING! not enough points')}
        ## temp[nrow(temp)+1-rev(seq(1,nrow(temp),by=nskip)[1:(ntotal/npar)]),]
}
, realCovs, integerCovs, binaryCovs)


outfile <- paste0(dirname,'/','results_saga.txt')
nclusters <- ncol(parmList$q)
nFsamples <- nrow(parmList$q)
##
otherCovs <- setdiff(covNames, maincov)
maincovnames <- c('0', '1')


datafile <- 'Cortical_myelination_faux.csv'
alldata <- fread(datafile, sep=',')
alldata <- alldata[,..covNames]
## Grids for plots
grids <- foreach(acov=covNames)%do%{
    rg <- range(alldata[[acov]], na.rm=T)
    if(acov %in% realCovs){
        rg <- rg+c(-1,1)*IQR(alldata[[acov]],type=8,na.rm=T)/2
        Xgrid <- seq(rg[1], rg[2], length.out=256)
    }else if(acov %in% integerCovs){
            rg <- round(c((covMins[acov]+7*rg[1])/8, (covMaxs[acov]+7*rg[2])/8))
            Xgrid <- rg[1]:rg[2]
    }else{Xgrid <- 0:1}
    matrix(Xgrid, ncol=1, dimnames=list(NULL,acov))
}
names(grids) <- covNames
##
xcond <- matrix(0:1,ncol=2,dimnames=list(NULL,rep(maincov,2)))


## Frequencies of each feature given main-variate state
distsFA <- foreach(acov=otherCovs, .inorder=T)%dopar%{
    dists <- rbind(samplesF(Y=grids[[acov]], X=rbind(xcond[,1]), parmList=parmList, inorder=T),
                   samplesF(Y=grids[[acov]], X=rbind(xcond[,2]), parmList=parmList, inorder=T)
                   )
    dim(dists) <- c(length(grids[[acov]]), 2, ncol(dists))
    dimnames(dists) <- list(NULL, maincovnames, NULL)
    aperm(dists, c(3,1,2))
}
names(distsFA) <- otherCovs

## quantiles
qdistsFA <- foreach(acov=otherCovs)%do%{
    apply(distsFA[[acov]],c(2,3),function(x)c(mean(x,na.rm=T), quant(x=x, probs=c(1,15)/16, na.rm=T)))
}
names(qdistsFA) <- otherCovs

test <- 'rh_WM_postcentral'
##
tplot(grids[[test]], t(distsFA[[test]][seq(1,nrow(distsFA[[test]]),length.out=128),,1]),lty=1,col=blue,alpha=0.75)

## plot of frequencies of features given main variate
pdff(paste0(dirname,'/','plots_features_given_Group'), 'a4')
for(acov in otherCovs){
    agrid <- grids[[acov]]
    ymax <- quant(apply(qdistsFA[[acov]],2,function(x){quant(x,99/100)}),99/100)
ylim <- c(0,ymax)#max(qdistsFA[[acov]]))
for(i in 1:2){
    tplot(x=agrid, y=qdistsFA[[acov]][1,,i],
          col=i, lty=i, lwd=4, alpha=0.25, ylim=ylim,
          xlab=acov, ylab='frequency of feature for patients in different Group', add=(i==2))
    ##     tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    ## if(!any(is.na(tpar))){
    ##     Ogrid <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
    ##     axis(3,at=(log(Ogrid)-tpar['transfM'])/tpar['transfW'],labels=Ogrid,lwd=0,lwd.ticks=1,col.ticks='#bbbbbb80')
    ## }
    polygon(x=c(agrid,rev(agrid)), y=c(qdistsFA[[acov]][2,,i], rev(qdistsFA[[acov]][3,,i])), col=paste0(palette()[i],'40'), border=NA)
}
legend(x=agrid[1], y=ylim[2]*1.2, legend=c(paste0('distribution for patients in Group ',subgroupnames), '87.5% uncertainty'),
       col=palette()[c(1,2,7)], lty=c(1,2,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
                       )
}
dev.off()


pdff('testhists')
for(i in names(alldata)){
    if(!(covTypes[i]=='real')){nn <- 'i'}else{nn <- NULL}
    hh <- thist(alldata[[i]], n=nn)
    tplot(hh$breaks,hh$density,xlab=i,ylim=c(0,NA))
}
dev.off()


##
grids <- foreach(acov=covNames)%do%{
    rg <- range(alldata[[acov]], na.rm=T)
    if(acov %in% realCovs){
        rg <- rg+c(-1,1)*IQR(alldata[[acov]],type=8,na.rm=T)/2
        Xgrid <- seq(rg[1], rg[2], length.out=256)
    }else if(acov %in% integerCovs){
            rg <- round(c((covMins[acov]+3*rg[1])/4, (covMaxs[acov]+3*rg[2])/4))
            Xgrid <- rg[1]:rg[2]
    }else{Xgrid <- 0:1}
    matrix(Xgrid, ncol=1, dimnames=list(NULL,acov))
}
names(grids) <- covNames
##
xcond <- matrix(0:1,ncol=2,dimnames=list(NULL,rep(maincov,2)))

## Frequencies of each feature given main-variate value
distsFA <- foreach(acov=otherCovs)%dopar%{
    dists <- rbind(samplesF(Y=grids[[acov]], X=rbind(xcond[,1]), parmList=parmList, inorder=T),
                   samplesF(Y=grids[[acov]], X=rbind(xcond[,2]), parmList=parmList, inorder=T)
                   )
    dim(dists) <- c(length(grids[[acov]]), 2, ncol(dists))
    dimnames(dists) <- list(NULL, maincovnames, NULL)
    aperm(dists, c(3,1,2))
}
names(distsFA) <- otherCovs

## quantiles
qdistsFA <- foreach(acov=otherCovs)%dopar%{
    apply(distsFA[[acov]],c(2,3),function(x)c(mean(x,na.rm=T), quant(x=x, probs=c(1,15)/16, na.rm=T)))
}
names(qdistsFA) <- otherCovs

## plot of frequencies of features given AD state f(F|AD)
pdff(paste0(dirname,'/','allplots_features_given_Group'),'a5')
##par(mfrow=c(3,4), mai=c(0,0,0,0),oma=c(0,0,0,0))
#par(mfrow=c(3,4), mai=c(0,0,0,0),oma=c(0,0,0,0),mar=c(0,0,0,0))
for(acov in otherCovs){
    agrid <- grids[[acov]]
    ymax <- quant(apply(qdistsFA[[acov]],2,function(x){quant(x,99/100)}),99/100)
    ylim <- c(0,ymax)#max(qdistsFA[[acov]]))
    xlim <- c(NA,NA)
    ## tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    ## if(!any(is.na(tpar))){
    ##     xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
    ##     xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    ## }else{xticks <- NULL
    ##     xlabels <- TRUE}
    xticks <- NULL
    xlabels <- TRUE
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    for(i in 1:2){
        tplot(x=agrid, y=qdistsFA[[acov]][1,,i], yticks=NULL, xlim=xlim,
              col=i, lty=i, lwd=3, alpha=0.25, ylim=ylim, xticks=xticks, xlabels=xlabels,ylabels=F,cex.axis=1.25,cex.lab=1.5,ly=0,
              xlab=acov, ylab=NA, add=(i==2))
        polygon(x=c(agrid,rev(agrid)), y=c(qdistsFA[[acov]][2,,i], rev(qdistsFA[[acov]][3,,i])), col=paste0(palette()[i],'40'), border=NA)
    }
    ## legend(x=agrid[1], y=ylim[2]*1.2, legend=c(paste0('distribution for patients who will have ',diseasenames), '87.5% uncertainty'),
    ##        col=palette()[c(1,2,7)], lty=c(1,2,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
    ##        )
    ## for(i in 1:2){
    ##     histo <- thist(c(data.matrix(alldata[get(maincov)==(i-1),..acov])),n=32)
    ##     tplot(x=histo$breaks,y=histo$density,col=i,border=NA,alpha=0.75,xgrid=F,ygrid=F,add=T)
    ## }
}
dev.off()





tplot(y=distsFA[[4]][,200,1])






maincov <- 'Group'
realCovs <- dimnames(parmList$meanR)[[2]]
integerCovs <- dimnames(parmList$probI)[[2]]
binaryCovs <- dimnames(parmList$probB)[[2]]
covNames <- c(realCovs, integerCovs, binaryCovs)
##



if(!setequal(variateinfo$variate, covNames)){
    print('WARNING: mismatch in variate sets')
}
temp <- data.table();
for(i in covNames){temp <- rbind(temp,variateinfo[variate==i])}
variateinfo <- temp
rm('temp')
##
nrcovs <- length(realCovs)
nicovs <- length(integerCovs)
nbcovs <- length(binaryCovs)
ncovs <- length(covNames)
covTypes <- variateinfo$type
covMins <- variateinfo$min
covMaxs <- variateinfo$max
names(covTypes) <- names(covMins) <- names(covMaxs) <- variateinfo$variate
otherCovs <- setdiff(covNames, maincov)
subgroupnames <- c('0', '1')
##

frequenciesfile <- '_mcmctest-V13-D200-K64-I2048/_frequencies-R_mcmctest_149_1-V13-D200-K64-I2048.rds'
parmList <- readRDS(frequenciesfile)
nclusters <- ncol(parmList$q)
nFsamples <- nrow(parmList$q)









datafile <- 'Cortical_myelination_faux.csv'
alldata <- fread(datafile, sep=',')
alldata <- alldata[,..covNames]
## alldata <- alldata[Usage_ == 'train']
##
