## Author: PGL  Porta Mana
## Created: 2022-10-07T12:13:20+0200
## Last-Updated: 2023-01-19T21:12:52+0100
################
## Combine multiple Monte Carlo chains
################
if(!exists('tplot')){source('~/work/pglpm_plotfunctions.R')}

#rm(list=ls())

outputdir <- 'inferencep4'
extratitle <- 'longrunprobs'
## totsamples <- 4096L
datafile <- 'ingrid_data_nogds6.csv'
predictorfile <- 'predictors.csv'
predictandfile <- NULL # 'predictors.csv'
functionsfile <- 'functionsmcmc_2212120902.R'
showdata <- TRUE # 'histogram' 'scatter' FALSE
plotmeans <- TRUE
quantilestoshow <- c(1,7)/8# c(1,31)/32
ndata <- NULL # set this if you want to use fewer data
shuffledata <- FALSE # useful if subsetting data
family <- 'Palatino'


## load customized plot functions
if(!exists('tplot')){source('~/work/pglpm_plotfunctions.R')}
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
if(Sys.info()['nodename']=='luca-HP-Z2-G9'){
    ncores <- 20}else{
    ncores <- 6}
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

#### Load files with variate information and data ####
if(!exists('outputdir') || is.na(outputdir) || is.null(outputdir)){
    outputdir <- as.character(commandArgs(trailingOnly=TRUE))[1]
    if(is.na(outputdir)){
        outputdir <- './'
    }
}
outputdir <- paste0(sub('(.+)/','\\1',outputdir),'/')
setwd(outputdir)
origdir <- '../'
source(paste0(origdir,functionsfile)) # load functions for post-MCMC

varinfofile <- list.files(pattern='^_varinfo.*\\.rds')
if(length(varinfofile) !=1){stop('Problems with varinfo file')}
varinfo <- readRDS(varinfofile)
##
variate <- lapply(variatetypes, function(x)names(varinfo[['type']])[varinfo[['type']]==x])
len <- lapply(variate,length)
names(variate) <- names(len) <- variatetypes

variatenames <- unlist(variate)
if(!is.null(predictorfile)){predictorfile <- paste0(origdir,predictorfile) }
if(!is.null(predictandfile)){predictandfile <- paste0(origdir,predictandfile) }
if(!is.null(predictorfile) && !is.null(predictandfile)){
    predictors <- as.vector(unlist(read.csv(predictorfile, header=F)))
    predictands <- as.vector(unlist(read.csv(predictandfile, header=F)))
}else if(!is.null(predictorfile) && is.null(predictandfile)){
    predictors <- as.vector(unlist(read.csv(predictorfile, header=F)))
    predictands <- setdiff(unlist(variate), predictors)
}else if(is.null(predictorfile) && !is.null(predictandfile)){
    predictands <- as.vector(unlist(read.csv(predictandfile, header=F)))
    predictors <- setdiff(unlist(variate), predictands)
}else{warning('predictors and predictands both missing')}

mcsamplesfile <- list.files(pattern='^_jointmcsamples.*\\.rds')
if(length(mcsamplesfile) !=1){stop('Problems with mcsamples file')}
mcsamples <- readRDS(paste0(mcsamplesfile))

## data0 <- fread(paste0(origdir,datafile), sep=',')
## if(!all(unlist(variate) %in% colnames(data0))){cat('\nERROR: variates missing from datafile')}
## data0 <- data0[, unlist(variate), with=F]
## ## shuffle data
## if(exists('shuffledata') && shuffledata){data0 <- data0[sample(1:nrow(data0))]}
## if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(data0)}
## data0 <- data0[1:ndata]

cogvars <- c('ANARTERR_neuro', 'AVDEL30MIN_neuro', 'AVDELTOT_neuro', 'CATANIMSC_neuro', 'GDTOTAL_gds', 'RAVLT_immediate', 'TRAASCOR_neuro', 'TRABSCOR_neuro')
genvars <- c('Apoe4_', 'LRHHC_n_long')
demovars <- c('Gender_num_', 'AGE')

predictandvalues <- cbind(seq(varinfo[['min']][predictands],varinfo[['max']][predictands],length.out=varinfo[['n']][predictands]))
colnames(predictandvalues) <- predictands
rownames(predictandvalues) <- c('sMCI', 'cAD')


starttime <- Sys.time()
k <- 0
set.seed(6001)
npsamples <- 1024
probsamples <- foreach(asample=seq(1,nrow(mcsamples),length.out=npsamples), .combine=cbind)%do%{
    print(k)
    print('Remaining time:')
    print((Sys.time()-starttime)/k*(npsamples-k))
    k <- k+1
    datasamples <- t(generateVariates(Ynames=c(predictors,predictands), X=NULL,
                                      mcsamples=mcsamples[asample,,drop=F], varinfo=varinfo,
                                      n=npsamples)[,,1])
    yvals <- datasamples[,predictands,drop=F]
    out <- c(samplesFDistribution(Y=yvals,
                                X=datasamples[,predictors,drop=F],
                                mcsamples=mcsamples, varinfo=varinfo,
                                jacobian=TRUE, fn=identity))
    cbind(out,
          ## c(samplesFDistribution2(Y=yvals,
          ##                       X=NULL,
          ##                       mcsamples=asample, varinfo=varinfo,
          ##                       jacobian=TRUE, fn=identity)),
          (c(yvals)==1) * out + (c(yvals)==0) * (1-out)
          )
}
attr(probsamples, 'rng') <- NULL
attr(probsamples, 'doRNG_version') <- NULL
dim(probsamples) <- c(nrow(probsamples),2,ncol(probsamples)/2)
dimnames(probsamples) <- list(NULL,c('p(Y|X)','p(1|X)'),NULL)
probsamples <- aperm(probsamples)
## dim1: mcsamples, dim2:p/Yvalue, dim3: samples given lim-freq
saveRDS(probsamples, paste0('_probsamplesB_',npsamples,'.rds'))

print('done')


probsamples <- readRDS('_probsamplesB_1024.rds')

newpatientprobso <- readRDS('_newpatientMI_32768.rds')
newpatients <- readRDS('_points_newpatientMI_32768.rds')
origp <- newpatientprobso[,'all']*(newpatients[,predictands]==1)+(1-newpatientprobso[,'all'])*(newpatients[,predictands]==0)

pgrid <- seq(0,1,by=0.02)
pgridm <- pgrid[-1]-diff(pgrid)/2
histos <- apply(probsamples,1,function(xxx){
    thist(xxx['p(Y|X)',], n=pgrid)$density
})

pgrid <- seq(0,1,by=0.02)
pgridm <- pgrid[-1]-diff(pgrid)/2
histos <- apply(probsamples,1,function(xxx){
    thist(xxx['p(Y|X)',], n=pgrid)$density
})

subhist <- seq(1,ncol(histos),length.out=200)
pdff('plotnextpatientprob',paper='a4p')
tplot(x=pgridm,y=histos[,subhist], lty=1, lwd=0.5, alpha=0.75, col=7,
      ylab='density',
      xlab=paste0('probability of correct outcome'), cex.lab=1.4,
      ylim=c(0,NA))
tplot(x=pgridm,thist(newpatientprobso[,'all'],n=pgrid)$density, 
col=1,
lty=1, lwd=4, alpha=0.5, add=T)
dev.off()

normal <- function(x){x/sum(x)}

pgrid <- seq(0,1,by=0.02)
pgridm <- pgrid[-1]-diff(pgrid)/2
histos <- apply(probsamples,1,function(xxx){
    normal(thist(xxx['p(1|X)',], n=pgrid)$counts)
})
statlongrunprobs <- apply(histos,1,function(xxx){
c(tquant(xxx, c(8*2.5/100,1,2,4,6,7,8*97.5/100)/8), mean=mean(xxx))
})
##
subhist <- seq(1,ncol(histos),length.out=200)
pdff('plotnextpatientconvprob')#,paper='a4p')
## tplot(x=pgridm,y=histos[,subhist], lty=1, lwd=0.5, alpha=0.75, col=7,
##       ylab='density',
##       xlab=paste0('probability of conversion'), cex.lab=1.4,
##       ylim=c(0,NA))
tplot(pgridm,statlongrunprobs['97.5%',], col='white',
      xlab=expression(italic(P):~'probability of conversion'~'(2% bins)'),
      ylab=expression('fraction of population prognosed with prob.'~italic(P)),
      ly=3,mar=c(NA,4.5,2,NA), cex.axis=1.25,
      ## yticks=seq(0,ceiling(max(statlongrunprobs)),by=1),
      xticks=seq(0,1,by=0.1),#ceiling(max(statlongrunprobs)),by=1),
      xlabels=paste0(seq(0,1,by=0.1)*100,'%')#ceiling(max(statlongrunprobs)),by=1),'%')
      )
plotquantiles(pgridm,t(statlongrunprobs[c('2.5%','97.5%'),]), col=5,alpha=0.5)
##plotquantiles(pgridm,t(statlongrunprobs[c('25%','75%'),]), col=5,alpha=0.75)
## tplot(pgridm,statlongrunprobs['mean',],col=2,lty=2,add=T)
tplot(x=pgridm,normal(thist(origp,n=pgrid)$counts), col=1,
lty=1, lwd=4, alpha=0.25, add=T)
legend('topleft',legend=c(##'median',
                          ##'mean',
                          ##'50% uncertainty','75% uncertainty',
                          '95% uncertainty'),
       col=c(##1,1,
             alpha2hex(5,c(##0.25,0.5,
                             0.5))),
       lwd=c(##2,2,5,10,
             15),
       lty=c(##1,2,1,1,
           1),
       bty='n'
       )
dev.off()

choicefn <- function(pv,umx){
    out <- sapply(c(pv), function(p){
        which.max(c(umx %*% c(1-p,p)))
    })
    dim(out) <- dim(pv)
    dimnames(out) <- dimnames(pv)
    out
}
##
um <- matrix(c(1,0.9,0.8,0,
               0,0.3,0.5,1), 4,2)
rownames(um) <- c(1:4)
##
alltreatdistr <- apply(probsamples,1,function(xxx){
    tabulate(choicefn(xxx['p(1|X)',], um), nbins=4)/ncol(xxx)
})

apply(alltreatdistr,1, tquant, prob=c(1,4,7)/8)*100
##          [,1]    [,2]    [,3]    [,4]
## 12.5% 21.8872 2.53906 33.2153 32.0312
## 50%   24.5117 3.32031 37.0117 34.9609
## 87.5% 27.0508 4.39453 41.1133 38.0859

apply(alltreatdistr,1, tquant, prob=c(5,95)/100)*100
##        [,1]    [,2]    [,3]    [,4]
## 5%  20.5078 2.24609 31.2500 30.5664
## 95% 28.1006 5.56641 42.3828 39.8438
