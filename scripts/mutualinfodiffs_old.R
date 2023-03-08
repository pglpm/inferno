## Author: PGL  Porta Mana
## Created: 2022-10-07T12:13:20+0200
## Last-Updated: 2022-12-21T14:12:40+0100
################
## Combine multiple Monte Carlo chains
################
if(!exists('tplot')){source('~/work/pglpm_plotfunctions.R')}

#rm(list=ls())

outputdir <- 'mutualinforesults2'
extratitle <- 'mutualinfo'
## totsamples <- 4096L
datafile <- 'ingrid_data_nogds6.csv'
predictorfile <- 'predictors.csv'
predictandfile <- NULL # 'predictors.csv'
mainfilelocation <- '_inference1/'
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

#### Load files with MC information and data ####
if(!exists('outputdir') || is.na(outputdir) || is.null(outputdir)){
    outputdir <- as.character(commandArgs(trailingOnly=TRUE))[1]
    if(is.na(outputdir)){
        outputdir <- './'
    }
}
outputdir <- paste0(sub('(.+)/','\\1',outputdir),'/')
dir.create(outputdir)
origdir <- paste0(getwd(), '/')
setwd(outputdir)
source(paste0(origdir,functionsfile)) # load functions for post-MCMC

varinfofile <- paste0(origdir, mainfilelocation,
                      list.files(path=paste0(origdir,mainfilelocation), pattern='^_varinfo.*\\.rds'))
mcsamplesfile <-  paste0(origdir, mainfilelocation,
                         list.files(path=paste0(origdir,mainfilelocation), pattern='^_jointmcsamples.*\\.rds'))

varinfo <- readRDS(paste0(varinfofile))

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

mcsamples <- readRDS(paste0(mcsamplesfile))

data0 <- fread(paste0(origdir,datafile), sep=',')
if(!all(unlist(variate) %in% colnames(data0))){cat('\nERROR: variates missing from datafile')}
data0 <- data0[, unlist(variate), with=F]
## shuffle data
if(exists('shuffledata') && shuffledata){data0 <- data0[sample(1:nrow(data0))]}
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(data0)}
data0 <- data0[1:ndata]


#### CALCULATION OF MUTUAL INFO ####

predictandvalues <- cbind(seq(varinfo[['min']][predictands],varinfo[['max']][predictands],length.out=varinfo[['n']][predictands]))
colnames(predictandvalues) <- predictands

set.seed(101)
nsamplesMI <- 4096
xsamples2 <- t(generateVariates(Ynames=c(predictors,predictands), X=NULL,
                                mcsamples=mcsamples, varinfo=varinfo,
                                    n=nsamplesMI)[,,1])
saveRDS(xsamples2,paste0('xsamples2-',nsamplesMI,'.rds'))
##
condprobsx <- samplesFDistribution(Y=xsamples2[,predictands,drop=F],
                                   X=NULL,
                                       mcsamples=mcsamples, varinfo=varinfo,
                                   jacobian=FALSE, fn=mean)
saveRDS(condprobsx,paste0('condprobsx-',nsamplesMI,'.rds'))
##
MIdata <- foreach(v=c('',predictors), .combine=cbind)%do%{
    print(v)
    predictors0 <- setdiff(predictors,v)
    ##
    condprobsy <- samplesFDistribution(Y=xsamples2[,predictors0,drop=F],
                                       X=NULL,
                                       mcsamples=mcsamples, varinfo=varinfo,
                                       jacobian=FALSE, fn=mean)
    ##
    condprobsj <- samplesFDistribution(Y=xsamples2[,c(predictors0,predictands),drop=F],
                                       X=NULL,
                                       mcsamples=mcsamples, varinfo=varinfo,
                                       jacobian=FALSE, fn=mean)
    ##
    saveRDS(condprobsy,paste0('condprobsy-',v,'-',nsamplesMI,'.rds'))
    saveRDS(condprobsj,paste0('condprobsj-',v,'-',nsamplesMI,'.rds'))
    ##
    log2(condprobsj)-log2(condprobsx)-log2(condprobsy)
}
colnames(MIdata) <- c('all',predictors)
saveRDS(MIdata,paste0('MIdata-',nsamplesMI,'.rds'))

stop('End of script')

nsamplesMI <- nrow(MIdata)
decreases <- (colMeans(MIdata)/mean(MIdata[,'all'])-1)*100
sds <- apply(MIdata,2,sd)/sqrt(nsamplesMI)
variances <- (
    abs(apply(MIdata,2,sd)/mean(MIdata[,'all']) -
    colMeans(MIdata)*sd(MIdata[,'all'])/mean(MIdata[,'all'])^2)
        ) *100/sqrt(nsamplesMI)
sorto <- order(decreases)
signif(cbind(
    decreases[sorto],
    variances[sorto]*4
),4)
cbind(
    round(decreases[sorto],signif(-log10(variances[sorto]*4),1)),
    signif(variances[sorto]*4,2)
)

abs(apply(MIdata,2,sd))/sqrt(nsamplesMI)/colMeans(MIdata)*100

testsamps <- apply(MIdata,2,function(xxx){rnorm(1e4,mean=mean(xxx),sd=sd(xxx)/sqrt(nsamplesMI))})

thist((testsamps[,'GDTOTAL_gds']/testsamps[,'all']-1)*100, plot=T)


signif(cbind(colMeans(MIdata), apply(MIdata,2,sd)/sqrt(nsamplesMI)),3)



### tests

ysamples <- t(generateVariates(Ynames=sample(variatenames), X=NULL,
                                mcsamples=mcsamples, varinfo=varinfo,
                               n=4096*4)[,,1])
##
pdff('testgeneratingfunction')
for(v in variatenames){
    if(varinfo[['type']][v] %in% c('I','B')){
        xgrid <- seq(varinfo[['min']][v], varinfo[['max']][v],
                     length.out=varinfo[['n']][v])
        nh <- (varinfo[['max']][v]-varinfo[['min']][v])/(varinfo[['n']][v]-1)
        nh <- seq(varinfo[['min']][v]-nh/2, varinfo[['max']][v]+nh/2, length.out=varinfo[['n']][v]+1)
    }else{
        xgrid <- seq(varinfo[['plotmin']][v], varinfo[['plotmax']][v],
                     length.out=256)
        nh <- max(10,round(length(ysamples[,v])/64))
    }
    xgrid <- cbind(xgrid)
    colnames(xgrid) <- v
    testpdf <- samplesFDistribution(Y=xgrid,
                                    X=NULL,
                                    mcsamples=mcsamples, varinfo=varinfo,
                                    jacobian=TRUE, fn=mean)
    histo <- thist(ysamples[,v],n=nh)
    tplot(list(histo$mids, xgrid),list(histo$density,testpdf[,1]),ylim=c(0,NA),xlab=v)
}
dev.off()
