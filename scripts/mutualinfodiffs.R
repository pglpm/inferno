## Author: PGL  Porta Mana
## Created: 2022-10-07T12:13:20+0200
## Last-Updated: 2022-12-30T19:50:52+0100
################
## Combine multiple Monte Carlo chains
################
if(!exists('tplot')){source('~/work/pglpm_plotfunctions.R')}

#rm(list=ls())

nsamplesMI <- 4096*4
outputdir <- NA
extratitle <- 'mutualinfo'
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
sink('MIoutput.txt')

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


#### CALCULATION OF MUTUAL INFO ####

predictandvalues <- cbind(seq(varinfo[['min']][predictands],varinfo[['max']][predictands],length.out=varinfo[['n']][predictands]))
colnames(predictandvalues) <- predictands

set.seed(101)
xsamples2 <- t(generateVariates(Ynames=c(predictors,predictands), X=NULL,
                                mcsamples=mcsamples, varinfo=varinfo,
                                    n=nsamplesMI)[,,1])
saveRDS(xsamples2,paste0('xsamples2-',nsamplesMI,'.rds'))
##
condprobsx <- samplesFDistribution(Y=xsamples2[,predictands,drop=F],
                                   X=NULL,
                                       mcsamples=mcsamples, varinfo=varinfo,
                                   jacobian=FALSE, fn=mean)
colnames(condprobsx) <- predictands
saveRDS(condprobsx,paste0('condprobsx-',nsamplesMI,'.rds'))
##
ii <- 0
time0 <- Sys.time()
allcondp <- foreach(v=c('',predictors), .combine=cbind)%do%{
    cat(paste0('\n',v,' - est. remaining time: '));print((Sys.time()-time0)/ii*(length(predictors)+1-ii))
    ii <- ii+1
    predictors0 <- setdiff(predictors,v)
    ##
    condprobsxgy <- samplesFDistribution(Y=xsamples2[,predictands,drop=F],
                                       X=xsamples2[,predictors0,drop=F],
                                       mcsamples=mcsamples, varinfo=varinfo,
                                       jacobian=FALSE, fn=mean)
    ##
    colnames(condprobsxgy) <- (if(v!=''){v}else{'all'})
    ## saveRDS(condprobsxgy,paste0('condprobsxgallminus-',v,'-',nsamplesMI,'.rds'))
    condprobsxgy
}
saveRDS(allcondp,paste0('condprobsxgiveny-',nsamplesMI,'.rds'))


mism <- apply(log2(allcondp)-log2(c(condprobsx)),2,function(xxx){
    c(mean(xxx,na.rm=T),
      sd(xxx,na.rm=T)/sqrt(sum(is.finite(xxx)))
      )
})
rownames(mism) <- c('mean','sd','O1','O7')

midiff <- apply(mism,2,function(xxx){c(mism[1,1]-xxx[1], abs(mism[2,1]/mism[1,1]+xxx[2]/xxx[1])*xxx[1])})
rownames(midiff) <- c('mean','sd')

cat('\n\nMI:','\n')
print(signif(mism[,order(mism[1,])],c(3,1)))

cat('\n\nMI differences:','\n')
print(signif(midiff[,order(midiff[1,])],c(3,1)))

sink()
stop('None. End of script')

## colnames(MIdata) <- c('all',predictors)
## saveRDS(MIdata,paste0('MIdata-',nsamplesMI,'.rds'))


## decreases <- (colMeans(MIdata)/mean(MIdata[,1])-1)*100
## variances <- (
##     abs(apply(MIdata,2,sd)/mean(MIdata[,1]) -
##     colMeans(MIdata)*sd(MIdata[,1])/mean(MIdata[,1])^2)
##         ) *100/sqrt(nsamplesMI)
## sorto <- order(decreases)
## signif(cbind(
##     decreases[sorto],
##     variances[sorto]
## ),4)
## cbind(
##     round(decreases[sorto],signif(-log10(variances[sorto]),1)),
##     signif(variances[sorto],2)
## )


## stop('End of script')
