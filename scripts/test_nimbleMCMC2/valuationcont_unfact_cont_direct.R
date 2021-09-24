## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-09-18T07:18:46+0200
################
## Script for evaluation of regression:
## Unfactorizable prior
## r|x exchangeable
## x not exchangeable
## r continuous
################
#### Custom setup ####
## Colour-blind friendly palettes, from https://personal.sron.nl/~pault/
library('khroma')
palette(colour('bright')())
## palette(colour('muted')())
library('data.table')
library('ggplot2')
library('ggthemes')
theme_set(theme_bw(base_size=18))
scale_colour_discrete <- scale_colour_bright
#library('cowplot')
library('png')
library('foreach')
library('doFuture')
library('doRNG')
registerDoFuture()
#library('LaplacesDemon') # used for Dirichlet generator
library('ash')
library('extraDistr')
library('PReMiuM')
library('mvtnorm')
options(bitmapType='cairo')
library('nimble')

load(file='_nimbleoutputv2.RData')

##################################################################
##################################################################
## Function to calculate \sum_t F^{(t)}_{r|x}
##
source('functions_rmsdregr_nimble.R')
boundaries <- log(c(2,3))
rangeR <- range(alldata$log_RMSD, na.rm=T)
rangeR <- (rangeR - mean(rangeR)) * 1.1 + mean(rangeR)
lRgrid <- 101
Rgrid <- seq(rangeR[1],rangeR[2],length.out=lRgrid)

ndata <- 6000
seldata <- 1:ndata
dset <- 'tail'
basepriors <- rbind(asdata=normalize(as.vector(table(alldata$bin_RMSD))), uniform=rep(1,3)/3)
## Compare with a model based only on prior frequencies of r
unseldata <- setdiff(1:nrow(alldata), seldata)
nTest <- 2000
## construct test set
testdata <- alldata[do.call(dset,list(x=unseldata,n=3*nTest)), c('bin_RMSD',covNames), with=F]
## testd <- data[unseldata, covNames, with=F]
## testdata <- data.table()
## for(val in rmsdVals){
##     testdata <- rbind(testdata, tail(testd[bin_RMSD==val],n=round(3*nTest*priorR[val])))
## }
## rm(testd)
##
gc()
plan(sequential)
plan(multisession, workers = 6L)
condfreqs <- expeRgivenX(as.matrix(testdata), parmList)
plan(sequential)

##
## Utility functions
## function to compute utilities in test set
lossfunction <- function(truevalues,choices){-(truevalues-choices)^2}
##
gains <- lossfunction(testdata$log_RMSD, condfreqs)
chancegains <- lossfunction(testdata$log_RMSD, log(2.5))
##
basescores <- data.table(model=gains, chance=chancegains)
##
pdff(paste0('quadraticscores_direct_nonfact_conditional-run',version,'_P_',dset))
##
print(
    ggplot(dt <- melt(basescores[,1:2])) +
    geom_histogram(aes(x=value, fill=variable), color='white', bins=10, alpha=0.5, position='identity') + scale_fill_bright() + #ylim(c(0,nTest*3)) +
    xlab('quadratic gain') + theme(legend.position='top') +
    geom_vline(data=dt[,.(value=mean(value)), by=variable], aes(xintercept=value, color=variable, linetype=variable), alpha=0.5, size=2)
)
dev.off()

stop('End')


####
## Utility matrices
dgain <- diag(1,3)
cgain <- 1-sapply(1:3,function(x){abs(x-1:3)})/2
## function for randomly choosing ties
resample <- function(x, ...) x[sample.int(length(x), ...)]
## function to compute utilities in test set
utilities <- function(truevalues,probX,utilitym){
    if(is.null(dim(probX))){probX <- t(matrix(probX,nrow=length(probX),ncol=length(truevalues)))}
    y <- apply(t(tcrossprod(utilitym, normalizem(probX))), 1,
               function(z){resample(which(z==max(z)))})
    diag(utilitym[truevalues,y])
}
## function to compute various prob-averages in test set
probscores <- function(truevalues,probX,meanfunction){
    if(is.null(dim(probX))){probX <- t(matrix(probX,nrow=length(probX),ncol=length(truevalues)))}
    y <- normalizem(probX)
    meanfunction(diag(y[,truevalues]))
}
## function to construct data table with results
metrics <- function(truevalues, probX, chanceprior){
    list(
        delta_gain = data.table(
                 model=utilities(truevalues,probX,dgain),
                 chance=utilities(truevalues,chanceprior,dgain),
                 range=range(dgain)
             ),
         ##
        contig_gain = data.table(
                 model=utilities(truevalues,probX,cgain),
                 chance=utilities(truevalues,chanceprior,cgain),
                 range=range(cgain)
             ),
         ##
        log_score = data.table(
                 model=probscores(truevalues,probX,log),
                 chance=probscores(truevalues,chanceprior,log),
                 range=c(-Inf,0)
             ),
         ##
        mean_score = data.table(
                 model=probscores(truevalues,probX,identity),
                 chance=probscores(truevalues,chanceprior,identity),
                 range=c(0,1)
             )
        )}
##
##
#### Calculation of utilities and scores for test set
## source(file='calibration_plots.R')
##

case <- 1 # 'fdata'
## Frequencies of RMSD categories in test set are same as full dataset.
## It does not make sense to change them, because the conditional frequencies r|x
## then also change - we can then make the model have a score as bad as we please
##
scores <- list()
basepriors <- rbind(asdata=normalize(as.vector(table(data$bin_RMSD))), uniform=rep(1,3)/3)
for(dset in c('tail', 'head')){
    for(ii in 1:nrow(basepriors)){
        priorR <- basepriors[ii,]
    ## Compare with a model based only on prior frequencies of r
    unseldata <- setdiff(1:nrow(data), seldata)
    nTest <- 2000
    ## construct test set
    testdata <- data[do.call(dset,list(x=unseldata,n=3*nTest)), c('bin_RMSD',covNames), with=F]
    ## testd <- data[unseldata, covNames, with=F]
    ## testdata <- data.table()
    ## for(val in rmsdVals){
    ##     testdata <- rbind(testdata, tail(testd[bin_RMSD==val],n=round(3*nTest*priorR[val])))
    ## }
    ## rm(testd)
    ##
    gc()
    plan(sequential)
    plan(multisession, workers = 6L)
        condfreqscont <- predictRmeanX(MCMCdata[[case]], testdata)[,1]
        condfreqs <- 
        
        plan(sequential)
    ##
    set.seed(247)
    scores[[dset]][[rownames(basepriors)[ii]]] <- metrics(truevalues=testdata[,bin_RMSD], probX=condfreqs, chanceprior=priorR)
    print(dset)
    pdff(paste0('scores_direct_nonfact_conditional_P',rownames(basepriors)[ii],'_',dset))
    basescores <- scores[[dset]][[ii]]
    ##
        print(
            ggplot(dt <- melt(basescores$delta_gain[,1:2])) +
        geom_histogram(aes(x=value, fill=variable), color='white', bins=length(unique(c(dgain))), alpha=0.5, position='identity') + scale_fill_bright() + ylim(c(0,nTest*3)) +
        xlab('delta utility') +
        geom_vline(data=dt[,.(value=mean(value)), by=variable], aes(xintercept=value, color=variable, linetype=variable), alpha=0.5, size=2)
        )
    ##
print(    ggplot(dt <- melt(basescores$contig_gain[,1:2])) +
        geom_histogram(aes(x=value, fill=variable), color='white', bins=length(unique(c(cgain))), alpha=0.5, position='identity') + scale_fill_bright() + ylim(c(0,nTest*3)) + 
        xlab('contiguous utility') +
        geom_vline(data=dt[,.(value=mean(value)), by=variable], aes(xintercept=value, color=variable, linetype=variable), alpha=0.5, size=2))
    ##
        print(
            ggplot(dt <- melt(basescores$log_score[,1:2])) +
        geom_histogram(aes(x=value, y=..density.., fill=variable), color='white', bins=10, alpha=0.5, position='identity') + scale_fill_bright() + ylim(c(0,2.5)) +
        xlab('log-probability') + xlim(c(-4,0)) +
        geom_vline(data=dt[,.(value=mean(value)), by=variable], aes(xintercept=value, color=variable, linetype=variable), alpha=0.5, size=2)
        )
    ##
print(    ggplot(dt <- melt(basescores$mean_score[,1:2])) +
        geom_histogram(aes(x=value, y=..density.., fill=variable), color='white', bins=10, alpha=0.5, position='identity') + scale_fill_bright() + ylim(c(0,7.5+2.5*3/4))+
        xlab('probability') + xlim(c(0,1)) +
        geom_vline(data=dt[,.(value=mean(value)), by=variable], aes(xintercept=value, color=variable, linetype=variable), alpha=0.5, size=2))
    ##
    dev.off()
    }}
save(list=c('scores'), file=paste0('scores_direct_continuousr_conditional_nonfact_T',nTest*3,'_N',ndata,'_',length(covNums),'covs.RData'))



