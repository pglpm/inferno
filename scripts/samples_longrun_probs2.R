## Author: PGL  Porta Mana
## Created: 2022-10-07T12:13:20+0200
## Last-Updated: 2023-01-21T18:24:29+0100
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

data0 <- fread(paste0(origdir,datafile), sep=',')
if(!all(unlist(variate) %in% colnames(data0))){cat('\nERROR: variates missing from datafile')}
data0 <- data0[, unlist(variate), with=F]
## shuffle data
if(exists('shuffledata') && shuffledata){data0 <- data0[sample(1:nrow(data0))]}
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(data0)}
data0 <- data0[1:ndata]

cogvars <- c('ANARTERR_neuro', 'AVDEL30MIN_neuro', 'AVDELTOT_neuro', 'CATANIMSC_neuro', 'GDTOTAL_gds', 'RAVLT_immediate', 'TRAASCOR_neuro', 'TRABSCOR_neuro')
genvars <- c('Apoe4_', 'LRHHC_n_long')
demovars <- c('Gender_num_', 'AGE')

predictandvalues <- cbind(seq(varinfo[['min']][predictands],varinfo[['max']][predictands],length.out=varinfo[['n']][predictands]))
colnames(predictandvalues) <- predictands
rownames(predictandvalues) <- c('sMCI', 'cAD')

um <- matrix(c(1,0.9,0.8,0,
               0,0.3,0.5,1), 4,2)
rownames(um) <- c(1:4)

Xlist <- NULL
##
temp <- as.list(predictors)
names(temp) <- predictors
Xlist <- c(Xlist, temp)
for(avar in predictors){
temp <- list(setdiff(predictors,avar))
names(temp) <- paste0('all_minus_',avar)
Xlist <- c(Xlist, temp)
}
temp <- list('cog'=c(cogvars,demovars), 'noncog'=c(genvars,demovars))
Xlist <- c(Xlist, temp, list(all=predictors))

set.seed(6001)
points <- rbind(t(generateVariates(Ynames=c(predictors,predictands),
                             X=NULL,
                             mcsamples=mcsamples,
                             varinfo=varinfo,
                             n=nrow(mcsamples))[,,1]) ,
                t(generateVariates(Ynames=c(predictors,predictands),
                             X=NULL,
                             mcsamples=mcsamples,
                             varinfo=varinfo,
                             n=nrow(mcsamples))[,,1])
                )

meanrm <- function(x){mean(x, na.rm=T)}
##
starttime <- Sys.time()
k <- 0
newpatientMIs <- foreach(preds=Xlist, .combine=cbind)%do%{
    print(preds)
    print('Remaining time:')
    print((Sys.time()-starttime)/k*(length(Xlist)-k))
    k <- k+1
    samplesFDistribution(Y=points[,predictands,drop=F],
                                          X=points[,preds,drop=F],
                                          mcsamples=mcsamples, varinfo=varinfo,
                                          jacobian=TRUE, fn=meanrm)
}
colnames(newpatientMIs) <- names(Xlist)

saveRDS(newpatientMIs, paste0('newpatientMI_',nrow(points)))


## ## rows: patients, cols: longrunfreqs
## longrunprobsamples <-samplesFDistribution(Y=points[,predictands,drop=F],
##                                           X=points[,predictors,drop=F],
##                                           mcsamples=mcsamples, varinfo=varinfo,
##                                           jacobian=TRUE)

## saveRDS(longrunprobsamples, paste0('longrunprobsamplesB',npsamples,'.rds'))


## ## npsamples <- 1024
## ## longrunprobsamples <- readRDS(paste0('_longrunprobsamples',npsamples,'.rds'))

## longrunprobs2 <-samplesFDistribution(Y=points[,predictands,drop=F],
##                                           X=points[,predictors,drop=F],
##                                           mcsamples=mcsamples, varinfo=varinfo,
##                                           jacobian=TRUE, fn=meanrm)


print('done')


shortnames <- c(
'TMTA',
'TMTB',
'Age',
'HV',
'APOE4',
'Sex',
'ANART',
'RAVLT-del',
'RAVLT-rec',
'CFT',
'GDS',
'RAVLT-imm'
)

shortnamesacc <- c(
    shortnames,
    paste0('all minus ', shortnames),
'cognitive+Age+Sex',
'APOE4+HC+Age+Sex',
    'all'
)
names(shortnamesacc) <- names(Xlist)

allshortnames <- c(Subgroup_num_='cAD', shortnamesacc)

saveRDS(allshortnames, 'variatenames.rds')

newpatientprobso <- readRDS('_newpatientMI_32768.rds')
newpatients <- readRDS('_points_newpatientMI_32768.rds')
newpatientprobs <- newpatientprobso

meansnpp <- colMeans(newpatientprobs)
accus <- apply(newpatientprobs,2,function(xxx){
    mean((xxx>0.5)+0.5*(xxx==0.5))
})
orderv <- accus
ord <- order(orderv, decreasing=T)
ord <- ord[-which(ord==which(names(orderv)=='all'))]
##
pgrid <- seq(0,1,by=0.02)
histos <- apply(newpatientprobs,2,function(probs){
    thist(probs,n=pgrid)$density
})
##
pdff('newpatient_probcorrect')
for(apred in colnames(histos)[ord]){
    apredn <- shortnamesacc[apred]
ymax <- max(histos[,c('all',apred)])*1.05
tplot(x=pgrid[-1]-diff(pgrid)/2,y=histos[,'all'], 
      ylab='density',
      xlab=paste0('probability for true predictand, given predictor: ',apredn), cex.lab=1.4,
ylim=c(0,ymax), col=1,
      lty=1, lwd=6, alpha=0.5)
tplot(x=pgrid[-1]-diff(pgrid)[1]/2,y=histos[,apred],
      col=6, lwd=3, lty=5, alpha=0.5, add=T)
    abline(v=c(accus['all'],accus[apred]),
           col=alpha2hex2(0.25,c(5,2)), lty=c(1,5), lwd=c(3,2))
    legend('topleft',legend=c(
                         paste0('all (acc.: ',round(round(100*accus['all']*1,1)/1,1),'%)'),
                         paste0(apredn,' (acc.: ',round(round(100*accus[apred]*1,1)/1,1),'%)')),
       bty='n', col=c(1,6,7), lwd=c(6,3), box.lwd=, lty=c(1,2,NA), density=c(0,0,20), border=c(NA,NA,1), angle=90)
}
dev.off()

origp <- newpatientprobso[,'all']*(newpatients[,predictands]==1)+(1-newpatientprobso[,'all'])*(newpatients[,predictands]==0)
##
pdff('plotnextpatientprob',paper='a4p')
tplot(x=pgrid[-1]-diff(pgrid)/2,y=thist(origp,n=pgrid)$density, 
      ylab='density',
      xlab=paste0('probability of conversion'), cex.lab=1.4,
ylim=c(0,NA), col=1,
lty=1, lwd=4, alpha=0)
dev.off()



phist <- t(apply(probAD2,1,function(xxx){
    thist(xxx,n=seq(0,1,by=0.025))$density
}))

xgrid <- thist(probAD2[1,],n=seq(0,1,by=0.025))$mids

Ophist <- apply(phist,2,function(xxx){tquant(xxx, c(1,4,7)/8)})

subsett <- round(seq(1,nrow(phist),length.out=100))
pdff('histog_future_probs')
tplot(x=xgrid,y=t(phist[subsett,]), xlim=0:1,ylim=c(0,NA),
      col=7,lty=1,lwd=1,
      ylab='frequency density',
      xlab='prognostic probability for future patients')
dev.off()













meanprobAD2 <- newpatientprobs

choicefn <- function(pv,umx){
    out <- sapply(c(pv), function(p){
        which.max(c(umx %*% c(1-p,p)))
    })
    dim(out) <- dim(pv)
    dimnames(out) <- dimnames(pv)
    out
}

treat <- choicefn(probAD2, um)
treatmean <- choicefn(meanprobAD2, um)


ordering <- order(meanprobAD2)
subsett <- round(seq(1,length(meanprobAD2),length.out=50))
##syms <- (c(1,3,2,4))
##syms <- as.character(1:4)
syms <- c(expression(italic(alpha)), expression(italic(beta)), expression(italic(gamma)), expression(italic(delta)))
pdff('plotnextpatientprob')
tplot(x=-10,y=-10,xlim=c(0,length(meanprobAD2)),ylim=0:1,
            xlab='sample patient #', ylab='probability of conversion to AD')
## plotquantiles(x=1:length(meanprobAD2),y=OprobAD2[ordering,], col=6,alpha=0.75)
tplot(y=meanprobAD2[ordering], ylim=0:1, col=2, lwd=3,
      xlab='sample patient #', ylab='probability of conversion to AD',
      add=T)
text(x=subsett, y=meanprobAD2[ordering][subsett],
     syms[treatmean[ordering][subsett]], cex=1.25, col=alpha2hex(8,0.25))
## tplot(x=as.list(subsett),y=as.list(meanprobAD2[ordering][subsett]), type='p',
##       pch=syms[treatmean[ordering][subsett]],
##       col=1, cex=1.25,
##       add=T)
## plotquantiles(x=1:nrow(datapoints2),y=CprobAD2[ordering,], col=5,alpha=0.75)
dev.off()




probAD2 <- newpatientprobs

meanprobAD2 <- rowMeans(probAD2,na.rm=T)
OprobAD2 <- t(apply(probAD2,1,function(xxx)tquant(xxx,c(1,7)/8)))
CprobAD2 <- t(apply(probAD2,1,function(xxx)tquant(xxx,c(5,95)/100)))
accprobAD2 <- sapply(meanprobAD2,function(xxx)max(xxx,1-xxx))

ordering <- order(meanprobAD2)
pdff('plotnextpatientprob')
tplot(y=meanprobAD2[ordering], ylim=0:1, col=2, lwd=3,
      xlab='sample patient #', ylab='probability of conversion to AD')
## plotquantiles(x=1:nrow(datapoints2),y=CprobAD2[ordering,], col=5,alpha=0.75)
plotquantiles(x=1:length(ordering),y=OprobAD2[ordering,], col=6,alpha=0.75)
dev.off()

