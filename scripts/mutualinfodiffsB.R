## Author: PGL  Porta Mana
## Created: 2022-10-07T12:13:20+0200
## Last-Updated: 2023-01-01T20:38:56+0100
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

if(!file.exists(paste0('condprobsx-',nsamplesMI,'.rds'))){
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
}else{
    condprobsx <- readRDS(paste0('condprobsx-',nsamplesMI,'.rds'))
    xsamples2 <- readRDS(paste0('xsamples2-',nsamplesMI,'.rds'))
}

Xlist <- NULL
##
temp <- as.list(setdiff(variatenames,predictands))
names(temp) <- setdiff(variatenames,predictands)
Xlist <- c(Xlist, temp)
##
cogvars <- c('ANARTERR_neuro', 'AVDEL30MIN_neuro', 'AVDELTOT_neuro', 'CATANIMSC_neuro', 'GDTOTAL_gds', 'RAVLT_immediate', 'TRAASCOR_neuro', 'TRABSCOR_neuro')
genvars <- c('Apoe4_', 'LRHHC_n_long')
demovars <- c('Gender_num_', 'AGE')
temp <- list('cog'=c(cogvars,demovars), 'noncog'=c(genvars,demovars))
Xlist <- c(Xlist, temp)

##
ii <- 0
time0 <- Sys.time()
allcondp2 <- foreach(v=names(Xlist), .combine=cbind)%do%{
    cat(paste0('\n',v,' - est. remaining time: '));print((Sys.time()-time0)/ii*(length(predictors)+1-ii))
    ii <- ii+1
    predictors0 <- Xlist[[v]]
    ##
    condprobsxgy <- samplesFDistribution(Y=xsamples2[,predictands,drop=F],
                                       X=xsamples2[,predictors0,drop=F],
                                       mcsamples=mcsamples, varinfo=varinfo,
                                       jacobian=FALSE, fn=mean)
    ##
    colnames(condprobsxgy) <- v
    ## saveRDS(condprobsxgy,paste0('condprobsxgallminus-',v,'-',nsamplesMI,'.rds'))
    condprobsxgy
}

if(!file.exists(paste0('condprobsxgiveny-',nsamplesMI,'.rds'))){
    saveRDS(allcondp2,paste0('morecondprobsxgiveny2-',nsamplesMI,'.rds'))
}else{
    allcondp <- readRDS(paste0('condprobsxgiveny-',nsamplesMI,'.rds'))
    allcondp <- cbind(allcondp2,allcondp)
    saveRDS(allcondp, paste0('allcondprobsxgiveny-',nsamplesMI,'.rds'))
}

test <- apply(cbind(allcondp,condprobsx),1,function(xxx){all(is.finite(xxx))})
allcondp <- allcondp[test,]
condprobsx <- condprobsx[test,]

mism <- apply(log2(allcondp)-log2(c(condprobsx)),2,function(xxx){
    c(mean(xxx,na.rm=T),
      sd(xxx,na.rm=T)/sqrt(sum(is.finite(xxx)))
      )
})
rownames(mism) <- c('mean','sd')

midiff <- apply(mism,2,function(xxx){c(mism[1,1]-xxx[1], abs(mism[2,1]/mism[1,1]+xxx[2]/xxx[1])*xxx[1])})
rownames(midiff) <- c('mean','sd')

cat('\n\nMI:','\n')
print(t(signif(mism[,order(mism[1,])],c(5,1))))

cat('\n\nMI differences:','\n')
print(signif(midiff[,order(midiff[1,])],c(3,1)))

sink()
stop('None. End of script')

sequ <- 1:100
pdff('plotMI',paper='a4')
mism2 <- sort(mism[1,])
tplot(x=rep(sequ,ncol(mism))[1:ncol(mism)],y=mism2,pch=16,cex=1,col=8,type='p');
text(sequ,mism2,colnames(mism), adj=c(0), cex=0.5)
dev.off()


shortnames <- c(
'GDS',
'Sex',
'APOE4',
'Age',
'ANART',
'TMTA',
'TMTB',
'CFT',
'HC',
'APOE4+HC+Age+Sex',
'RAVLT-rec',
'RAVLT-imm',
'RAVLT-del',
'all minus RAVLT-del',
'all minus RAVLT-imm',
'all minus TMTB',
'all minus RAVLT-rec',
'all minus HC',
'cognitive+Age+Sex',
'all minus TMTA',
'all minus CFT',
'all minus Age',
'all minus GDS',
'all minus ANART',
'all minus Sex',
'all minus APOE4',
'all'
)

shortnames2 <- shortnames
shortnames2[order(mism[1,])] <- shortnames

cbind(names(sort(mism[1,])),
      shortnames)

cbind(colnames(mism), shortnames2)


sequ <- 1:ncol(mism)
mord <- order(mism[1,])
mism2 <- mism[1,mord]
errs <- mism[2,mord]
dist <- max(errs)
labs <- shortnames2[mord]
pdff('plotMI',paper='a4p')
tplot(y=sequ,x=mism2,pch=16,cex=1,col=8,type='p',
      yticks=NA,ylab=NA,xlab='mutual information/Sh',
      xlim=c(0,0.2))#max(mism2+errs)+0.1))
tplot(x=rbind(pmax(0,mism2-errs),mism2+errs),
y=rbind(sequ,sequ),lty=1,lwd=2,col=8,
add=T)
text(mism2+dist,sequ,labs, adj=c(-0.05,0.5), cex=1.25)
dev.off()



accus <- apply(1*(allcondp>0.5)+0.5*(allcondp==0.5),2,function(xxx){
    c(mean(xxx,na.rm=T),
      sd(xxx,na.rm=T)/sqrt(sum(is.finite(xxx)))
      )
})
rownames(accus) <- c('mean','sd')

cat('\n\nAcc:','\n')
print(t(signif(accus[,ordwier(accus[1,])],c(5,1))))

sequ <- 1:ncol(mism)
mord <- order(mism[1,])
accus2 <- accus[1,mord]
errs <- accus[2,mord]
dist <- max(errs)
labs <- shortnames2[mord]
pdff('plotAcc',paper='a4p')
tplot(y=sequ[1:length(accus2)],x=accus2,pch=16,cex=1,col=8,type='p',
      yticks=NA,ylab=NA,xlab='accuracy',
      xlim=c(0.5,0.75))#max(accus2)+0.1))
tplot(x=rbind(pmax(0,accus2-errs),accus2+errs),
y=rbind(sequ,sequ),lty=1,lwd=2,col=8,
add=T)
text(accus2+dist,sequ,labs, adj=c(-0.05,0.5), cex=1.25)
dev.off()
m

cbind((1-0.53)mmm*mism[1,mord]+0.5, accus[1,mord])

(1-mi)*0.5+mi
0.5-0.5*mi+mi
0.5+

##testo <- (accus[1,]-min(accus[1,]))/diff(range(accus[1,]))+(mism[1,]-min(mism[1,]))/diff(range(mism[1,]))
tord <- order(mism[1,])
accus2 <- accus[1,tord]
sequ <- 1:100
errs <- accus[2,tord]
pdff('plotAcc',paper='a4p',asp=0.5/sqrt(2))
tplot(y=sequ[1:length(accus2)],x=accus2,pch=16,cex=1,col=8,type='p',
      yticks=NA,ylab=NA,
      xlim=c(NA,max(accus2)+0.1))
tplot(x=rbind(pmax(0,accus2-errs),accus2+errs),
y=rbind(sequ[1:length(accus2)],sequ[1:length(accus2)]),lty=1,lwd=2,col=8,
add=T)
text(accus2+errs,sequ[1:length(accus2)],names(accus2), adj=c(-0.05,0.5), cex=0.75)
dev.off()

cbind(sort((1-mism[1,]/mism[1,'all'])*100))
cbind(sort((1-accus[1,]/accus[1,'all'])*100))

cbind(
Acc=(1-accus[1,]/accus[1,'all'])*100
)

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
