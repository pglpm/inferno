## Author: PGL  Porta Mana
## Created: 2022-10-07T12:13:20+0200
## Last-Updated: 2023-01-25T10:41:37+0100
################
## Combine multiple Monte Carlo chains
################
if(!exists('tplot')){source('~/work/pglpm_plotfunctions.R')}

#rm(list=ls())

outputdir <- 'inferencep4'
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

varinfofile <- list.files(pattern='^varinfo\\.rds')
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

#### Read conditional probabilities for entropy quantities ####
nsamplesMI <- 4096*4

## samples of points in 13D space
datapoints <- readRDS(paste0('xsamples2-',nsamplesMI,'.rds'))

## marginal probabilities of predictand for the samples
probsY <- readRDS(paste0('condprobsx-',nsamplesMI,'.rds'))

## conditional probabilities of predictand given other variates, for the samples
probsYgivenX <- readRDS(paste0('allcondprobsxgiveny-',nsamplesMI,'.rds'))

shortnames <- c(
'APOE4',
'Sex',
'ANART',
'RAVLT-del',
'RAVLT-rec',
'CFT',
'GDS',
'RAVLT-imm',
'TMTA',
'TMTB',
'Age',
'HC',
'cognitive+Age+Sex',
'APOE4+HC+Age+Sex',
'all',
'all minus TMTA',
'all minus TMTB',
'all minus Age',
'all minus HC',
'all minus APOE4',
'all minus Sex',
'all minus ANART',
'all minus RAVLT-del',
'all minus RAVLT-rec',
'all minus CFT',
'all minus GDS',
'all minus RAVLT-imm'
)
##      [,1]                       
##  [1,] "Apoe4_"                   
##  [2,] "Gender_num_"              
##  [3,] "ANARTERR_neuro"           
##  [4,] "AVDEL30MIN_neuro"         
##  [5,] "AVDELTOT_neuro"           
##  [6,] "CATANIMSC_neuro"          
##  [7,] "GDTOTAL_gds"              
##  [8,] "RAVLT_immediate"          
##  [9,] "TRAASCOR_neuro"           
## [10,] "TRABSCOR_neuro"           
## [11,] "AGE"                      
## [12,] "LRHHC_n_long"             
## [13,] "cog"                      
## [14,] "noncog"                   
## [15,] "all"                      
## [16,] "allminus_TRAASCOR_neuro"  
## [17,] "allminus_TRABSCOR_neuro"  
## [18,] "allminus_AGE"             
## [19,] "allminus_LRHHC_n_long"    
## [20,] "allminus_Apoe4_"          
## [21,] "allminus_Gender_num_"     
## [22,] "allminus_ANARTERR_neuro"  
## [23,] "allminus_AVDEL30MIN_neuro"
## [24,] "allminus_AVDELTOT_neuro"  
## [25,] "allminus_CATANIMSC_neuro" 
## [26,] "allminus_GDTOTAL_gds"     
## [27,] "allminus_RAVLT_immediate" 

colnames(probsYgivenX) <- shortnames

accnextpatient <- apply(probsYgivenX,2,function(xxx){
    c(mean((xxx > 0.5) + 0.5*(xxx == 0.5), na.rm=T),
      sd((xxx > 0.5) + 0.5*(xxx == 0.5), na.rm=T)*2/sqrt(nsamplesMI)
      )
})


t(accnextpatient[,ordering])

ordering <- c(setdiff(order(accnextpatient[1,]),which(colnames(accnextpatient)=='all')),
              which(colnames(accnextpatient)=='all'))
pdff('accuracy_ranks', paper='a4p')
tplot(x=accnextpatient[1,ordering],y=1:ncol(accnextpatient), type='p', pch=16,
      xlim=c(0.5,0.785), cex=1,
      xlab='accuracy (expected utility)',
      col=8, mar=c(NA,1,1,1),
 xticks=seq(0.5,0.8, by=0.02),#cex.axis=1,
      yticks=NA,ylab=NA)
tplot(x=rbind(accnextpatient[1,ordering]-accnextpatient[2,ordering],
              accnextpatient[1,ordering]+accnextpatient[2,ordering]),
      y=rbind(1:ncol(accnextpatient), 1:ncol(accnextpatient)),
      col=8, lty=1, lwd=2,
      add=T)
text(x=accnextpatient[1,ordering]+accnextpatient[2,ordering]+0.002,
     y=1:ncol(accnextpatient),
     colnames(accnextpatient)[ordering], adj=c(0,0.5),
     cex=1.5)
dev.off()


minextpatient <- apply(probsYgivenX,2,function(xxx){
    xxx <- log2(xxx)-log2(probsY[,1])
    c(mean(xxx, na.rm=T),
      sd(xxx,na.rm=T)*2/sqrt(nsamplesMI))
})


t(minextpatient[,ordering])

ordering <- c(setdiff(order(minextpatient[1,]),which(colnames(minextpatient)=='all')),
              which(colnames(minextpatient)=='all'))
pdff('MI_ranks', paper='a4p')
tplot(x=minextpatient[1,ordering],y=1:ncol(minextpatient), type='p', pch=16,
      xlim=c(0,0.22), cex=1,
      xlab='mutual information/Sh',
      col=8, mar=c(NA,1,1,1),
 ##     xticks=seq(0.5,0.752),#cex.axis=1,
      yticks=NA,ylab=NA)
tplot(x=rbind(minextpatient[1,ordering]-minextpatient[2,ordering],
              minextpatient[1,ordering]+minextpatient[2,ordering]),
      y=rbind(1:ncol(minextpatient), 1:ncol(minextpatient)),
      col=8, lty=1, lwd=2,
      add=T)
text(x=minextpatient[1,ordering]+minextpatient[2,ordering]+0.002,
     y=1:ncol(minextpatient),
     colnames(minextpatient)[ordering], adj=c(0,0.5),
     cex=1.5)
dev.off()

## sequ <- 1:ncol(mism)
## mord <- order(mism[1,])
## accus2 <- accus[1,mord]
## errs <- accus[2,mord]
## dist <- max(errs)
## labs <- shortnames2[mord]
## pdff('plotAcc',paper='a4p')
## tplot(y=sequ[1:length(accus2)],x=accus2,pch=16,cex=1,col=8,type='p',
##       yticks=NA,ylab=NA,xlab='accuracy',
##       xlim=c(0.5,0.75))#max(accus2)+0.1))
## tplot(x=rbind(pmax(0,accus2-errs),accus2+errs),
## y=rbind(sequ,sequ),lty=1,lwd=2,col=8,
## add=T)
## text(accus2+dist,sequ,labs, adj=c(-0.05,0.5), cex=1.25)
## dev.off()


cbind(sort(accnextpatient[1,]))

tplot()



## choose subset of typical points
set.seed(101)
datapoints2 <- t(generateVariates(Ynames=c(predictors,predictands), X=NULL,
                                mcsamples=mcsamples, varinfo=varinfo,
                                    n=nrow(mcsamples))[,,1])


## conditional probability distributions of predictand=1 for subset
probAD2 <- samplesFDistribution(Y=predictandvalues['cAD',,drop=F],
                                X=datapoints2[,predictors,drop=F],
                                mcsamples=mcsamples, varinfo=varinfo,
                                jacobian=TRUE, fn=identity)

## saveRDS(probAD2,'probAD2.rds')
## saveRDS(datapoints2,'datapoints2.rds')

## probAD2bis <- samplesFDistribution2(Y=datapoints2[,predictors,drop=F],
##                                 X=NULL,
##                                 mcsamples=mcsamples, varinfo=varinfo,
##                                 jacobian=TRUE, fn=identity)

probAD2 <- probAD2[-which(is.na(probAD2), arr.ind = T)[,'row'],]


meanprobAD2 <- rowMeans(probAD2,na.rm=T)
OprobAD2 <- t(apply(probAD2,1,function(xxx)tquant(xxx,c(1,7)/8)))
CprobAD2 <- t(apply(probAD2,1,function(xxx)tquant(xxx,c(5,95)/100)))
accprobAD2 <- sapply(meanprobAD2,function(xxx)max(xxx,1-xxx))

ordering <- order(meanprobAD2)
pdff('plot1000patients')
tplot(y=meanprobAD2[ordering], ylim=0:1, col=2, lwd=3,
      xlab='sample patient #', ylab='probability of conversion to AD')
## plotquantiles(x=1:nrow(datapoints2),y=CprobAD2[ordering,], col=5,alpha=0.75)
plotquantiles(x=1:length(ordering),y=OprobAD2[ordering,], col=6,alpha=0.75)
dev.off()

## tplot(y=list(meanprobAD2[ordering],accprobAD2[ordering]), ylim=0:1)
## plotquantiles(x=1:nrow(datapoints2),y=cbind(O1probAD2,O7probAD2)[ordering,])

## tplot(y=accprobAD2[ordering], ylim=0:1)

## tplot(y=probAD2[ordering,1:100], ylim=0:1,
##       lty=1,lwd=1,col=5,alpha=0.9)
## tplot(y=meanprobAD2[ordering], col=1,lwd=3, add=T)

mean(pmax(meanprobAD2,1-meanprobAD2))

nsamplesperpoint <- 2
samplesresults <- sapply(meanprobAD2,function(xxx){
    sample(c(0,1), nsamplesperpoint, prob=c(1-xxx,xxx), replace=T)
})

samplesresultsnoise <- samplesresults+runif(n=length(samplesresults), -0.1,0.1)

tplot(y=t(samplesresultsnoise[,ordering]), type='p', pch=16, cex=0.5, col=1,
      ylim=c(-0.5,1.5), yticks=0:1,ylabels=c('sMCI','cAD'))

## selection of maximum exp. utility
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
pdff('plot1000patients')
tplot(x=-10,y=-10,xlim=c(0,length(meanprobAD2)),ylim=0:1,
            xlab='sample patient #', ylab='probability of conversion to AD')
plotquantiles(x=1:length(meanprobAD2),y=OprobAD2[ordering,], col=6,alpha=0.75)
tplot(y=meanprobAD2[ordering], ylim=0:1, col=2, lwd=3,
      xlab='sample patient #', ylab='probability of conversion to AD',
      add=T)
text(x=subsett, y=meanprobAD2[ordering][subsett],
     syms[treatmean[ordering][subsett]], cex=1.25, col=colalpha2hex(8,0.25))
## tplot(x=as.list(subsett),y=as.list(meanprobAD2[ordering][subsett]), type='p',
##       pch=syms[treatmean[ordering][subsett]],
##       col=1, cex=1.25,
##       add=T)
## plotquantiles(x=1:nrow(datapoints2),y=CprobAD2[ordering,], col=5,alpha=0.75)
dev.off()

rangestreat <- foreach(xxx=treat,.combine=cbind)%do%{
    tabulate(xxx, nbins=4)/length(xxx)
}

tplot(x=1:4,y=rangestreat[,1:1000],lty=1,lwd=0.5,alpha=0.75,col=7)

summary(t(rangestreat))

##        [,1]  [,2]     [,3]  [,4]
## 12.5% 0.220 0.026 0.258000 0.261
## 50%   0.247 0.037 0.357000 0.357
## 87.5% 0.269 0.062 0.465875 0.442


## > apply(rangestreat,1,function(xx){tquant(xx,c(1,4,7)/8)})
##           [,1]      [,2]     [,3]     [,4]
## 12.5% 0.227406 0.0214949 0.261877 0.249176
## 50%   0.248901 0.0293112 0.367855 0.351246
## 87.5% 0.265511 0.0532487 0.477956 0.447240

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


pdff('plot4096patients',paper='a4p')
tplot(y=-10,x=-10,ylim=c(0,length(meanprobAD2)),xlim=0:1,
            ylab='sample patient #', xlab='probability of conversion to AD')
#plotquantiles(y=1:length(meanprobAD2),x=OprobAD2[ordering,], col=6,alpha=0.75)
tplot(x=meanprobAD2[ordering], xlim=0:1, col=2, lwd=3,
      ylab='sample patient #', xlab='probability of conversion to AD',
      add=T)
text(y=subsett, x=meanprobAD2[ordering][subsett],
     syms[treatmean[ordering][subsett]], cex=1.25, col=colalpha2hex(8,0.25))
## tplot(x=as.list(subsett),y=as.list(meanprobAD2[ordering][subsett]), type='p',
##       pch=syms[treatmean[ordering][subsett]],
##       col=1, cex=1.25,
##       add=T)
## plotquantiles(x=1:nrow(datapoints2),y=CprobAD2[ordering,], col=5,alpha=0.75)
dev.off()




set.seed(500)
npsamples <- 1024
probsamples2 <- foreach(asample=t(mcsamples), .combine=cbind)%dorng%{
    asample <- t(asample)
    datasamples <- t(generateVariates(Ynames=c(predictors,predictands), X=NULL,
                                      mcsamples=asample, varinfo=varinfo,
                                      n=npsamples)[,,1])
    yvals <- datasamples[,predictands,drop=F]
    out <- c(samplesFDistribution2(Y=yvals,
                                X=datasamples[,predictors,drop=F],
                                mcsamples=asample, varinfo=varinfo,
                                jacobian=TRUE, fn=identity))
    cbind(out,
          c(samplesFDistribution2(Y=yvals,
                                X=NULL,
                                mcsamples=asample, varinfo=varinfo,
                                jacobian=TRUE, fn=identity)),
          (c(yvals)==1) * out + (c(yvals)==0) * (1-out)
          )
}
attr(probsamples, 'rng') <- NULL
attr(probsamples, 'doRNG_version') <- NULL
dim(probsamples) <- c(nrow(probsamples),3,ncol(probsamples)/3)
dimnames(probsamples) <- list(NULL,c('p(Y|X)','p(Y)','p(1|X)'),NULL)
probsamples <- aperm(probsamples)
## dim1: mcsamples, dim2:p/Yvalue, dim3: samples given lim-freq
##saveRDS(probsamples, '_probsamples.rds')







set.seed(500)
npsamples <- 4096
probsamples <- foreach(asample=t(mcsamples), .combine=cbind)%dorng%{
    asample <- t(asample)
    datasamples <- t(generateVariates(Ynames=c(predictors,predictands), X=NULL,
                                      mcsamples=asample, varinfo=varinfo,
                                      n=npsamples)[,,1])
    yvals <- datasamples[,predictands,drop=F]
    out <- c(samplesFDistribution2(Y=yvals,
                                X=datasamples[,predictors,drop=F],
                                mcsamples=asample, varinfo=varinfo,
                                jacobian=TRUE, fn=identity))
    cbind(out,
          c(samplesFDistribution2(Y=yvals,
                                X=NULL,
                                mcsamples=asample, varinfo=varinfo,
                                jacobian=TRUE, fn=identity)),
          (c(yvals)==1) * out + (c(yvals)==0) * (1-out)
          )
}
attr(probsamples, 'rng') <- NULL
attr(probsamples, 'doRNG_version') <- NULL
dim(probsamples) <- c(nrow(probsamples),3,ncol(probsamples)/3)
dimnames(probsamples) <- list(NULL,c('p(Y|X)','p(Y)','p(1|X)'),NULL)
probsamples <- aperm(probsamples)
## dim1: mcsamples, dim2:p/Yvalue, dim3: samples given lim-freq
##saveRDS(probsamples, '_probsamples.rds')

probsamples <- readRDS('_probsamples.rds')

allmutualinfodistr <- apply(probsamples,1,function(xxx){
    mean(log2(xxx['p(Y|X)',])) -
        mean(log2(xxx['p(Y)',]))
})

tquant(allmutualinfodistr, c(1,4,7)/8)
##    12.5%      50%    87.5% 
## 0.127781 0.158739 0.192135 

allcondentrdistr <- apply(probsamples,1,function(xxx){
    -mean(log2(xxx['p(Y|X)',]))
})

tquant(allcondentrdistr, c(1,4,7)/8)
##    12.5%      50%    87.5% 
## 0.803513 0.836117 0.867244 

allaccdistr <- apply(probsamples,1,function(xxx){
    mean((xxx['p(Y|X)',] > 0.5) + 0.5*(xxx['p(Y|X)',] == 0.5))
})

tquant(allaccdistr, c(1,4,7)/8)
##    12.5%      50%    87.5% 
## 0.666992 0.689941 0.712646 


alltreatdistr <- apply(probsamples,1,function(xxx){
    tabulate(choicefn(xxx['p(Y|X)',], um), nbins=4)/ncol(xxx)
})

apply(alltreatdistr,1, tquant, prob=c(1,4,7)/8)*100
##            [,1]      [,2]     [,3]     [,4]
## 12.5% 0.0319824 0.0129395 0.270752 0.420227
## 50%   0.0454102 0.0483398 0.387207 0.512207
## 87.5% 0.0742188 0.1061707 0.493652 0.589813

apply(alltreatdistr,1, tquant, prob=c(5,95)/100)*100
##        [,1]     [,2]    [,3]    [,4]
## 5%  2.51465  0.65918 22.0703 34.2493
## 95% 8.34961 13.06152 59.0649 63.0164



pgrid <- seq(0,1,by=0.02)
pgrid2 <- pgrid[-1]-diff(pgrid)/2
longrunprobs <- apply((probsamples),1,function(xxx){
thist(xxx['p(1|X)',], n=pgrid)$density*diff(pgrid)[1]
})

statlongrunprobs <- apply(longrunprobs,1,function(xxx){
c(tquant(xxx, c(8*2.5/100,1,2,4,6,7,8*97.5/100)/8), mean=mean(xxx))
})

pdff('distribution_prognostic_probs2',paper='a4p')
tplot(pgrid2,t(statlongrunprobs), col='white',
      xlab=expression(italic(P):~'probability of conversion'~'(2% bins)'),
      ylab=expression('fraction of population prognosed with'~italic(P)),
      ly=3,mar=c(NA,4.5,2,NA), cex.axis=1.25,
      ## yticks=seq(0,ceiling(max(statlongrunprobs)),by=1),
      xticks=seq(0,1,by=0.1),#ceiling(max(statlongrunprobs)),by=1),
      xlabels=paste0(seq(0,1,by=0.1)*100,'%')#ceiling(max(statlongrunprobs)),by=1),'%')
      )
plotquantiles(pgrid2,t(statlongrunprobs[c(1,7),]), col=5,alpha=0.75)
plotquantiles(pgrid2,t(statlongrunprobs[c(2,6),]), col=5,alpha=0.75)
plotquantiles(pgrid2,t(statlongrunprobs[c(3,5),]), col=5,alpha=0.75)
tplot(pgrid2,statlongrunprobs['50%',],add=T)
tplot(pgrid2,thist(origp,n=pgrid)$density*diff(pgrid)[1],lty=2,lwd=2,add=T)
##
tplot(pgrid2,statlongrunprobs['mean',],lty=3,lwd=2,add=T)
legend('topleft',legend=c('median','mean',
                           '50% uncertainty','75% uncertainty','95% uncertainty'),
       col=c(1,1,
             alpha2hex(5,c(0.25,0.5,0.75))),
       lwd=c(2,2,
             5,10,15),
       lty=c(1,2,
             1,1,1),
       bty='n'
       )
dev.off()

subsett <- round(seq(1,ncol(longrunprobs),length.out=5))
pdff('samples_prognostic_probs',paper='a4p')
tplot(pgrid2,longrunprobs[,subsett], lty=1, alpha=0.5,
      xlab=expression(italic(P):~'probability of conversion'~'(2% bins)'),
      ylab=expression('fraction of population prognosed with'~italic(P)),
      ly=3,mar=c(NA,4.5,2,NA), cex.axis=1.25,
      ## yticks=seq(0,ceiling(max(statlongrunprobs)),by=1),
      xticks=seq(0,1,by=0.1),#ceiling(max(statlongrunprobs)),by=1),
      xlabels=paste0(seq(0,1,by=0.1)*100,'%')#ceiling(max(statlongrunprobs)),by=1),'%')
      )
## plotquantiles(pgrid2,t(statlongrunprobs[c(1,5),]), col=5,alpha=0.75)
## plotquantiles(pgrid2,t(statlongrunprobs[c(2,4),]), col=5,alpha=0.75)
## tplot(pgrid2,statlongrunprobs['50%',],add=T)
#tplot(pgrid2,statlongrunprobs['mean',],lty=2,lwd=1,add=T)
legend('topright',legend=c('median',#'mean',
                           '50% uncertainty','75% uncertainty'),
       col=c(1,#1,
             alpha2hex(5,c(0.5,0.75))),
       lwd=c(2,#1,
             5,10),
       lty=c(1,#2,
             1,1),
       bty='n'
       )
dev.off()






pdff('single_predictors')
for(v in predictors){
    xgrid <- cbind(seq(varinfo[['plotmin']][v], varinfo[['plotmax']][v],
                       length.out=min(varinfo[['n']][v],256)))
    colnames(xgrid) <- v
    ppp <- samplesFDistribution(Y=predictandvalues['cAD',,drop=F],
                                X=xgrid,
                                mcsamples=asample, varinfo=varinfo,
                                jacobian=TRUE, fn=mean)[,1]
    tplot(xgrid,ppp,xlab=v, ylab='p of conversion',ylim=0:1)
}
dev.off()

tplot(x=pgrid2, longrunprobs[,round(seq(1,ncol(longrunprobs),length.out=100))], col=7, lty=1,alpha=0.75)

Olongrunprobs <- apply(probsamples,1,function(xxx){
c(tquant(xxx['p(1|X)',], c(1,4,7)/8), mean=mean(xxx))
})

thist(Olongrunprobs['mean',],plot=T)




set.seed(673)
basepointsa <- rbeta(n=10, 2,2)-0.5
basepointsb <- rbeta(n=10, 3,3)-0.5
basepointsc <- rbeta(n=10, 2,2)-0.5
basepointsd <- rbeta(n=10, 3,3)-0.5
##
impdatax <- cbind(
    c(basepointsa-1, basepointsb+0.25),
    c(basepointsa-1, basepointsb+0.25)
)
impdatay <- cbind(
    c(basepointsc-1, (basepointsd)+0.25),
    c(basepointsd+0.25, (basepointsc)-1)
)
##
pdff('example_importance_context3',asp=1)
tplot(x=impdatax,y=impdatay, type='p', pch=c(1,2),
      xlim=c(-2,0.75), ylim=c(-2,0.75),
      xticks=NA, yticks=NA,xlab=expression(italic(X)[1]), ylab=expression(italic(X)[2]),
      mar=c(NA,3.5,1,1), ly=2)
tplot(x=impdatax,y=t(t(impdatay*0)-c(1.95,2.05)), type='p',, pch=c(1,2),
      cex=1,
      add=T)
tplot(y=impdatay,x=t(t(impdatax*0)-c(1.95,2.05)), type='p',, pch=c(1,2),
      cex=1,
      add=T)
abline(h=-1.75, col=7)
abline(v=-1.75, col=7)
## text(-0.5,-2.25,'marginals',adj=c(0.5,0.5),xpd=NA,col=8)
## text(-2.25,0,'marginals',adj=c(1,0.5),xpd=NA,srt=90,col=8)
dev.off()




savedataOABC <- readRDS('datapatientsOABC.rds')

var <- 'LRHHC_n_long'
xgrid <- cbind(seq(varinfo$plotmin[var], varinfo$plotmax[var], by=max(varinfo$delta[var], (varinfo$plotmax[var]-varinfo$plotmin[var])/256)))
colnames(xgrid) <- var
yvals <- (samplesFDistribution(Y=xgrid, X=rbind(savedataOABC['curtis',-which(colnames(savedataOABC)==var)]),
                             mcsamples=mcsamples, varinfo=varinfo,
                             jacobian=TRUE))
yvals3 <- (samplesFDistribution(Y=xgrid, X=NULL,
                             mcsamples=mcsamples, varinfo=varinfo,
                             jacobian=TRUE))

##
samplesyvals <- c(generateVariates(Ynames='LRHHC_n_long', X=rbind(savedataOABC['curtis',-which(colnames(savedataOABC)==var)]),
                                mcsamples=mcsamples, varinfo=varinfo,
                                n=nrow(mcsamples)))
cove <- round(tquant(samplesyvals, c(2.5,50,97.5)/100)*1e3,1)

pdff(paste0('curtis_distr_',allshortnames[var]))
tplot(x=xgrid*1e3,
      y=cbind(rowMeans(yvals),rowMeans(yvals3)),
      ylim=c(0,NA),
      lty=c(1,5), col=c(3,7), alpha=c(0,0.125), lwd=c(4,4),
      xlab=paste0("Curtis's HV"),
      ylab=paste0('probability',if(varinfo$n[var]==Inf){' density'})
      )
legend('topright',
       legend=c('Curtis', 'Whole population'),
       ## legend=c(paste0('Curtis (median: ',cove['50%'],
       ##                     ',\n95% coverage: [',cove[1],', ',cove[3],'])'),
       ##                     'Whole population'),
       lty=c(1,2), col=c(3,7), lwd=c(4,3),
       bty='n')
dev.off()

yvals2 <- (samplesFDistribution(Y=predictandvalues['cAD',,drop=F], X=rbind(savedataOABC['olivia',]),
                             mcsamples=mcsamples, varinfo=varinfo,
                             jacobian=TRUE))

tquant(yvals2, c(2.5,50,97.5)/100)
##     2.5%      50%    97.5% 
## 0.193682 0.299063 0.426756 

var <- 'LRHHC_n_long'
xgrid <- cbind(seq(varinfo$plotmin[var], varinfo$plotmax[var], by=max(varinfo$delta[var], (varinfo$plotmax[var]-varinfo$plotmin[var])/256)))
colnames(xgrid) <- var
xvalues <- t(sapply(xgrid,function(xxx){temp <- savedataOABC['curtis',]; temp[var] <- xxx; temp}))
yvals1 <- (samplesFDistribution(Y=predictandvalues['cAD',,drop=F],
                                X=xvalues,
                             mcsamples=mcsamples, varinfo=varinfo,
                             jacobian=TRUE))
##

plotss <- round(seq(1,ncol(yvals1),length.out=100))
pdff(paste0('curtis_prob_conversion_',allshortnames[var]))
tplot(x=xgrid*(if(var=='LRHHC_n_long'){1e3}else{1}),
      y=yvals1[,plotss],
      ylim=c(0,1),
      lty=1, col=7, alpha=0.5, lwd=1,
      xlab=paste0("Curtis's unknown HV"),
      ylab=paste0("Curtis's probability of conversion to AD")
      )
tplot(x=xgrid*(if(var=='LRHHC_n_long'){1e3}else{1}),
      y=rowMeans(yvals1),
      lty=1, col=8, alpha=0.5, lwd=3,
      add=T)
## legend('topright', legend=c('will convert to AD', 'stable MCI'),
##        lty=c(2,1), col=c(2,5), lwd=3,
##        bty='n')
dev.off()
















temp <- thist(yvals2,n=seq(0,1,by=0.02))
tplot(x=temp$breaks,y=temp$density)









tplot(x=xgrid*1e3,
      y=cbind(rowMeans(yvals),rowMeans(yvals2),rowMeans(yvals3)),
      ylim=c(0,NA),
      lty=1,  alpha=0, lwd=4,
      xlab=paste0("Curtis's HV"),
      ylab=paste0('probability',if(varinfo$n[var]==Inf){' density'})
      )





predictand0 <- cbind(0)
predictand1 <- cbind(1)
colnames(predictand0) <- colnames(predictand1) <- predictands



test <- savedataOABC['curtis',-which(colnames(savedataOABC)=='LRHHC_n_long'),drop=F]
morecurtis <- c(generateVariates(Ynames='LRHHC_n_long', X=test,
                                mcsamples=mcsamples, varinfo=varinfo,
                                n=nrow(mcsamples)))

newcurtis <- cbind('LRHHC_n_long'=morecurtis)
for(i in 1:ncol(test)){
    newcurtis <- cbind(test[,i],newcurtis)
    colnames(newcurtis)[1] <- colnames(test)[i]
}

otherprobcurtis <- samplesFDistribution(Y=predictand1,
                                  X=newcurtis,
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=mean)

um <- matrix(c(1,0.9,0.8,0,
               0,0.3,0.5,1), 4,2)
rownames(um) <- c(1:4)

exputcurtis <- sapply(otherprobcurtis, function(xxx){
    um %*% c(1-xxx, xxx)
})

apply(exputcurtis,1,tquant, prob=c(0.5/100,5/100,1/8,7/8,95/100,99.5/100))*10
##          [,1]    [,2]    [,3]    [,4]
## 0.5%  2.97225 4.78335 5.89167 7.01188
## 5%    2.97226 4.78336 5.89168 7.01766
## 12.5% 2.97235 4.78341 5.89170 7.02096
## 87.5% 2.97904 4.78743 5.89371 7.02765
## 95%   2.98234 4.78940 5.89470 7.02774
## 99.5% 2.98812 4.79287 5.89643 7.02775













################### OLD #########################
tplot(x=1,y=1,type='p',pch=expression('beta'))
text(c(1,1.2),c(1,1.2),expression(alpha))

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




































pYXsamples <- readRDS()









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

testp <- t(generateVariates(Ynames=variatenames, X=NULL, mcsamples=mcsamples[1,,drop=F],
                          varinfo=varinfo, n=4096*4)[,,1])

v <- variatenames[13]
vgrid <- cbind(seq(varinfo[['plotmin']][v], varinfo[['plotmax']][v], length.out=min(256, varinfo[['n']][v])))
colnames(vgrid) <- v
testd <-samplesFDistribution(Y=vgrid, X=NULL,
                             mcsamples=mcsamples[1,,drop=F], varinfo=varinfo,
                             jacobian=TRUE, fn=mean)[,1]
##
testh <- thist(testp[,v],n=40)#,n=seq(varinfo[['plotmin']][v]-1e-4, varinfo[['plotmax']][v]+1e-4, length.out=min(128, varinfo[['n']][v])))
tplot(x=list(vgrid,testh$breaks), y=list(testd, testh$density/max(testh$density)*max(testd)),
      ylim=c(0,NA),xlab=v)


allshortnames <- readRDS('variatenames.rds')

var <- names(allshortnames)[allshortnames=='RAVLT-del']
xgrid <- cbind(seq(varinfo$plotmin[var], varinfo$plotmax[var], length.out=min(varinfo$n[var],256)))
colnames(xgrid) <- var
yvals1 <- (samplesFDistribution(Y=xgrid, X=predictandvalues['cAD',,drop=F],
                             mcsamples=mcsamples, varinfo=varinfo,
                             jacobian=TRUE))
yvals0 <- (samplesFDistribution(Y=xgrid, X=predictandvalues['sMCI',,drop=F],
                             mcsamples=mcsamples, varinfo=varinfo,
                             jacobian=TRUE))
##
plotss <- round(seq(1,ncol(yvals1),length.out=100))
pdff(paste0('population_distr_',allshortnames[var]))
tplot(x=xgrid,
      y=matrix(rbind(yvals1[,plotss],yvals0[,plotss]), nrow=length(xgrid)),
      ylim=c(0,NA),
      lty=1, col=c(2,5), alpha=0.9, lwd=1,
      xlab=allshortnames[var],
      ylab='population frequency'
      )
tplot(x=xgrid,
      y=cbind(rowMeans(yvals1),rowMeans(yvals0)),
      lty=c(2,1), col=c(6,1), alpha=0.5, lwd=3,
      add=T)
dev.off()


var <- names(allshortnames)[allshortnames=='HV']


#### plot marginal distributions

for(var in variatenames){
    print(var)
xgrid <- cbind(seq(varinfo$plotmin[var], varinfo$plotmax[var], by=max(varinfo$delta[var], (varinfo$plotmax[var]-varinfo$plotmin[var])/256)))
colnames(xgrid) <- var
yvals1 <- (samplesFDistribution(Y=xgrid, X=predictandvalues['cAD',,drop=F],
                             mcsamples=mcsamples, varinfo=varinfo,
                             jacobian=TRUE))
yvals0 <- (samplesFDistribution(Y=xgrid, X=predictandvalues['sMCI',,drop=F],
                             mcsamples=mcsamples, varinfo=varinfo,
                             jacobian=TRUE))
##
plotss <- round(seq(1,ncol(yvals1),length.out=100))
    breaks <- if(varinfo$n[var]==Inf){
                  seq(min(xgrid),max(xgrid),length.out=32)
              }else{c(xgrid[1]-varinfo$delta[var]/2, xgrid+varinfo$delta[var]/2)}
    histog0 <- thist(c(data.matrix(data0[Subgroup_num_==0,..var])), n=breaks)
    histog1 <- thist(c(data.matrix(data0[Subgroup_num_==1,..var])), n=breaks)
ymax <- max(histog0$density, histog1$density,yvals1[,plotss],yvals0[,plotss])
##
    pdff(paste0('population_distr_scat_',allshortnames[var]))
    tplot(x=histog0$mids*(if(var=='LRHHC_n_long'){1e3}else{1}), y=histog0$density, type='b', cex=0.75, col=3, alpha=0.65, lwd=1,
                ylim=c(0,ymax),
      xlab=allshortnames[var],
      ylab=paste0('population frequency',if(varinfo$n[var]==Inf){' density'})
)
    tplot(x=histog1$mids*(if(var=='LRHHC_n_long'){1e3}else{1}), y=histog1$density, type='b', pch=2, cex=0.75, col=4, alpha=0.75, lwd=1,
          add=T)
    tplot(x=xgrid*(if(var=='LRHHC_n_long'){1e3}else{1}),
      y=matrix(rbind(yvals1[,plotss],yvals0[,plotss]), nrow=length(xgrid)),
      lty=1, col=c(2,5), alpha=0.9, lwd=1,
      add=T)
tplot(x=xgrid*(if(var=='LRHHC_n_long'){1e3}else{1}),
      y=cbind(rowMeans(yvals1),rowMeans(yvals0)),
      lty=c(2,1), col=c(6,1), alpha=0.5, lwd=3,
      add=T)
legend('topright', legend=c('will convert to AD','data histogram','', 'stable MCI','data histogram'),
       lty=c(2,NA,NA,1,NA), pch=c(NA,2,NA,NA,1), col=c(2,4,NA,5,3), lwd=c(3,1,NA,3,1),
       bty='n')
dev.off()
}


for(var in variatenames){
    print(var)
xgrid <- cbind(seq(varinfo$plotmin[var], varinfo$plotmax[var], by=max(varinfo$delta[var], (varinfo$plotmax[var]-varinfo$plotmin[var])/256)))
colnames(xgrid) <- var
yvals1 <- (samplesFDistribution(X=xgrid, Y=predictandvalues['cAD',,drop=F],
                             mcsamples=mcsamples, varinfo=varinfo,
                             jacobian=TRUE))
##
plotss <- round(seq(1,ncol(yvals1),length.out=100))
pdff(paste0('prob_conversion_',allshortnames[var]))
tplot(x=xgrid*(if(var=='LRHHC_n_long'){1e3}else{1}),
      y=yvals1[,plotss],
      ylim=c(0,1),
      lty=1, col=7, alpha=0.5, lwd=1,
      xlab=allshortnames[var],
      ylab=paste0('population frequency of conversion to AD')
      )
tplot(x=xgrid*(if(var=='LRHHC_n_long'){1e3}else{1}),
      y=rowMeans(yvals1),
      lty=1, col=8, alpha=0.5, lwd=3,
      add=T)
## legend('topright', legend=c('will convert to AD', 'stable MCI'),
##        lty=c(2,1), col=c(2,5), lwd=3,
##        bty='n')
dev.off()
}


## varinfo$delta <- sapply(1:length(varinfo$n),function(i){
##         if(varinfo$n[i] < Inf){
##             c((varinfo$max[i]-varinfo$min[i])/(varinfo$n[i]-1))
##             }else{1/varinfo$n[i]}
##         })



countna <- apply(data0,1,function(xxx){sum(is.na(xxx))})
countna[countna > 0]
sum(countna > 0)
which(countna > 0)

data0[367,]

data0[which(countna > 0)]

