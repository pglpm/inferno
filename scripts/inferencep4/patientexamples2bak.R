## Author: PGL  Porta Mana
## Created: 2022-10-07T12:13:20+0200
## Last-Updated: 2022-12-29T00:59:39+0100
################
## Combine multiple Monte Carlo chains
################
if(!exists('tplot')){source('~/work/pglpm_plotfunctions.R')}

#rm(list=ls())

outputdir <- 'patientexamples2'
extratitle <- 'patientexamples'
## totsamples <- 4096L
datafile <- 'ingrid_data_nonangds6.csv'
predictorfile <- 'predictors.csv'
predictandfile <- NULL # 'predictors.csv'
mainfilelocation <- 'inferencep4/'
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
                      list.files(path=paste0(origdir,mainfilelocation), pattern='^varinfo.*\\.rds'))
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


predictand0 <- cbind(0)
predictand1 <- cbind(1)
colnames(predictand0) <- colnames(predictand1) <- predictands

givens <- c('Apoe4_', 'Gender_num_', 'AGE', 'LRHHC_n_long')
nogivens <- setdiff(variatenames, c(givens,predictands))
dataprior1 <- sum(data0[[predictands]]==1)/ndata


olivia <- ariel <- bianca <- readRDS('OCBdata.rds')
curtis <- readRDS('ABdata.rds')

Oll <- sapply(list(predictand0, predictand1), function(xxx){
    samplesFDistribution(Y=olivia[,nogivens,drop=F],
                                   X=cbind(olivia[,givens,drop=F], xxx),
                                       mcsamples=mcsamples, varinfo=varinfo,
                         jacobian=FALSE)
}, simplify='array')
colnames(Oll) <- c(0,1)
##
All <- sapply(list(predictand0, predictand1), function(xxx){
    samplesFDistribution(Y=ariel[,nogivens,drop=F],
                                   X=cbind(ariel[,givens,drop=F], xxx),
                                       mcsamples=mcsamples, varinfo=varinfo,
                         jacobian=FALSE)
}, simplify='array')
colnames(All) <- c(0,1)
##
Bll <- sapply(list(predictand0, predictand1), function(xxx){
    samplesFDistribution(Y=bianca[,nogivens,drop=F],
                                   X=cbind(bianca[,givens,drop=F], xxx),
                                       mcsamples=mcsamples, varinfo=varinfo,
                         jacobian=FALSE)
}, simplify='array')
colnames(All) <- c(0,1)
##
Cll <- sapply(list(predictand0, predictand1), function(xxx){
    samplesFDistribution(Y=curtis[,nogivens,drop=F],
                                   X=cbind(curtis[,givens,drop=F], xxx),
                                       mcsamples=mcsamples, varinfo=varinfo,
                         jacobian=FALSE)
}, simplify='array')
colnames(All) <- c(0,1)

Oprob1 <- Oll[,'1']*dataprior1/(Oll[,'1']*dataprior1 + Oll[,'0']*(1-dataprior1))
Oprob1m <- mean(Oll[,'1'])*dataprior1/(mean(Oll[,'1'])*dataprior1 + mean(Oll[,'0'])*(1-dataprior1))






##
set.seed(101)
## seed 101:
## OC patient: 182
## A patient: 170
## prior B: 0.3
ntrypatients <- 512
ngender <- 'Gender_num_'
trypatientsOC <- t(generateVariates(Ynames=setdiff(variatenames,c(ngender,predictands)), X=cbind('Gender_num_'=1),
                                mcsamples=mcsamples, varinfo=varinfo,
                                n=ntrypatients)[,,1])
trypatientsOC <- cbind(trypatientsOC, 'Gender_num_'=1)
##
trypatientsA <- t(generateVariates(Ynames=setdiff(variatenames,c(ngender,predictands)), X=cbind('Gender_num_'=0),
                                mcsamples=mcsamples, varinfo=varinfo,
                                n=ntrypatients)[,,1])
trypatientsA <- cbind(trypatientsA, 'Gender_num_'=0)


predictand0 <- cbind(0)
predictand1 <- cbind(1)
colnames(predictand0) <- colnames(predictand1) <- predictands
##
givens <- c('Apoe4_', 'Gender_num_', 'AGE', 'LRHHC_n_long')
nogivens <- setdiff(variatenames, c(givens,predictands))
dataprior1 <- sum(data0[[predictands]]==1)/ndata

llpatientsOC0 <- samplesFDistribution(Y=trypatientsOC[,nogivens,drop=F],
                                      X=cbind(trypatientsOC[,givens,drop=F],
                                              'Subgroup_num_'=0),
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=identity)
##
llpatientsOC1 <- samplesFDistribution(Y=trypatientsOC[,nogivens,drop=F],
                                      X=cbind(trypatientsOC[,givens,drop=F],
                                              'Subgroup_num_'=1),
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=identity)
##
llpatientsOCd <- samplesFDistribution(Y=predictand1,
                                  X=trypatientsOC[,predictors,drop=F],
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=identity)

ipatientsOC0 <- samplesFDistribution(Y=trypatientsOC[,predictors,drop=F],
                                      X=predictand0,
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=identity)
##
ipatientsOC1 <- samplesFDistribution(Y=trypatientsOC[,predictors,drop=F],
                                      X=predictand1,
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=identity)

##
##
##
## ratepatientsOC0 <- samplesFDistribution(Y=predictand0,
##                                       X=trypatientsOC[,givens,drop=F],
##                                   mcsamples=mcsamples, varinfo=varinfo,
##                                   jacobian=TRUE, fn=identity)
##
ratepatientsOC1 <- samplesFDistribution(Y=predictand1,
                                      X=trypatientsOC[,givens,drop=F],
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=identity)
##


## find Olivia
probs <- rowMeans(llpatientsOCd)
orderO <- order(abs(probs-0.3))
## 85
##                         [,1]
## Apoe4_             0.00000000
## ANARTERR_neuro    18.00000000
## AVDEL30MIN_neuro   5.00000000
## AVDELTOT_neuro    10.00000000
## CATANIMSC_neuro   21.00000000
## GDTOTAL_gds        3.00000000
## RAVLT_immediate   36.00000000
## TRAASCOR_neuro    21.00000000
## TRABSCOR_neuro   114.00000000
## AGE               75.38823609
## LRHHC_n_long       0.00425972
## Gender_num_        1.00000000
##
## posterior  0.301679


## find Ariel
ool <- 1
xgrid <- seq(0,1,length.out=256)
tplot(x=xgrid,
      y=sapply(ool,function(oo)mean(ipatientsOC1[orderO[oo],])*xgrid/(mean(ipatientsOC1[orderO[oo],])*xgrid + mean(ipatientsOC0[orderO[oo],])*(1-xgrid)))
      )
tplot(x=0:1,y=0:1,col=7,alpha=0.5,add=T)
abline(h=c(0.3,0.5,0.7))

oo <- 1
probs[orderO[oo]]
xgrid <- 0.65
mean(ipatientsOC1[orderO[oo],])*xgrid/(mean(ipatientsOC1[orderO[oo],])*xgrid + mean(ipatientsOC0[orderO[oo],])*(1-xgrid))
## posterior 0.47294


## find Curtis

##
llpatientsAd <- samplesFDistribution(Y=predictand1,
                                  X=trypatientsA[,predictors,drop=F],
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=identity)

probs <- rowMeans(llpatientsAd)
orderC <- order(abs(probs-0.8))
probs[orderC[1]]
## 98
## Apoe4_             1.00000000
## ANARTERR_neuro    24.00000000
## AVDEL30MIN_neuro   0.00000000
## AVDELTOT_neuro     4.00000000
## CATANIMSC_neuro   10.00000000
## GDTOTAL_gds        3.00000000
## RAVLT_immediate   27.00000000
## TRAASCOR_neuro    24.00000000
## TRABSCOR_neuro   145.00000000
## AGE               65.41445772
## LRHHC_n_long       0.00368181
## Gender_num_        0.00000000
## posterior 0.702336

choicefn <- function(pv,umx){
    sapply(pv, function(p){
        as.numeric(rownames(umx)[which.max(c(umx %*% c(1-p,p)))])
        })
}
##
aa <- function(a){1-a/4}
bb <- function(a){3/5+3*a/4}
b <- function(a){3/5-a/2}
a <- 5/15
matrix(c(1,aa(a),bb(a),0,
               0,a,b(a),1), 4,2)

um <- matrix(c(1,0.95,0.85,0,
               0,0.3,0.45,1), 4,2)
rownames(um) <- c(1:4)
um2 <- matrix(c(1,0.95,0.7,0,
               0,0.3,0.45,1), 4,2)
## um2 <- matrix(c(1,0.7,0.5,0,
##                 0,0.5,0.8,1),4,2)
rownames(um2) <- c(1:4)
##
pseq <- seq(0,1,length.out=256)
pdff('decisionthresholds3')
tplot(x=list(pseq,pseq), y=list(choicefn(pseq,um),choicefn(pseq,um2)),
      xlab='posterior probability')
abline(v=c(0.301679,0.47294,0.702336),col=3,lwd=2)
dev.off()







## for(i in 1:10){
i <- 5
##takewith1 <- (trypatientsOC[,predictands]==1)
orderOC <- order(abs(rowMeans(llpatientsOCd)-0.6))# + takewith1*1e12)
closeOC <- orderOC[i]
trypatientsOC[closeOC,]
## rowMeans(llpatientsOCd)[closeOC]
##
ll1 <- rowMeans(llpatientsOC1)[closeOC]
ll0 <- rowMeans(llpatientsOC0)[closeOC]
xgrid <- seq(0,1,length.out=256)
tplot(x=list(xgrid,xgrid),
      y=list(choicefn(xgrid,um),
             choicefn(ll1*xgrid/(ll1*xgrid+ll0*(1-xgrid)),um)),
      )
##
xgrid <- 0.3
postp <- ll1*xgrid/(ll1*xgrid+ll0*(1-xgrid))
    if(length(unique(c(
        choicefn(postp,um),choicefn(rowMeans(llpatientsOCd)[closeOC],um),choicefn(rowMeans(llpatientsOCd)[closeOC],um2)
    ))) >2 ){
        print(i)
        print(rowMeans(llpatientsOCd)[closeOC])
        print(postp)
print('B')
print(choicefn(postp,um))
print('O')
print(choicefn(rowMeans(llpatientsOCd)[closeOC],um))
print('C')
print(choicefn(rowMeans(llpatientsOCd)[closeOC],um2))
}
## }

llpatientsA0 <- samplesFDistribution(Y=trypatientsA[,nogivens,drop=F],
                                      X=cbind(trypatientsA[,givens,drop=F],
                                              'Subgroup_num_'=0),
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=identity)
##
llpatientsA1 <- samplesFDistribution(Y=trypatientsA[,nogivens,drop=F],
                                      X=cbind(trypatientsA[,givens,drop=F],
                                              'Subgroup_num_'=1),
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=identity)
##
llpatientsAd <- samplesFDistribution(Y=predictand1,
                                  X=trypatientsA[,predictors,drop=F],
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=identity)

orderA <- order(abs(rowMeans(llpatientsAd)-0.8))# + takewith1*1e12)
closeA <- orderA[1]
trypatientsA[closeA,]
postp <- rowMeans(llpatientsAd)[closeA]
postp
print(choicefn(postp,um))

saveRDS(trypatientsOC[closeOC,,drop=F],'OCBdata.rds')
saveRDS(llpatientsOC1[closeOC,],'OCBlikelihood1.rds')
saveRDS(llpatientsOC0[closeOC,],'OCBlikelihood0.rds')
saveRDS(llpatientsOCd[closeOC,],'OCBprobability1.rds')
##
saveRDS(trypatientsA[closeA,,drop=F],'ABdata.rds')
saveRDS(llpatientsA1[closeA,],'ABlikelihood1.rds')
saveRDS(llpatientsA0[closeA,],'ABlikelihood0.rds')
saveRDS(llpatientsAd[closeA,],'ABprobability1.rds')



#### Consistency checks for Bayes's theorem

llpatientsOCG0 <- samplesFDistribution(Y=predictand0,
                                  X=trypatientsOC[,givens,drop=F],
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=identity)
##
llpatientsOCG1 <- samplesFDistribution(Y=predictand1,
                                  X=trypatientsOC[,givens,drop=F],
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=identity)


llpatientsAG0 <- samplesFDistribution(Y=predictand0,
                                  X=trypatientsA[,givens,drop=F],
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=identity)
##
llpatientsAG1 <- samplesFDistribution(Y=predictand1,
                                  X=trypatientsA[,givens,drop=F],
                                  mcsamples=mcsamples, varinfo=varinfo,
                                  jacobian=TRUE, fn=identity)

checkposts <- (rowMeans(llpatientsOC1)*rowMeans(llpatientsOCG1)/(rowMeans(llpatientsOC1)*rowMeans(llpatientsOCG1)+rowMeans(llpatientsOC0)*(rowMeans(llpatientsOCG0))))/rowMeans(llpatientsOCd) - 1
summary(checkposts)
checkposts[closeOC]*100
checkposts <- rowMeans(llpatientsOC1*llpatientsOCG1/(llpatientsOC1*llpatientsOCG1+llpatientsOC0*llpatientsOCG0))/rowMeans(llpatientsOCd) - 1
summary(checkposts)
checkposts[closeOC]*100


checkposts <- (rowMeans(llpatientsA1)*rowMeans(llpatientsAG1)/(rowMeans(llpatientsA1)*rowMeans(llpatientsAG1)+rowMeans(llpatientsA0)*(rowMeans(llpatientsAG0))))/rowMeans(llpatientsAd) - 1
summary(checkposts)
checkposts[closeA]*100
checkposts <- rowMeans(llpatientsA1*llpatientsAG1/(llpatientsA1*llpatientsAG1+llpatientsA0*llpatientsAG0))/rowMeans(llpatientsAd) - 1
summary(checkposts)
checkposts[closeA]*100

## Find prior for B: decision threshold max 30%


um <- matrix(c(1,0.75,0, 0,0.75,1),3,2)
rownames(um) <- c(1:3)
um2 <- matrix(c(1,0.75,0,0.5, 0,0.5,1,0.75),4,2)
rownames(um2) <- c(1:4)
testres <- foreach(pat=1:ntrypatients, .combine=rbind)%do%{
    logl0 <- mean(llpatients0[pat,])
    logl1 <- mean(llpatients1[pat,])
    prior1 <- 0.15
    postp0 <- logl1*dataprior/(logl1*dataprior + logl0*(1-dataprior))
    postp1 <- logl1*prior1/(logl1*prior1 + logl0*(1-prior1))
    dec0 <- choicefn(postp0,um)
    dec1 <- choicefn(postp1,um)
    dec0b <- choicefn(postp0,um2)
    dec1b <- choicefn(postp1,um2)
    c(dec0,dec1,dec0b)
}

totry <- which(apply(testres,1,function(xxx)length(unique(xxx))==3))


##
choicefn <- function(pv,umx){
    sapply(pv, function(p){
        as.numeric(rownames(umx)[which.max(c(umx %*% c(1-p,p)))])
        })
}
##
##um <- matrix(c(1,0.6,0, 0,0.7,1),3,2)
um <- matrix(c(1,0.7,0.5,0,
               0,0.5,0.6,1),4,2)
rownames(um) <- c(1:4)
um2 <- matrix(c(1,0.7,0.5,0,
                0,0.5,0.8,1),4,2)
rownames(um2) <- c(1:4)
##
pseq <- seq(0,1,length.out=256)
pdff('decisionthresholds2')
tplot(x=list(pseq,pseq), y=list(choicefn(pseq,um),choicefn(pseq,um2)),
      xlab='posterior probability')
dev.off()

aa <- 0.8
um <- matrix(c(1,aa,7/3*(1-aa),0,
               0,7/3*(1-aa),aa,1),4,2)
rownames(um) <- c(1:4)
um
##
aa <- 0.65
um2 <- matrix(c(1,aa,16/15*(7-10*aa),0,
               0,7/3*(1-aa),8*aa/3-13/15,1),4,2)
## um2 <- matrix(c(1,aa,21/5-6*aa,0,
##                0,7/3*(1-aa),2/15*(5*aa+4),1),4,2)
rownames(um2) <- c(1:4)
um2
##
pseq <- seq(0,1,length.out=256)
pdff('decisionthresholds2')
tplot(x=list(pseq,pseq), y=list(choicefn(pseq,um),choicefn(pseq,um2)),
      xlab='posterior probability')
dev.off()


###
um <- matrix(c(1,0.8,0.5,0,
               0,0.5,0.8,1),4,2)
rownames(um) <- c(1:4)
um
um2 <- matrix(c(1,0.7,0.4,0,
                0,0.8,0.9,1),4,2)
## um2 <- matrix(c(1,0.65,0.55,0,
##                 0,0.85,0.85,1),4,2)
rownames(um2) <- c(1:4)
um2
##
pseq <- seq(0,1,length.out=256)
pdff('decisionthresholds2')
tplot(x=list(pseq,pseq), y=list(choicefn(pseq,um),choicefn(pseq,um2)),
      xlab='posterior probability')
dev.off()


#####
nn <- 1e6
aa <- runif(nn)
bb <- runif(nn)
cc <- runif(nn)
dd <- runif(nn)
which(8 > 8*aa+2*bb & 8 > 8*cc+2*dd &
      4*cc+6*dd >)



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


decreases <- (colMeans(MIdata)/mean(MIdata[,1])-1)*100
variances <- (
    abs(apply(MIdata,2,sd)/mean(MIdata[,1]) -
    colMeans(MIdata)*sd(MIdata[,1])/mean(MIdata[,1])^2)
        ) *100/sqrt(nsamplesMI)
sorto <- order(decreases)
signif(cbind(
    decreases[sorto],
    variances[sorto]
),4)
cbind(
    round(decreases[sorto],signif(-log10(variances[sorto]*5),1)),
    signif(variances[sorto]*5,1)
)



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


ysamples <- t(generateVariates(Ynames=sample(predictors), X=cbind('Subgroup_num_'=0),
                                mcsamples=mcsamples, varinfo=varinfo,
                               n=4096*4)[,,1])
##
pdff('testgeneratingfunction_conditional0')
for(v in predictors){
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
                                    X=cbind('Subgroup_num_'=0),
                                    mcsamples=mcsamples, varinfo=varinfo,
                                    jacobian=TRUE, fn=mean)
    histo <- thist(ysamples[,v],n=nh)
    tplot(list(histo$mids, xgrid),list(histo$density,testpdf[,1]),ylim=c(0,NA),xlab=v)
}
dev.off()

ysamples <- t(generateVariates(Ynames=sample(predictors), X=cbind('Subgroup_num_'=1),
                                mcsamples=mcsamples, varinfo=varinfo,
                               n=4096*4)[,,1])
##
pdff('testgeneratingfunction_conditional1')
for(v in predictors){
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
                                    X=cbind('Subgroup_num_'=1),
                                    mcsamples=mcsamples, varinfo=varinfo,
                                    jacobian=TRUE, fn=mean)
    histo <- thist(ysamples[,v],n=nh)
    tplot(list(histo$mids, xgrid),list(histo$density,testpdf[,1]),ylim=c(0,NA),xlab=v)
}
dev.off()


separate <- function(x){
    dx <- min(diff(sort(unique(x))))
    cts <- table(x)
}
