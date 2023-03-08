## Author: PGL  Porta Mana
## Created: 2022-10-07T12:13:20+0200
## Last-Updated: 2022-12-20T20:39:35+0100
################
## Combine multiple Monte Carlo chains
################
if(!exists('tplot')){source('~/work/pglpm_plotfunctions.R')}

#rm(list=ls())

nsamplesMI <- 1024
outputdir <- 'examplenewpatients1'
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


## tiemax <- function(xxx){ sample(rep( which(xxx==max(xxx)), 2), 1, replace=F) }



#### CALCULATION OF MUTUAL INFO ####

predictandvalues <- cbind(seq(varinfo[['min']][predictands],varinfo[['max']][predictands],length.out=varinfo[['n']][predictands]))
colnames(predictandvalues) <- predictands

nsamplesMI <- 8
MIdata <- foreach(v=c('',predictors), .combine=cbind)%do%{
    print(v)
    predictors0 <- setdiff(predictors,v)
    xsamples2 <- t(generateVariates(Ynames=c(predictors0,predictands), X=NULL,
                                    mcsamples=mcsamples, varinfo=varinfo,
                                    n=nsamplesMI)[,,1])
    saveRDS(xsamples2,paste0('xsamples2-',v,'-',nsamplesMI,'.rds'))
    ##
    condprobsy <- samplesFDistribution(Y=xsamples2[,predictors0,drop=F],
                                       X=NULL,
                                       mcsamples=mcsamples, varinfo=varinfo,
                                       jacobian=FALSE, fn=mean)
    ##
    condprobsx <- samplesFDistribution(Y=xsamples2[,predictands,drop=F],
                                       X=NULL,
                                       mcsamples=mcsamples, varinfo=varinfo,
                                       jacobian=FALSE, fn=mean)
    ##
    condprobsj <- samplesFDistribution(Y=xsamples2,
                                       X=NULL,
                                       mcsamples=mcsamples, varinfo=varinfo,
                                       jacobian=FALSE, fn=mean)
    ##
    saveRDS(condprobsx,paste0('condprobsx-',v,'-',nsamplesMI,'.rds'))
    saveRDS(condprobsy,paste0('condprobsy-',v,'-',nsamplesMI,'.rds'))
    saveRDS(condprobsj,paste0('condprobsj-',v,'-',nsamplesMI,'.rds'))
    ##
    mutualinfos2 <- log2(condprobsj)-log2(condprobsx)-log2(condprobsy)
    saveRDS(mutualinfos2,paste0('mutualinfos2-',v,'-',nsamplesMI,'.rds'))
    mutualinfos2
}
colnames(MIdata) <- c('all',predictors)


        print(summary(mutualinfos2))
mutualinfo2 <- mean(mutualinfos2)
mutualinfosd2 <- sd(mutualinfos2)/sqrt(nsamplesMI)
names(mutualinfo2) <- names(mutualinfosd2) <- 'Sh'
c(mutualinfo2,mutualinfosd2)


#### CALCULATION OF MUTUAL INFO (CORRECTED BY BASE RATE) ####

priorp <- 0.5
priorpvect <- c(1-priorp, priorp)
predictandvalues <- cbind(seq(varinfo[['min']][predictands],varinfo[['max']][predictands],length.out=varinfo[['n']][predictands]))
colnames(predictandvalues) <- predictands

xsamples <- aperm(
    sapply(predictandvalues, function(yyy){
    generateVariates(Ynames=predictors, X=cbind('Subgroup_num_'=yyy),
                     mcsamples=mcsamples, varinfo=varinfo,
                     n=nsamplesMI)
    },simplify='array')[,,1,]
)
dimnames(xsamples)[[1]] <- predictandvalues
saveRDS(xsamples,paste0('xsamples',nsamplesMI,'.rds'))

condprobs <- sapply(predictandvalues, function(xxx){
    sapply(predictandvalues, function(yyy){
    samplesFDistribution(Y=xsamples[xxx+1,,],
                         X=cbind('Subgroup_num_'=yyy),
                         mcsamples=mcsamples, varinfo=varinfo,
                         jacobian=FALSE, fn=mean)[,1]
    }, simplify='array')
}, simplify='array')
saveRDS(condprobs,paste0('condprobs',nsamplesMI,'.rds'))

mutualinfos <- rowSums(sapply(predictandvalues,function(yyy){
    priorpvect[yyy+1] * (
        log2(condprobs[,yyy+1,yyy+1]) -
        log2(rowSums(sapply(predictandvalues,function(xxx){
            priorpvect[xxx+1] * condprobs[,xxx+1,yyy+1]
        })))
    )
}))
print(summary(mutualinfos))
mutualinfo <- mean(mutualinfos)
mutualinfosd <- sd(mutualinfos)/sqrt(nsamplesMI)
names(mutualinfo) <- names(mutualinfosd) <- 'Sh'
c(mutualinfo,mutualinfosd)
saveRDS(mutualinfo,paste0('mutualinfo',nsamplesMI,'.rds'))


stop('end')

pdff('check_samples')
for(v in predictors){
        vn <- min(varinfo[['n']][v],256)
        vmin <- varinfo[['plotmin']][v]
        vmax <- varinfo[['plotmax']][v]
        vd <- (vmax-vmin)/(vn-1)
    for(i in 0:1){
        histo <- thist(xsamples[,v,i+1],
                       n=seq(vmin-vd/2,vmax+vd/2,length.out=vn+1))
        xgrid <- cbind(histo$mids)
        colnames(xgrid) <- v
    distri <- samplesFDistribution(Y=xgrid,
                         X=cbind('Subgroup_num_'=i),
                         mcsamples=mcsamples, varinfo=varinfo,
                         jacobian=TRUE, fn=mean)
    tplot(x=list(xgrid,xgrid),y=list(histo$density,distri))
    }
}
dev.off()

test <-     samplesFDistribution(Y=xsamples[1:10,,0+1],
                         X=cbind('Subgroup_num_'=0),
                         mcsamples=mcsamples, varinfo=varinfo, fn=mean)


test2 <-     samplesFDistribution(Y=xsamples[1:10,,0+1],
                         X=cbind('Subgroup_num_'=0),
                         mcsamples=mcsamples, varinfo=varinfo)








                X=cbind('Subgroup_num_'=yyy,
                    if(length(givens) > 0){data.matrix(datapatients[,..givens])} ),





npatients <- 1024
set.seed(999)
datapatients <- t(generateVariates(n=npatients, Ynames=variatenames, X=NULL, mcsamples=mcsamples, varinfo=varinfo)[,,1])
rownames(datapatients) <- paste0('patient',1:npatients)

llpatients <- aperm(sapply(1:length(predictandvalues), function(v2){
    samplesFDistribution(Y=datapatients[,predictors,drop=F],X=predictandvalues[v2,,drop=F],mcsamples=mcsamples,varinfo=varinfo)
},simplify='array'), c(2,1,3))
dimnames(llpatients) <- list(NULL, rownames(datapatients), predictandvalues)

neworder <- order(colMeans(llpatients)[,2]/(colMeans(llpatients)[,2]+colMeans(llpatients)[,1]))
llpatients <- llpatients[,neworder,]

## tplot(y=colMeans(llpatients))

pset <- sum(data0[[predictands]]==1)/length(data0[[predictands]])
ptry <- c(0.15,pset,0.65)
subs <- round(seq(1,nrow(mcsamples),length.out=100))
cola <- c(5,7,2)
colb <- c(1,8,6)
pdff('dependence_on_baserate')
for(i in 1:length(ptry)){
    p <- ptry[i]
    tplot(y=t(p*(llpatients)[subs,,2]/(p*(llpatients)[subs,,2]+(1-p)*(llpatients)[subs,,1])),
          xlab='feature-value set ID', ylab='probability of conversion to AD',
      ylim=c(0,1),lty=1,lwd=1,alpha=7/8,col=cola[i],add=(i>1))
    tplot(y=(p*colMeans(llpatients)[,2]/(p*colMeans(llpatients)[,2]+(1-p)*colMeans(llpatients)[,1])),
          ylim=c(0,1),lty=1,lwd=3,col=colb[i], add=T)
}
legend('topleft',paste0('prior prob.: ',signif(ptry,2)),
       col=colb,lty=1,lwd=3,
       bty='n')
dev.off()

ppatients <- c(0.5,0.5,0.5)

testres <- sapply(1:npatients, function(v1){
    mean(llpatients[,v1,2])*ppatients[v1]/(
        mean(llpatients[,v1,2])*ppatients[v1] +
        mean(llpatients[,v1,1])*(1-ppatients[v1])
    )
})

testres <- apply(llpatients,2,function(xxx){
    colMeans(xxx)[]
})

data0 <- fread(paste0(origdir,datafile), sep=',')
if(!all(unlist(variate) %in% colnames(data0))){cat('\nERROR: variates missing from datafile')}
data0 <- data0[, unlist(variate), with=F]
## shuffle data
if(exists('shuffledata') && shuffledata){data0 <- data0[sample(1:nrow(data0))]}
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(data0)}
data0 <- data0[1:ndata]

proportion <- 1/4 # 0/1 in learning set
dataad0 <- which(data0$Subgroup_num_==0)
dataad1 <- which(data0$Subgroup_num_==1)
nlearn0 <- round(length(dataad0)*proportion)
nlearn1 <- round(length(dataad1)*(1-proportion))
print(c(nlearn0,nlearn1))
sellearn0 <- sample(dataad0,nlearn0,replace=F)
sellearn1 <- sample(dataad1,nlearn1,replace=F)
selact0 <- setdiff(dataad0,sellearn0)
selact1 <- setdiff(dataad1,sellearn1)
##
dataexlearn <- data0[c(sellearn0,sellearn1)]
dataexact <- data0[c(selact0,selact1)]
print(c(nrow(dataexlearn),nrow(dataexact)))
##
saveRDS(dataexlearn,'dataexample_learn.rds')
saveRDS(dataexact,'dataexample_act.rds')
##
fwrite(dataexlearn,'dataexample_learn.csv')
fwrite(dataexact,'dataexample_act.csv')



samplefiles <- list.files(pattern='^_mcsamples-.*-F\\.rds')
## print(samplefiles)
## setwd('..')
nchains <- length(samplefiles)
## nchains <- max(as.integer(sub('.*(\\d+)-F.rds', '\\1', samplefiles)))
print(paste0('Found ',nchains,' chains'))
##
## nstages <- max(as.integer(sub('.*--(\\d+)-\\d+.rds', '\\1', samplefiles)))
## print(paste0('Found ',nstages,' stages'))
##
## basename <- samplefiles[!(sub(paste0('^_[^-]+(-.+--)','\\d+','-F','\\.rds$'), '\\1', samplefiles)==samplefiles)][]
basename <- sub(paste0('^_[^-]+(-.+--)','\\d+','-F','\\.rds$'), '\\1', samplefiles[1])
##
##
##
print(paste0('Reading from files ',basename,'...'))
chainlist <- sort(as.integer(sub(paste0('^_.*--','(\\d+)','-F','\\.rds$'), '\\1', samplefiles)))
## print(chainlist)
remainingsamples <- totsamples
print(paste0('need ESS > ', totsamples/nchains))
traces <- mcsamples <- NULL
donechains <- 0
for(achain in chainlist){
    donechains <- donechains+1
    totake <- round(remainingsamples/(nchains-donechains+1))
    temptrace <- readRDS(file=paste0('_mctraces',basename,achain,'-F.rds'))
    minESS <- floor(min(LaplacesDemon::ESS(temptrace * (abs(temptrace) < Inf))))
    lsamples <- nrow(temptrace)
    if(minESS >= totake){
        cat(paste0('chain ',achain,': ESS = ',minESS))
    }else{
        cat(paste0('WARNING chain ',achain,': insufficient ESS = ',minESS))
    }
    cat(paste0(' (req. ',totake,')\n'))
        topick <- round(seq(from=1, to=lsamples, length.out=totake))
    ##
    traces <- rbind(traces, temptrace[topick,])
    mcsamples <- rbind(mcsamples, readRDS(file=paste0('_mcsamples',basename,achain,'-F.rds'))[topick,])
    remainingsamples <- remainingsamples - totake
}
##
## parmList <- mcsamples2parmlist(mcsamples, realVars, integerVars, categoryVars, binaryVars)
saveRDS(mcsamples,file=paste0('_jointmcsamples-',basename,'-',totsamples,'.rds'))
saveRDS(traces,file=paste0('_jointtraces-',basename,'-',totsamples,'.rds'))
traces2 <- traces[apply(traces,1,function(x){all(is.finite(x))}),]
flagll <- nrow(traces) != nrow(traces2)

##
## Traces
## Data (standardized for real variates)
## dat <- list()
## if(nrvars>0){ dat$Real=t((t(data.matrix(alldata[, ..realVars])) - variateparameters[realVars,'location'])/variateparameters[realVars,'scale'])}
## ## if(nivars>0){ dat$Integer=data.matrix(alldata[, ..integerVars])}
## if(nivars>0){ dat$Integer=t((t(data.matrix(alldata[, ..integerVars])) - variateparameters[integerVars,'location']))}
## ## if(ncvars>0){ dat$Category=data.matrix(alldata[, ..categoryVars])}
## if(ncvars>0){ dat$Category=t((t(data.matrix(alldata[, ..categoryVars])) - variateparameters[categoryVars,'location']))}
## ## if(nbvars>0){ dat$Binary=data.matrix(alldata[, ..binaryVars])}
## if(nbvars>0){ dat$Binary=t((t(data.matrix(alldata[, ..binaryVars])) - variateparameters[binaryVars,'location']))}
##
##
## print('Calculating traces...')
## ll <- llSamplesmc(dat, parmList)
## condprobsd <- logsumsamplesF(Y=do.call(cbind,dat)[, mainvar, drop=F],
##                              X=do.call(cbind,dat)[, setdiff(varNames, mainvar),
##                                                   drop=F],
##                              parmList=parmList, inorder=T)
## condprobsi <- logsumsamplesF(Y=do.call(cbind,dat)[,
##                                                   setdiff(varNames, mainvar), drop=F],
##                              X=do.call(cbind,dat)[,
##                                                   mainvar, drop=F],
##                              parmList=parmList, inorder=T)
##
## traces <- cbind(loglikelihood=ll, 'mean of direct logprobabilities'=condprobsd, 'mean of inverse logprobabilities'=condprobsi)*10/log(10)/ndata #medians, iqrs, Q1s, Q3s
## badcols <- foreach(i=1:ncol(traces), .combine=c)%do%{if(all(is.na(traces[,i]))){i}else{NULL}}
## if(!is.null(badcols)){traces <- traces[,-badcols]}
## saveRDS(traces,file=paste0(outputdir,'/_jointprobtraces-',outputdir,'-',totsamples,'.rds'))
##
funMCSE <- function(x){LaplacesDemon::MCSE(x, method='batch.means')$se}
diagnESS <- LaplacesDemon::ESS(traces2)
diagnIAT <- apply(traces2, 2, function(x){LaplacesDemon::IAT(x)})
diagnBMK <- LaplacesDemon::BMK.Diagnostic(traces2[1:(4*trunc(nrow(traces2)/4)),], batches=4)[,1]
diagnMCSE <- 100*apply(traces2, 2, function(x){funMCSE(x)/sd(x)})
diagnStat <- apply(traces2, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
diagnBurn <- apply(traces2, 2, function(x){LaplacesDemon::burnin(matrix(x[1:(10*trunc(length(x)/10))], ncol=1))})
diagnBurn2 <- proposeburnin(traces2, batches=10)
diagnThin <- proposethinning(traces2)
##
cat(paste0('\nESSs (',totsamples/nchains,'): ',paste0(round(diagnESS), collapse=', ')))
cat(paste0('\nIATs: ',paste0(round(diagnIAT), collapse=', ')))
cat(paste0('\nBMKs: ',paste0(round(diagnBMK,3), collapse=', ')))
cat(paste0('\nMCSEs: ',paste0(round(diagnMCSE,2), collapse=', ')))
cat(paste0('\nStationary: ',paste0(diagnStat, collapse=', ')))
cat(paste0('\nBurn-in I: ',paste0(diagnBurn, collapse=', ')))
cat(paste0('\nBurn-in II: ',diagnBurn2))
cat('\n')
## cat(paste0('\nProposed thinning: ',paste0(diagnThin, collapse=', ')))
##       )
##         }
colpalette <- c(7,2,1)
names(colpalette) <- colnames(traces)
##

## Plot various info and traces
print('Plotting traces and marginal samples')

##
graphics.off()
pdff(paste0('jointmcsummary',extratitle,'-',basename,'-',totsamples),'a4')
pdf1 <- dev.cur()
pdff(paste0('marginals',extratitle,'-',basename,'-',totsamples),'a4')
pdf2 <- dev.cur()
## Traces of likelihood and cond. probabilities
dev.set(pdf1)
for(avar in colnames(traces)){
    tplot(y=traces[,avar], type='l', lty=1, col=colpalette[avar],
          main=paste0(avar,
                      '\nESS = ', signif(diagnESS[avar], 3),
                      ' | IAT = ', signif(diagnIAT[avar], 3),
                      ' | BMK = ', signif(diagnBMK[avar], 3),
                      ' | MCSE = ', signif(diagnMCSE[avar], 3),
                      ' | stat: ', diagnStat[avar],
                      ' | burn I: ', diagnBurn[avar],
                      ' | burn II: ', diagnBurn2
                      ),
          ylab=paste0(avar,'/dHart'), xlab='sample', family=family
          )
}
par(mfrow=c(1,1))
## Samples of marginal frequency distributions
if(plotmeans){nfsamples <- totsamples}else{nfsamples <- 100}
subsample <- round(seq(1,nfsamples, length.out=100))
for(v in unlist(variate)){
    contvar <- varinfo[['type']][v] %in% c('R','O','D')
    rg <- c(varinfo[['plotmin']][v], varinfo[['plotmax']][v])
    if(contvar){
        Xgrid <- cbind(seq(rg[1], rg[2], length.out=256))
    }else{
        Xgrid <- seq(varinfo[['min']][v], varinfo[['max']][v], length.out=varinfo[['n']][v])
        Xgrid <- cbind(Xgrid[Xgrid >= rg[1] & Xgrid <= rg[2]])
    }
    colnames(Xgrid) <- v
    plotsamples <- samplesFDistribution(Y=Xgrid, X=NULL, mcsamples=mcsamples, varinfo=varinfo, subsamples=round(seq(1,nrow(mcsamples),length.out=nfsamples)), jacobian=TRUE)
    ymax <- tquant(apply(plotsamples[,subsample],2,function(x){tquant(x,31/32)}),31/32, na.rm=T)
    if((showdata=='histogram' || showdata==TRUE) && !contvar){
        datum <- data0[[v]]
        datum <- datum[!is.na(datum)]
        nh <- (varinfo[['max']][v]-varinfo[['min']][v])/(varinfo[['n']][v]-1)
        nh <- seq(varinfo[['min']][v]-nh/2, varinfo[['max']][v]+nh/2, length.out=varinfo[['n']][v]+1)
        histo <- thist(datum, n=nh)
        ymax <- max(ymax,histo$counts/sum(histo$counts))
    }
    ##
    dev.set(pdf1)
    if(!(varinfo[['type']][v] %in% c('O','D'))){
        ##
        tplot(x=Xgrid, y=plotsamples[,subsample], type='l', col=5, alpha=7/8, lty=1, lwd=2,
              xlab=paste0(v, (if(varinfo[['type']][v] %in% c('I','B','C')){' (discrete)'}else{' (continuous)'})),
              ylab=paste0('frequency', (if(varinfo[['type']][v] %in% c('R','O','D')){' density'}else{''})),
              ylim=c(0, ymax), family=family)
        ##
        if(plotmeans){
            tplot(x=Xgrid, y=rowMeans(plotsamples, na.rm=T), type='l', col=1, alpha=0.25, lty=1, lwd=4, add=T)
        }
    }else{ # plot of a continuous doubly-bounded variate
        interior <- which(Xgrid > varinfo[['min']][v] & Xgrid < varinfo[['max']][v])
        tplot(x=Xgrid[interior], y=plotsamples[interior,subsample], type='l', col=5, alpha=7/8, lty=1, lwd=2,
              xlab=paste0(v, ' (continuous with deltas)'),
              ylab=paste0('frequency (density)'),
              ylim=c(0, ymax), family=family)
        if(length(interior) < length(Xgrid)){
            tplot(x=Xgrid[-interior], y=plotsamples[-interior,subsample,drop=F]*ymax, type='p', pch=2, cex=2, col=5, alpha=7/8, lty=1, lwd=2, xlab=paste0(v), ylab=paste0('frequency'), ylim=c(0, ymax), family=family,add=T)
        }
        if(plotmeans){
            tplot(x=Xgrid[interior], y=rowMeans(plotsamples, na.rm=T)[interior], type='l', col=1, alpha=0.25, lty=1, lwd=3, add=T)
            if(length(interior) < length(Xgrid)){
                tplot(x=Xgrid[-interior], y=rowMeans(plotsamples, na.rm=T)[-interior]*ymax, type='p', pch=2, cex=2, col=1, alpha=0.25, lty=1, lwd=3, add=T)
            }
        }
    }
    ##
    if((showdata=='histogram')||(showdata==TRUE && !contvar)){
        datum <- data0[[v]]
        datum <- datum[!is.na(datum)]
        ## fiven <- varinfo[v,c('min','Q1','Q2','Q3','max')]
        if(!(varinfo[['type']][v] %in% c('O','D'))){
            if(contvar){
                nh <- max(10,round(length(datum)/64))
            }else{
                nh <- (varinfo[['max']][v]-varinfo[['min']][v])/(varinfo[['n']][v]-1)
                nh <- seq(varinfo[['min']][v]-nh/2, varinfo[['max']][v]+nh/2, length.out=varinfo[['n']][v]+1)
            }
            histo <- thist(datum, n=nh)
            if(contvar){
                histomax <- max(rowMeans(plotsamples))/max(histo$density)
                tplot(x=histo$mids, y=histo$density*histomax, col=yellow, lty=1, alpha=2/4, border=darkgrey, border.alpha=3/4, lwd=4, family=family, ylim=c(0,NA), add=TRUE)
            }else{
                tplot(x=histo$mids, y=histo$counts/sum(histo$counts), col=yellow, alpha=2/4, border=darkgrey, border.alpha=3/4, lty=1, lwd=4, family=family, ylim=c(0,NA), add=TRUE)
            }
        }else{ # histogram for censored variate
            interior <- which(datum > varinfo[['min']][v] & datum < varinfo[['max']][v])
            histo <- thist(datum[interior], n=max(10,round(length(interior)/64)))
            interiorgrid <- which(Xgrid > varinfo[['min']][v] & Xgrid < varinfo[['max']][v])
            histomax <- 1#max(rowMeans(plotsamples)[interiorgrid])/max(histo$density)
            tplot(x=histo$mids, y=histo$density*histomax, col=yellow, alpha=2/4, border=darkgrey, border.alpha=3/4, lty=1, lwd=4, family=family, ylim=c(0,NA), add=TRUE)
            ##
            pborder <- sum(datum <= varinfo[['min']][v])/length(datum)
            if(pborder > 0){
                tplot(x=varinfo[['min']][v], y=pborder*ymax, type='p', pch=0, cex=2, col=7, alpha=0, lty=1, lwd=5, family=family, ylim=c(0,NA), add=TRUE)
            }
            ##
            pborder <- sum(datum >= varinfo[['max']][v])/length(datum)
            if(pborder > 0){
                tplot(x=varinfo[['max']][v], y=pborder*ymax, type='p', pch=0, cex=2, col=7, alpha=0, lty=1, lwd=5, family=family, ylim=c(0,NA), add=TRUE)
            }
        }
    }else if((showdata=='scatter')|(showdata==TRUE & contvar)){
        datum <- data0[[v]]
        datum <- datum[!is.na(datum)]
        diffdatum <- c(apply(cbind(c(0,diff(datum)),c(diff(datum),0)),1,min))/2
        scatteraxis(side=1, n=NA, alpha='88', ext=5,
                    x=datum+runif(length(datum),
                                  min=-min(diff(sort(unique(datum))))/4,
                                  max=min(diff(sort(unique(datum))))/4),
                                  col=yellow)
    }
        fiven <- fivenum(datum)
        abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'),lwd=4)
    ##
    dev.set(pdf2)
    marguncertainty <- t(apply(plotsamples, 1, function(x){tquant(x, quantilestoshow)}))
    tplot(x=Xgrid, y=rowMeans(plotsamples), type='l', col=1, lty=1, lwd=3, xlab=paste0(v), ylab=paste0('frequency', (if(varinfo[['type']][v] %in% c('R','O','D')){' density'}else{''})),
          ylim=c(0, ymax), family=family)
    ##  93.75% marginal credibility intervals
    plotquantiles(x=Xgrid, y=marguncertainty, col=5, alpha=0.75)
    abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'),lwd=4)
    ##
    if((showdata=='histogram')||(showdata==TRUE && !contvar)){
        if(contvar){
            tplot(x=histo$mids, y=histo$density*histomax, col=yellow, alpha=0.5, border=darkgrey, border.alpha=3/4, lty=1, lwd=4, family=family, ylim=c(0,NA), add=TRUE)
        }else{
            tplot(x=histo$mids, y=histo$counts/sum(histo$counts), col=yellow, alpha=0.5, border=darkgrey, border.alpha=3/4, lty=1, lwd=4, family=family, ylim=c(0,NA), add=TRUE)
}
    }else if((showdata=='scatter')|(showdata==TRUE & contvar)){
        scatteraxis(side=1, n=NA, alpha='88', ext=8,
                    x=datum+runif(length(datum),
                                  min=-min(diff(sort(unique(datum))))/4,
                                  max=min(diff(sort(unique(datum))))/4),
                                  col=yellow)
    }
##
}
dev.off(pdf1)
dev.off(pdf2)
cat('\nDone\n\n')
setwd(origdir)
