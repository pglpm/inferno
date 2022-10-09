## Author: PGL  Porta Mana
## Created: 2022-10-07T12:13:20+0200
## Last-Updated: 2022-10-09T14:14:48+0200
################
## Combine multiple Monte Carlo chains
################

outputdirectory <- '_test3NI33-V55-D40-K64-I1024'
totsamples <- 1024
variateinfofile <- 'metadata_noint33_noSW.csv' #***
datafile <- 'data_ep.csv' #***
mainvar <- 'group'
showdata <- 'histogram' # 'histogram' 'scatter' FALSE
plotmeans <- TRUE

rm(list=ls())
source('~/.Rprofile')
## load customized plot functions
if(!exists('tplot')){source('~/work/pglpm_plotfunctions.R')}
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
if(file.exists("/cluster/home/pglpm/R")){
    ncores <- availableCores()}else{
    ncores <- 4}
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

variateinfo <- fread(variateinfofile)
varNames <- variateinfo$variate
varTypes <- variateinfo$type
names(varTypes) <- varNames
realVars <- varNames[varTypes=='real']
integerVars <- varNames[varTypes=='integer']
categoryVars <- varNames[varTypes=='categorical']
binaryVars <- varNames[varTypes=='binary']
#varNames <- c(realVars, integerVars, categoryVars, binaryVars)
nrvars <- length(realVars)
nivars <- length(integerVars)
ncvars <- length(categoryVars)
nbvars <- length(binaryVars)
nvars <- length(varNames)
##
variateparameters <- data.matrix(read.csv(paste0(outputdirectory,'/variateparameters.csv'), row.names='variate'))
alldata <- fread(datafile, sep=',')
if(!all(varNames %in% names(alldata))){print('ERROR: variates missing from datafile')}
alldata <- alldata[, ..varNames]
## shuffle data
if(exists('shuffledata') && shuffledata){alldata <- alldata[sample(1:nrow(alldata))]}
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(alldata)}
alldata <- alldata[1:ndata]
##
source('functions_mcmc.R')
setwd(outputdirectory)
samplefiles <- list.files(pattern='^_mcsamples-')
setwd('..')
ncores <- max(as.integer(sub('.*-(\\d+).rds', '\\1', samplefiles)))
print(paste0('Found ',ncores,' chains'))
##
nstages <- max(as.integer(sub('.*--(\\d+)-\\d+.rds', '\\1', samplefiles)))
print(paste0('Found ',nstages,' stages'))
##
basename <- samplefiles[!(sub(paste0('^_[^-]+(-.+--)',nstages,'-\\d+\\.rds$'), '\\1', samplefiles)==samplefiles)][]
basename <- sub(paste0('^_[^-]+(-.+--)',nstages,'-\\d+\\.rds$'), '\\1', basename[1])
##
filesamples <- totsamples/ncores
print(paste0('need ESS > ',filesamples))
mcsamples <- foreach(i=1:ncores, .combine=rbind, .inorder=TRUE)%do%{
    traces <- readRDS(file=paste0(outputdirectory,'/_probtraces',basename,nstages,'-',i,'.rds'))
    minESS <- floor(min(LaplacesDemon::ESS(traces * (abs(traces) < Inf))))
    ## print(paste0('chain ',i,' ess=',minESS))
    traces <- readRDS(file=paste0(outputdirectory,'/_mcsamples',basename,nstages,'-',i,'.rds'))
    lsamples <- nrow(traces)
    if(minESS >= filesamples){
        print(paste0('chain ',i,': ESS = ',minESS))
        topick <- rev(seq(
        from=lsamples, length.out=filesamples, by=-ceiling(lsamples/minESS)
        ))
        }else{
        print(paste0('WARNING chain ',i,': insufficient sample size, ESS = ',minESS))
        topick <- rev(round(seq(
        from=lsamples, to=1, length.out=filesamples
        )))
        }
    traces[topick, ]
}
##
parmList <- mcsamples2parmlist(mcsamples, realVars, integerVars, categoryVars, binaryVars)
saveRDS(parmList,file=paste0(outputdirectory,'/_jointfrequencies-',totsamples,'.rds'))
##
## Traces
## Data (standardized for real variates)
dat <- list()
if(nrvars>0){ dat$Real=t((t(data.matrix(alldata[, ..realVars])) - variateparameters[realVars,'location'])/variateparameters[realVars,'scale'])}
## if(nivars>0){ dat$Integer=data.matrix(alldata[, ..integerVars])}
if(nivars>0){ dat$Integer=t((t(data.matrix(alldata[, ..integerVars])) - variateparameters[integerVars,'location']))}
## if(ncvars>0){ dat$Category=data.matrix(alldata[, ..categoryVars])}
if(ncvars>0){ dat$Category=t((t(data.matrix(alldata[, ..categoryVars])) - variateparameters[categoryVars,'location']))}
## if(nbvars>0){ dat$Binary=data.matrix(alldata[, ..binaryVars])}
if(nbvars>0){ dat$Binary=t((t(data.matrix(alldata[, ..binaryVars])) - variateparameters[binaryVars,'location']))}
##
##
print('Calculating traces...')
ll <- llSamples(dat, parmList)
condprobsd <- logsumsamplesF(Y=do.call(cbind,dat)[, mainvar, drop=F],
                             X=do.call(cbind,dat)[, setdiff(varNames, mainvar),
                                                  drop=F],
                             parmList=parmList, inorder=T)
condprobsi <- logsumsamplesF(Y=do.call(cbind,dat)[,
                                                  setdiff(varNames, mainvar), drop=F],
                             X=do.call(cbind,dat)[,
                                                  mainvar, drop=F],
                             parmList=parmList, inorder=T)
##
traces <- cbind(loglikelihood=ll, 'mean of direct logprobabilities'=condprobsd, 'mean of inverse logprobabilities'=condprobsi)*10/log(10)/ndata #medians, iqrs, Q1s, Q3s
badcols <- foreach(i=1:ncol(traces), .combine=c)%do%{if(all(is.na(traces[,i]))){i}else{NULL}}
if(!is.null(badcols)){traces <- traces[,-badcols]}
saveRDS(traces,file=paste0(outputdirectory,'/_jointprobtraces-',totsamples,'.rds'))
##
funMCSE <- function(x){LaplacesDemon::MCSE(x, method='batch.means')$se}
diagnESS <- LaplacesDemon::ESS(traces * (abs(traces) < Inf))
diagnIAT <- apply(traces, 2, function(x){LaplacesDemon::IAT(x[is.finite(x)])})
diagnBMK <- LaplacesDemon::BMK.Diagnostic(traces, batches=2)[,1]
diagnMCSE <- 100*apply(traces, 2, function(x){funMCSE(x)/sd(x)})
diagnStat <- apply(traces, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
diagnBurn <- apply(traces, 2, function(x){LaplacesDemon::burnin(matrix(x[1:(10*trunc(length(x)/10))], ncol=1))})
##
## tracegroups <- list(loglikelihood=1,
##                     'main given rest'=2,
##                     'rest given main'=3
##                     )
## grouplegends <- foreach(agroup=1:length(tracegroups))%do%{
##     c( paste0('-- STATS ', names(tracegroups)[agroup], ' --'),
##       paste0('min ESS = ', signif(min(diagnESS[tracegroups[[agroup]]]),6)),
##       paste0('max IAT = ', signif(max(diagnIAT[tracegroups[[agroup]]]),6)),
##       paste0('max BMK = ', signif(max(diagnBMK[tracegroups[[agroup]]]),6)),
##       paste0('max MCSE = ', signif(max(diagnMCSE[tracegroups[[agroup]]]),6)),
##       paste0('stationary: ', sum(diagnStat[tracegroups[[agroup]]]),'/',length(diagnStat[tracegroups[[agroup]]])),
##       paste0('burn: ', signif(max(diagnBurn[tracegroups[[agroup]]]),6))
##       )
##         }
colpalette <- c(7,2,1)
names(colpalette) <- colnames(traces)
##
## Plot various info and traces
print('Plotting traces and marginal samples')
family <- 'Palatino'
##
pdff(paste0(outputdirectory,'/jointmcsummary-',totsamples),'a4')
## Traces of likelihood and cond. probabilities
for(avar in colnames(traces)){
    tplot(y=traces[,avar], type='l', lty=1, col=colpalette[avar],
          main=paste0(avar,
                      '\nESS = ', signif(diagnESS[avar], 3),
                      ' | IAT = ', signif(diagnIAT[avar], 3),
                      ' | BMK = ', signif(diagnBMK[avar], 3),
                      ' | MCSE(6.27) = ', signif(diagnMCSE[avar], 3),
                      ' | stat: ', diagnStat[avar],
                      ' | burn: ', diagnBurn[avar]
                      ),
          ylab=paste0(avar,'/dHart'), xlab='sample', family=family
          )
}
## Samples of marginal frequency distributions
if(plotmeans){nfsamples <- totsamples}else{nfsamples <- 64}
for(avar in varNames){#print(avar)
    if(avar %in% realVars){
        rg <- signif((variateparameters[avar,c('thmin','thmax')] + 
                      7*variateparameters[avar,c('datamin','datamax')])/8, 2)
        if(!is.finite(rg[1])){rg[1] <- diff(variateparameters[avar,c('scale','datamin')])}
        if(!is.finite(rg[2])){rg[2] <- sum(variateparameters[avar,c('scale','datamax')])}
        Xgrid <- cbind(seq(rg[1], rg[2], length.out=256))
        ##histo <- thist(datum, n=32)#-exp(mean(log(c(round(sqrt(length(datum))), length(Xgrid))))))
    }else{
        rg <- round((variateparameters[avar,c('thmin','thmax')] + 
                     7*variateparameters[avar,c('datamin','datamax')])/8)
        Xgrid <- cbind(rg[1]:rg[2])
        ##histo <- thist(datum, n='i')
    }
    colnames(Xgrid) <- avar
    plotsamples <- samplesF(Y=Xgrid, parmList=parmList, nfsamples=min(nfsamples,nrow(mcsamples)), inorder=FALSE, rescale=variateparameters)
    ## ymax <- max(quant(apply(plotsamples,2,function(x){quant(x,99/100)}),99/100, na.rm=T), histo$density)
    ## tplot(x=histo$breaks, y=histo$density, col=yellow, lty=1, lwd=1, xlab=avar, ylab='probability density', ylim=c(0, ymax), family=family)
    ymax <- quant(apply(plotsamples,2,function(x){quant(x,99/100)}),99/100, na.rm=T)
    fiven <- variateparameters[avar,c('datamin','dataQ1','datamedian','dataQ2','datamax')]
    ##
        par(mfrow=c(1,1))
    tplot(x=Xgrid, y=plotsamples[,round(seq(1,nfsamples,length.out=64))], type='l', col=paste0(palette()[7], '44'), lty=1, lwd=2, xlab=paste0(avar,' (',variateinfo[variate==avar,type],')'), ylab=paste0('frequency',if(avar %in% realVars){' density'}), ylim=c(0, ymax), family=family)
    if(plotmeans){
        tplot(x=Xgrid, y=rowMeans(plotsamples), type='l', col=paste0(palette()[1], '88'), lty=1, lwd=3, add=T)
    }
    abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'),lwd=4)
    if(showdata=='histogram'){
        datum <- alldata[[avar]]
        histo <- thist(datum, n=(if(avar %in% realVars){NULL}else{'i'}))#-exp(mean(log(c(round(sqrt(length(datum))), length(Xgrid))))))
        tplot(x=histo$breaks, y=histo$density/max(histo$density)*max(rowMeans(plotsamples)), col=grey, alpha=0.75, border=NA, lty=1, lwd=1, family=family, ylim=c(0,NA), add=TRUE)
    }else if(showdata=='scatter'){
        datum <- alldata[[avar]]
        scatteraxis(side=1, n=NA, alpha='88', ext=8, x=datum+rnorm(length(datum),mean=0,sd=prod(variateparameters[avar,c('precision','scale')])/(if(avar %in% binaryVars){16}else{16})),col=yellow)
    }
}
dev.off()
print('Done')
