## Author: PGL  Porta Mana
## Created: 2022-10-07T12:13:20+0200
## Last-Updated: 2022-10-21T14:27:32+0200
################
## Combine multiple Monte Carlo chains
################

rm(list=ls())

outputdir <- '_testESS'
totsamples <- 1024*35
variateinfofile <- 'metadata_noint33_noSW.csv' #***
datafile <- 'data_ep.csv' #***
mainvar <- 'group'
showdata <- TRUE # 'histogram' 'scatter' FALSE
plotmeans <- TRUE

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
if(Sys.info()['nodename']=='luca-HP-Z2-G9'){
    ncores <- 20}else{
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

outputdir <- sub('(.+)/','\\1',outputdir)
variateinfo <- fread(variateinfofile)
if(!file.exists(paste0(variateinfofile))){ stop(paste0('ERROR: cannot find file ',variateinfofile)) }

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
variateparameters <- data.matrix(read.csv(paste0(outputdir,'/variateparameters.csv'), row.names='variate'))
alldata <- fread(datafile, sep=',')
if(!all(varNames %in% names(alldata))){print('ERROR: variates missing from datafile')}
alldata <- alldata[, ..varNames]
## shuffle data
if(exists('shuffledata') && shuffledata){alldata <- alldata[sample(1:nrow(alldata))]}
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(alldata)}
alldata <- alldata[1:ndata]
##
source('functions_mcmc.R')
setwd(outputdir)
samplefiles <- list.files(pattern='^_mcsamples-')
setwd('..')
#########################
## inverting nstages-nchains just for this script
#########################
nchains <- max(as.integer(sub('.*-(\\d+).rds', '\\1', samplefiles)))
print(paste0('Found ',nchains,' chains'))
##
nstages <- max(as.integer(sub('.*--(\\d+)-\\d+.rds', '\\1', samplefiles)))
print(paste0('Found ',nstages,' stages'))
##
basename <- samplefiles[!(sub(paste0('^_[^-]+(-.+--)',nstages,'-\\d+\\.rds$'), '\\1', samplefiles)==samplefiles)][]
basename <- sub(paste0('^_[^-]+(-.+--)',nstages,'-\\d+\\.rds$'), '\\1', basename[1])
##
filesamples <- totsamples/nchains
print(paste0('need ESS > ',filesamples))
mcsamples <- foreach(i=0:nstages, .combine=rbind, .inorder=TRUE)%do%{
    traces <- readRDS(file=paste0(outputdir,'/_probtraces',basename,i,'-',1,'.rds'))
    minESS <- floor(min(LaplacesDemon::ESS(traces * (abs(traces) < Inf))))
    ## print(paste0('chain ',i,' ess=',minESS))
    traces <- readRDS(file=paste0(outputdir,'/_mcsamples',basename,i,'-',1,'.rds'))
    lsamples <- nrow(traces)
    if(minESS >= filesamples){
        print(paste0('stage ',i,': ESS = ',minESS))
        topick <- rev(seq(
        from=lsamples, length.out=filesamples, by=-ceiling(lsamples/minESS)
        ))
        }else{
        print(paste0('WARNING stage ',i,': insufficient ESS = ',minESS))
        topick <- rev(round(seq(
        from=lsamples, to=1, length.out=filesamples
        )))
        }
    traces[topick, ]
}
##
parmList <- mcsamples2parmlist(mcsamples, realVars, integerVars, categoryVars, binaryVars)
saveRDS(parmList,file=paste0(outputdir,'/_jointfrequencies-',outputdir,'-',totsamples,'.rds'))
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
saveRDS(traces,file=paste0(outputdir,'/_jointprobtraces-',outputdir,'-',totsamples,'.rds'))
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
pdff(paste0(outputdir,'/jointmcsummary-',outputdir,'-',totsamples),'a4')
pdf1 <- dev.cur()
pdff(paste0(outputdir,'/marginals-',outputdir,'-',totsamples),'a4')
pdf2 <- dev.cur()
## Traces of likelihood and cond. probabilities
dev.set(pdf1)
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
par(mfrow=c(1,1))
## Samples of marginal frequency distributions
if(plotmeans){nfsamples <- totsamples}else{nfsamples <- 64}
subsample <- round(seq(1,nfsamples, length.out=63))
for(avar in varNames){#print(avar)
    if(avar %in% realVars){
        rg <- signif((variateparameters[avar,c('thmin','thmax')] + 
                      7*variateparameters[avar,c('datamin','datamax')])/8, 2)
        if(!is.finite(rg[1])){rg[1] <- diff(variateparameters[avar,c('scale','datamin')])}
        if(!is.finite(rg[2])){rg[2] <- sum(variateparameters[avar,c('scale','datamax')])}
        Xgrid <- cbind(seq(rg[1], rg[2], length.out=256))
    }else{
        rg <- round((variateparameters[avar,c('thmin','thmax')] + 
                     7*variateparameters[avar,c('datamin','datamax')])/8)
        Xgrid <- cbind(rg[1]:rg[2])
    }
    colnames(Xgrid) <- avar
    plotsamples <- samplesF(Y=Xgrid, parmList=parmList, nfsamples=min(nfsamples,nrow(mcsamples)), inorder=FALSE, transform=variateparameters)
    fiven <- variateparameters[avar,c('datamin','dataQ1','datamedian','dataQ2','datamax')]
    ##
    ymax <- quant(apply(plotsamples,2,function(x){quant(x,99/100)}),99/100, na.rm=T)
    dev.set(pdf1)
    tplot(x=Xgrid, y=plotsamples[,subsample], type='l', col=paste0(palette()[7], '44'), lty=1, lwd=2, xlab=paste0(avar,' (',variateinfo[variate==avar,type],')'), ylab=paste0('frequency',if(avar %in% realVars){' density'}), ylim=c(0, ymax), family=family)
    if(plotmeans){
        tplot(x=Xgrid, y=rowMeans(plotsamples), type='l', col=paste0(palette()[1], '88'), lty=1, lwd=3, add=T)
    }
    abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'),lwd=4)
    ##
    if((showdata=='histogram')|(showdata==TRUE & !(avar %in% realVars))){
        datum <- alldata[[avar]]
        histo <- thist(datum, n=(if(avar %in% realVars){min(max(10,sqrt(ndata)),100)}else{'i'}))#-exp(mean(log(c(round(sqrt(length(datum))), length(Xgrid))))))
        histomax <- (if(avar %in% realVars){max(rowMeans(plotsamples))/max(histo$density)}else{1L})
        tplot(x=histo$breaks, y=histo$density*histomax, col=grey, alpha=0.75, border=NA, lty=1, lwd=1, family=family, ylim=c(0,NA), add=TRUE)
    }else if((showdata=='scatter')|(showdata==TRUE & (avar %in% realVars))){
        datum <- alldata[[avar]]
        scatteraxis(side=1, n=NA, alpha='88', ext=8, x=datum+rnorm(length(datum),mean=0,sd=prod(variateparameters[avar,c('precision','scale')])/(if(avar %in% binaryVars){16}else{16})),col=yellow)
    }
    ##
    dev.set(pdf2)
    marguncertainty <- t(apply(plotsamples, 1, function(x){quant(x, c(1,31)/32)}))
    tplot(x=Xgrid, y=rowMeans(plotsamples), type='l', col=paste0(palette()[1]), lty=1, lwd=3, xlab=paste0(avar,' (',variateinfo[variate==avar,type],')'), ylab=paste0('frequency',if(avar %in% realVars){' density'}), ylim=c(0, ymax), family=family)
    ##  93.75% marginal credibility intervals
    plotquantiles(x=Xgrid, y=marguncertainty, col=4, alpha=0.75)
    abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'),lwd=4)
    ##
    if((showdata=='histogram')|(showdata==TRUE & !(avar %in% realVars))){
        tplot(x=histo$breaks, y=histo$density*histomax, col=grey, alpha=0.75, border=NA, lty=1, lwd=1, family=family, ylim=c(0,NA), add=TRUE)
    }else if((showdata=='scatter')|(showdata==TRUE & (avar %in% realVars))){
        scatteraxis(side=1, n=NA, alpha='88', ext=8, x=datum+rnorm(length(datum),mean=0,sd=prod(variateparameters[avar,c('precision','scale')])/(if(avar %in% binaryVars){16}else{16})),col=yellow)
    }
##
}
dev.off()
dev.off()
print('Done')


