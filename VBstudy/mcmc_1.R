## Author: PGL  Porta Mana
## Created: 2022-09-08T17:03:24+0200
## Last-Updated: 2022-09-10T19:49:15+0200
################
## Test script for VB's data analysis
################

## load customized plot functions
source('~/work/pglpm_plotfunctions.R')
##
## Read MCMC seed from command line
mcmcseed = as.integer(commandArgs(trailingOnly=TRUE))[1]
if(is.na(mcmcseed) | (!is.na(mcmcseed) & mcmcseed <=0)){mcmcseed <- 149}
print(paste0('MCMC seed = ',mcmcseed))
##
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
ncores <- 4#availableCores()-1
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
library('nimble')
## NB: also requires library('LaplacesDemon')
#### End custom setup ####

set.seed(707)
## Base name of directory where to save data and plots
baseversion <- '_mcmctest'
nclusters <- 64L
nsamples <- 1024L * 1L # 2L # number of samples AFTER thinning
niter0 <- 1024L * 2L # 3L # iterations burn-in
thin <- 4L #
nstages <- 1L # number of sampling stages beyond burn-in
maincov <- 'Group'
family <- 'Palatino'
#ndata <- 200 # ***set this if you want to use fewer data
posterior <- TRUE # if set to FALSE it samples and plots prior samples
##
## stagestart <- 0L # set this if continuing existing MC = last saved + 1
##

#### INFORMATION ABOUT THE VARIATES AND THEIR PRIOR PARAMETERS
## Pericalcarine r + h, postcentral r, cuneus r,
variateinfo <- data.table(
    variate=c('Site', 'SurfaceHoles', 'Age', 'Group', 'Sex',
              'lh_WM_pericalcarine',
              'lh_GM_pericalcarine',
              'rh_WM_pericalcarine',
              'rh_GM_pericalcarine',
              'rh_WM_postcentral',
              'rh_GM_postcentral',
              'rh_WM_cuneus',
              'rh_GM_cuneus'
              ),
    type=c('integer', 'real', 'real', 'binary', 'binary',
           'real', 'real', 'real', 'real', 'real', 'real', 'real', 'real'), # 'binary' or 'integer' or 'real'
    min=c(1, 0, 12, 0, 0,
          50, 50, 50, 50, 50, 50, 50, 50 ), # 'binary' should have 0
    max=c(8, 130, 100, 1, 0,
          100, 100, 100, 100, 100, 100, 100, 100 ), # 'binary' should have 1
    mean_mean=c(NA, 40, 30, NA, NA,
                75, 75, 75, 75, 75, 75, 75, 75), # only for 'real' variates, NA for others
    mean_sigma=c(NA, 60, 30, NA, NA,
                 20, 20, 20, 20, 20, 20, 20, 20 ),
    sigma_sqrtrate=c(NA, 1, 1, NA, NA,
                     0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001 )/0.3,
    sigma_shape=c(NA, 1, 1, NA, NA,
                  1, 1, 1, 1, 1, 1, 1, 1 )/8 # with 1/4, then min SD value is ~ 0.3 times the one in sigma_sqrtrate
)
## Effects of shape parameter:
## 1/8 (broader):
## > testdata <- -log10(rgamma(n=10^6, shape=1/8, rate=1^2))/2 ; 10^sort(c(quantile(testdata, c(1,7)/8), summary(testdata)))
##         Min.        12.5%      1st Qu.       Median         Mean      3rd Qu.        87.5%         Max. 
## 2.953894e-01 1.942542e+00 3.907964e+00 2.031920e+01 6.598714e+01 3.238969e+02 5.147364e+03 1.529446e+23 
##
## 1/4:
## > testdata <- -log10(rgamma(n=10^6, shape=1/4, rate=1^2))/2 ; 10^sort(c(quantile(testdata, c(1,7)/8), summary(testdata)))
##         Min.        12.5%      1st Qu.       Median         Mean      3rd Qu.        87.5%         Max. 
## 3.206804e-01 1.273500e+00 1.958767e+00 4.786957e+00 8.270111e+00 1.946847e+01 7.731127e+01 2.269262e+12 
##
## 1/2 (narrower):
## > testdata <- -log10(rgamma(n=10^6, shape=1/2, rate=1^2))/2 ; 10^sort(c(quantile(testdata, c(1,7)/8), summary(testdata)))
##         Min.        12.5%      1st Qu.       Median         Mean      3rd Qu.        87.5%         Max. 
## 2.822806e-01 9.213073e-01 1.228308e+00 2.097117e+00 2.669192e+00 4.443665e+00 8.989495e+00 8.112010e+05 



#### FILE WITH DATA
datafile <- 'Cortical_myelination_faux.csv'
##
covNames <- variateinfo$variate
covTypes <- variateinfo$type
covMins <- variateinfo$min
covMaxs <- variateinfo$max
names(covTypes) <- names(covMins) <- names(covMaxs) <- covNames
odata <- fread(datafile, sep=',')
if(!exists('stagestart')){stagestart <- 0L}

## Shuffle
## odata <- odata[sample(1:nrow(odata))]


#################################
## Setup for Monte Carlo sampling
#################################

realCovs <- covNames[covTypes=='real']
integerCovs <- covNames[covTypes=='integer']
binaryCovs <- covNames[covTypes=='binary']
covNames <- c(realCovs, integerCovs, binaryCovs)
nrcovs <- length(realCovs)
nicovs <- length(integerCovs)
nbcovs <- length(binaryCovs)
ncovs <- length(covNames)
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(odata)}
alldata <- odata[1:ndata, ..covNames]
if(stagestart>0){
    continue <- paste0('_finalstate-R',baseversion,'_',mcmcseed,'_',stagestart-1,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nsamples,'.rds')
}
##

##
for(obj in c('constants', 'dat', 'inits', 'bayesnet', 'model', 'Cmodel', 'confmodel', 'mcmcsampler', 'Cmcmcsampler')){if(exists(obj)){do.call(rm,list(obj))}}
gc()

dat <- list()
if(nrcovs>0){ dat$Real=data.matrix(alldata[, ..realCovs])}
if(nicovs>0){ dat$Integer=data.matrix(alldata[, ..integerCovs])}
if(nbcovs>0){ dat$Binary=data.matrix(alldata[, ..binaryCovs])}
##
dirname <- paste0(baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nsamples)
dir.create(dirname)
if(file.exists("/cluster/home/pglpm/R")){
initial.options <- commandArgs(trailingOnly = FALSE)
thisscriptname <- sub('--file=', "", initial.options[grep('--file=', initial.options)])
file.copy(from=thisscriptname, to=paste0(dirname,'/script-R',baseversion,'_',mcmcseed,'_','-V',length(covNames),'-D',ndata,'-K',nclusters,'.Rscript'))
}




####  CONSTANTS, PRIOR PARAMETERS, INITIAL VALUES
source('functions_mcmc.R') # load functions for post-MCMC calculations
##
## In previous versions some statistics of the data were computed
## to decide on the hyperparameters.
## Now this is not done, because wrong in principle
## and because it can lead to silly hyperparameters
##
## Find max integer value in data
if(nicovs > 0){
    ## maximum in data (for inital values)
    maxicovs <- apply(alldata[1:ndata,..integerCovs],2,function(x)max(x, na.rm=T))
    thmaxicovs <- covMaxs[integerCovs] # theoretical maximum
    matrixprobicovs <- matrix(0, nrow=nicovs, ncol=max(thmaxicovs), dimnames=list(integerCovs))
    for(avar in integerCovs){
        matrixprobicovs[avar,1:thmaxicovs[avar]] <- (1:thmaxicovs[avar])/sum(1:thmaxicovs[avar])
    }
}
## constants
constants <- list(nClusters=nclusters)
if(nrcovs>0){constants$nRcovs <- nrcovs}
if(nicovs>0){constants$nIcovs <- nicovs
    constants$maxIcovs <- ncol(matrixprobicovs)}
if(nbcovs>0){constants$nBcovs <- nbcovs}
if(posterior){constants$nData <- ndata}
##
initsFunction <- function(){
    c(list(
        qalpha0=rep(1/nclusters, nclusters) # cluster probabilities
    ),
    if(nrcovs>0){# real variates
        list(
            meanRmean0=variateinfo[variate %in% realCovs, mean_mean],
            meanRtau0=1/variateinfo[variate %in% realCovs, mean_sigma]^2, # dims = inv. variance
            tauRrate0=variateinfo[variate %in% realCovs, sigma_sqrtrate]^2, # dims = variance
            tauRshape0=variateinfo[variate %in% realCovs, sigma_shape]
        )},
    if(nicovs>0){# integer variates
        list(
            probIa0=rep(2, nicovs),
            probIb0=rep(1, nicovs),
            sizeIprob0=matrixprobicovs*matrixprobicovs,
            sizeI=matrix(maxicovs, nrow=nicovs, ncol=nclusters)
        )},
    if(nbcovs>0){# binary variates
        list(
            probBa0=rep(1,nbcovs),
            probBb0=rep(1,nbcovs)
        )},
    if(posterior){list(C=rep(1, ndata))} # cluster occupations: all in one cluster at first
)}


##
#### Mathematical form of the long-run frequency distribution
bayesnet <- nimbleCode({
    q[1:nClusters] ~ ddirch(alpha=qalpha0[1:nClusters])
    ##
    for(acluster in 1:nClusters){
        if(nrcovs>0){# real variates
            for(avar in 1:nRcovs){
                meanR[avar,acluster] ~ dnorm(mean=meanRmean0[avar], tau=meanRtau0[avar])
                tauR[avar,acluster] ~ dgamma(shape=tauRshape0[avar], rate=tauRrate0[avar])
            }
        }
        if(nicovs>0){# integer variates
            for(avar in 1:nIcovs){
                probI[avar,acluster] ~ dbeta(shape1=probIa0[avar], shape2=probIb0[avar])
                sizeI[avar,acluster] ~ dcat(prob=sizeIprob0[avar,1:maxIcovs])
            }
        }
        if(nbcovs>0){# binary variates
            for(avar in 1:nBcovs){
                probB[avar,acluster] ~ dbeta(shape1=probBa0[avar], shape2=probBb0[avar])
            }
        }
    }
    ##
    if(posterior){# cluster occupations
        for(adatum in 1:nData){
            C[adatum] ~ dcat(prob=q[1:nClusters])
        }            ##
        for(adatum in 1:nData){
            if(nrcovs>0){
                for(avar in 1:nRcovs){
                    Real[adatum,avar] ~ dnorm(mean=meanR[avar,C[adatum]], tau=tauR[avar,C[adatum]])
                }
            }
            if(nicovs>0){
                for(avar in 1:nIcovs){
                    Integer[adatum,avar] ~ dbinom(prob=probI[avar,C[adatum]], size=sizeI[avar,C[adatum]])
                }
            }
            if(nbcovs>0){
                for(avar in 1:nBcovs){
                    Binary[adatum,avar] ~ dbern(prob=probB[avar,C[adatum]])
                }
            }
        }
    }
})

##
##
timecount <- Sys.time()

if(posterior){
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=initsFunction(), data=dat, dimensions=list(q=nclusters, meanR=c(nrcovs,nclusters), tauR=c(nrcovs,nclusters), probI=c(nicovs,nclusters), probB=c(nbcovs,nclusters), C=ndata) )
}else{
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=initsFunction(), data=list(), dimensions=list(q=nclusters, meanR=c(nrcovs,nclusters), tauR=c(nrcovs,nclusters), probI=c(nicovs,nclusters), probB=c(nbcovs,nclusters)))
}
Cmodel <- compileNimble(model, showCompilerOutput=FALSE)
gc()


##
if(posterior){# Samplers for posterior sampling
    confmodel <- configureMCMC(Cmodel, nodes=NULL,
                               monitors=c('q',
                                          if(nrcovs > 0){c('meanR', 'tauR')},
                                          if(nicovs > 0){c('probI', 'sizeI')},
                                          if(nbcovs > 0){c('probB')}
                                          ),
                               monitors2=c('C')
                                           )
    ##
    for(adatum in 1:ndata){
        confmodel$addSampler(target=paste0('C[', adatum, ']'), type='categorical')
    }
    for(acluster in 1:nclusters){
        if(nrcovs>0){
            for(avar in 1:nrcovs){
                confmodel$addSampler(target=paste0('meanR[', avar, ', ', acluster, ']'), type='conjugate')
                confmodel$addSampler(target=paste0('tauR[', avar, ', ', acluster, ']'), type='conjugate')
            }
        }
        if(nicovs>0){
            for(avar in 1:nicovs){
                confmodel$addSampler(target=paste0('probI[', avar, ', ', acluster, ']'), type='conjugate')
                confmodel$addSampler(target=paste0('sizeI[', avar, ', ', acluster, ']'), type='categorical')
            }
        }
        if(nbcovs>0){
            for(avar in 1:nbcovs){
                confmodel$addSampler(target=paste0('probB[', avar, ', ', acluster, ']'), type='conjugate')
            }
        }
    }
    confmodel$addSampler(target=paste0('q[1:', nclusters, ']'), type='conjugate')
##
}else{# sampler for prior sampling
    confmodel <- configureMCMC(Cmodel, 
                               monitors=c('q',
                                          if(nrcovs>0){c('meanR', 'tauR')},
                                          if(nicovs>0){c('probI', 'sizeI')},
                                          if(nbcovs>0){c('probB')}
                                          ))
}
##
print(confmodel)

mcmcsampler <- buildMCMC(confmodel)
Cmcmcsampler <- compileNimble(mcmcsampler, resetFunctions = TRUE)
gc()

print('Setup time:')
print(Sys.time() - timecount)

##################################################
## Monte Carlo sampler and plots of MC diagnostics
##################################################
for(stage in stagestart+(0:nstages)){
    totalruntime <- Sys.time()

    print(paste0('==== STAGE ', stage, ' ===='))
    version <- paste0(baseversion,'_',mcmcseed,'_', stage)
    gc()
    if(stage==stagestart){# first sampling stage
        if(exists('continue') && is.character(continue)){# continuing previous
            initsc <- readRDS(paste0(dirname,'/',continue))
            inits0 <- initsFunction()
            for(aname in names(initsc)){inits0[[aname]] <- initsc[[aname]]}
            mcsamples <- runMCMC(Cmcmcsampler, nburnin=0, niter=nsamples*thin, thin=thin, thin2=nsamples*thin, inits=inits0, setSeed=mcmcseed)
        }else{# no previous script runs
            inits0 <- initsFunction
            mcsamples <- runMCMC(Cmcmcsampler, nburnin=1, niter=niter0+1, thin=1, thin2=niter0, inits=inits0, setSeed=mcmcseed)
        }
    }else{# subsequent sampling stages
        Cmcmcsampler$run(niter=nsamples*thin, thin=thin, thin2=nsamples*thin, reset=FALSE, resetMV=TRUE)
    }
    ##
    mcsamples <- as.matrix(Cmcmcsampler$mvSamples)
    print('MCMC time:')
    print(Sys.time() - totalruntime)
    ##
    if(any(is.na(mcsamples))){print('WARNING: SOME NA OUTPUTS')}
    if(any(!is.finite(mcsamples))){print('WARNING: SOME INFINITE OUTPUTS')}
    saveRDS(mcsamples,file=paste0(dirname,'/_mcsamples-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    ##
    ## save final state of MCMC chain
    finalstate <- as.matrix(Cmcmcsampler$mvSamples2)
    finalstate <- c(mcsamples[nrow(mcsamples),], finalstate[nrow(finalstate),])
    ##
    ## Check how many "clusters" were occupied. Warns if too many
    occupations <- finalstate[grepl('^C\\[', names(finalstate))]
    usedclusters <- length(unique(occupations))
    if(usedclusters > nclusters-5){print('WARNING: TOO MANY CLUSTERS OCCUPIED')}
    print(paste0('OCCUPIED CLUSTERS: ', usedclusters, ' OF ', nclusters))
    saveRDS(finalstate2list(finalstate),file=paste0(dirname,'/_finalstate-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    ##
    ## SAVE THE PARAMETERS FOR THE TRANSDUCER
    parmList <- mcsamples2parmlist(mcsamples)
    saveRDS(parmList,file=paste0(dirname,'/_frequencies-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(parmList$q),'.rds'))
    ##
    ## Diagnostics
    ## Log-likelihood
    ll <- llSamples(dat, parmList)
    flagll <- FALSE
    if(!posterior && !any(is.finite(ll))){
        flagll <- TRUE
        ll <- rep(0, length(ll))}
    condprobsd <- logsumsamplesF(Y=data.matrix(na.omit(alldata))[, maincov, drop=F],
                                 X=data.matrix(na.omit(alldata))[, setdiff(covNames, maincov), drop=F],
                                 parmList=parmList, inorder=T)
    condprobsi <- logsumsamplesF(Y=data.matrix(na.omit(alldata))[, setdiff(covNames, maincov), drop=F],
                                 X=data.matrix(na.omit(alldata))[, maincov, drop=F],               
                                 parmList=parmList, inorder=T)
    ##
    traces <- cbind(loglikelihood=ll, 'mean of direct logprobabilities'=condprobsd, 'mean of inverse logprobabilities'=condprobsi)*10/log(10)/ndata #medians, iqrs, Q1s, Q3s,
    badcols <- foreach(i=1:ncol(traces), .combine=c)%do%{if(all(is.na(traces[,i]))){i}else{NULL}}
    if(!is.null(badcols)){traces <- traces[,-badcols]}
    saveRDS(traces,file=paste0(dirname,'/_probtraces-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(parmList$q),'.rds'))
    ##
    if(nrow(traces)>=1000){
        funMCSE <- function(x){LaplacesDemon::MCSE(x, method='batch.means')$se}
    }else{
        funMCSE <- function(x){LaplacesDemon::MCSE(x)}
    }
    diagnESS <- LaplacesDemon::ESS(traces * (abs(traces) < Inf))
    diagnIAT <- apply(traces, 2, function(x){LaplacesDemon::IAT(x[is.finite(x)])})
    diagnBMK <- LaplacesDemon::BMK.Diagnostic(traces, batches=2)[,1]
    diagnMCSE <- 100*apply(traces, 2, function(x){funMCSE(x)/sd(x)})
    diagnStat <- apply(traces, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
    diagnBurn <- apply(traces, 2, function(x){LaplacesDemon::burnin(matrix(x[1:(10*trunc(length(x)/10))], ncol=1))})
    ##
    tracegroups <- list(loglikelihood=1,
                        'main given rest'=2,
                        'rest given main'=3
                        )
    grouplegends <- foreach(agroup=1:length(tracegroups))%do%{
        c( paste0('-- STATS ', names(tracegroups)[agroup], ' --'),
          paste0('min ESS = ', signif(min(diagnESS[tracegroups[[agroup]]]),6)),
          paste0('max IAT = ', signif(max(diagnIAT[tracegroups[[agroup]]]),6)),
          paste0('max BMK = ', signif(max(diagnBMK[tracegroups[[agroup]]]),6)),
          paste0('max MCSE = ', signif(max(diagnMCSE[tracegroups[[agroup]]]),6)),
          paste0('stationary: ', sum(diagnStat[tracegroups[[agroup]]]),'/',length(diagnStat[tracegroups[[agroup]]])),
          paste0('burn: ', signif(max(diagnBurn[tracegroups[[agroup]]]),6))
          )
    }
    colpalette <- c(7,2,1)
    names(colpalette) <- colnames(traces)
    ##
    ## Plot various info and traces
    pdff(paste0(dirname,'/mcsummary2-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(parmList$q)),'a4')
    matplot(1:2, type='l', col='white', main=paste0('Stats stage ',stage), axes=FALSE, ann=FALSE)
    legendpositions <- c('topleft','topright','bottomleft','bottomright')
    for(alegend in 1:length(grouplegends)){
        legend(x=legendpositions[alegend], bty='n', cex=1.5,
               legend=grouplegends[[alegend]] )
    }
    legend(x='center', bty='n', cex=1,
           legend=c(
               paste0('Occupied clusters: ', usedclusters, ' of ', nclusters),
               paste0('LL: ', signif(mean(ll),3), ' +- ', signif(sd(ll),3)),
               'WARNINGS:',
                    if(any(is.na(mcsamples))){'some NA MC outputs'},
                    if(any(!is.finite(mcsamples))){'some infinite MC outputs'},
                    if(usedclusters > nclusters-5){'too many clusters occupied'},
                    if(flagll){'infinite values in likelihood'}
           ))
    ##
    par(mfrow=c(1,1))
    for(avar in covNames){
        datum <- alldata[1:ndata][[avar]]
        if(avar %in% realCovs){
            rg <- c(covMins[avar], covMaxs[avar])
            if(!is.finite(rg[1])){rg[1] <- min(datum, na.rm=T) - IQR(datum, type=8, na.rm=T)}
            if(!is.finite(rg[2])){rg[2] <- max(datum, na.rm=T) + IQR(datum, type=8, na.rm=T)}
            Xgrid <- seq(rg[1], rg[2], length.out=256)
        }else{
            datum <- alldata[1:ndata][[avar]]
            rg <- range(datum, na.rm=T)
            rg <- round(c((covMins[avar]+7*rg[1])/8, (covMaxs[avar]+7*rg[2])/8))
            Xgrid <- rg[1]:rg[2]
        }
        Xgrid <- cbind(Xgrid)
        colnames(Xgrid) <- avar
        plotsamples <- samplesF(Y=Xgrid, parmList=parmList, nfsamples=min(64,nrow(mcsamples)), inorder=FALSE)
        ymax <- quant(apply(plotsamples,2,function(x){quant(x,99/100)}),99/100, na.rm=T)
        tplot(x=Xgrid, y=plotsamples, type='l', col=paste0(palette()[7], '44'), lty=1, lwd=2, xlab=avar, ylab='probability density', ylim=c(0, ymax), family=family)
    }
    ##
    par(mfrow=c(1,1))
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
    dev.off()
    ##
    print('Total runtime:')
    print(Sys.time() - totalruntime)
    ##
}

############################################################
## End MCMC
############################################################
plan(sequential)

