## Author: PGL  Porta Mana
## Created: 2022-09-08T17:03:24+0200
## Last-Updated: 2022-10-21T23:53:50+0200
################
## Exchangeable-probability calculation (non-parametric density regression)
################

#### USER INPUTS AND CHOICES ####
baseversion <- '_testESS' # *** ## Base name of output directory
datafile <- 'data_ep.csv' #***
mainvar <- 'group' # ***
variateinfofile <- 'metadata_noint33_noSW.csv' #***
nsamples <- 1024L * 1L # 2L # number of samples AFTER thinning
thin <- 1L #
nstages <- 35L # number of sampling stages beyond burn-in
niter0 <- 1024L * 1L # 3L # iterations burn-in
## ndata <- 5 # set this if you want to use fewer data
## shuffledata <- TRUE # useful if subsetting data
posterior <- TRUE # if set to FALSE it samples and plots prior samples
showdata <- TRUE # 'histogram' 'scatter' FALSE TRUE
plotdistributions <- TRUE # plot sampled distributions
plotmeans <- FALSE # plot frequency averages
##
nclusters <- 64L
alpha <- 1
compoundgamma <- TRUE # use beta-prime distribution for variance instead of gamma
compoundgammapars <- c(1,1)/2
categoryprior <- 1 # choices: 'Haldane' (1/n) or a number
casualinitvalues <- FALSE
## stagestart <- 3L # set this if continuing existing MC = last saved + 1
family <- 'Palatino'
####


#### Packages and setup ####
## load customized plot functions
if(!exists('tplot')){source('~/work/pglpm_plotfunctions.R')}
##
## Read MCMC seed from command line
mcmcseed = as.integer(commandArgs(trailingOnly=TRUE))[1]
if(is.na(mcmcseed) | (!is.na(mcmcseed) & mcmcseed <=0)){mcmcseed <- 1}
print(paste0('MCMC seed = ',mcmcseed))
##
set.seed(701+mcmcseed)
##
## Packages
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
    ncores <- 1}else{
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
library('nimble')
## NB: also requires libraries 'LaplacesDemon' and 'extraDistr'


#### EXTRACT INFORMATION ABOUT THE VARIATES AND THEIR PRIOR PARAMETERS

if(Sys.info()['nodename']=='luca-HP-Z2-G9'){origdir <- '../'}else{origdir <- ''}
variateinfo <- fread(paste0(origdir,variateinfofile))
## variateinfo <- do.call(rbind, list(
##     ##	 'variate',			'type',		'min',	'max',	'precision')
##     data.table('REANAME',		'real',		-100,	100,	NA)
## ,   data.table('BINNAME',		'binary',	0,	1,	NA)
## ,   data.table('INTNAME',		'integer',	0,	9,	NA)
## ,   data.table('CATNAME',		'category',	1,	10,	NA)
## ))
## colnames(variateinfo) <- c('variate', 'type', 'min', 'max', 'precision')

## Effects of shape parameter:
## 1/8 (broader):
## > testdata <- log10(rinvgamma(n=10^7, shape=1/8, scale=1^2))/2 ; 10^sort(c(quantile(testdata, c(1,7)/8), summary(testdata)))
##         Min.        12.5%      1st Qu.       Median         Mean      3rd Qu. 
## 2.878554e-01 1.937961e+00 3.903828e+00 2.034059e+01 6.641224e+01 3.260648e+02 
##        87.5%         Max. 
## 5.220388e+03 5.266833e+27 
##
## 1/4:
## > testdata <- log10(rinvgamma(n=10^7, shape=1/4, scale=1^2))/2 ; 10^sort(c(quantile(testdata, c(1,7)/8), summary(testdata)))
##         Min.        12.5%      1st Qu.       Median         Mean      3rd Qu. 
## 2.881895e-01 1.274590e+00 1.960271e+00 4.784803e+00 8.283664e+00 1.946374e+01 
##        87.5%         Max. 
## 7.795913e+01 4.670417e+14
##
## 1/2 (narrower):
## > testdata <- log10(rinvgamma(n=10^7, shape=1/2, scale=1^2))/2 ; 10^sort(c(quantile(testdata, c(1,7)/8), summary(testdata)))
##         Min.        12.5%      1st Qu.       Median         Mean      3rd Qu. 
## 2.571008e-01 9.218967e-01 1.229980e+00 2.098062e+00 2.670157e+00 4.440326e+00 
##        87.5%         Max. 
## 8.995125e+00 6.364370e+06 

##
varNames <- variateinfo$variate
varTypes <- variateinfo$type
varMins <- variateinfo$min
varMaxs <- variateinfo$max
names(varTypes) <- names(varMins) <- names(varMaxs) <- varNames
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


###########################################
## READ DATA AND SETUP SOME HYPERPARAMETERS
###########################################

alldata <- fread(paste0(origdir,datafile), sep=',')
if(!all(varNames %in% names(alldata))){print('ERROR: variates missing from datafile')}
alldata <- alldata[, ..varNames]
## shuffle data
if(exists('shuffledata') && shuffledata){alldata <- alldata[sample(1:nrow(alldata))]}
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(alldata)}
alldata <- alldata[1:ndata]
##
if(Sys.info()['nodename']=='luca-HP-Z2-G9'){
    dirname <- ''
}else{
    dirname <- paste0(baseversion,'-V',length(varNames),'-D',ndata,'-K',nclusters,'-I',nsamples,'/')
    dir.create(dirname)
}

## Copy this script to output directory
## if(Sys.info()['nodename']=='luca-HP-Z2-G9'){
##     initial.options <- commandArgs(trailingOnly = FALSE)
##     thisscriptname <- sub('--file=', "", initial.options[grep('--file=', initial.options)])
##     if(mcmcseed==1){file.copy(from=thisscriptname, to=paste0(dirname,'/script-R',baseversion,'-V',length(varNames),'-D',ndata,'-K',nclusters,'.Rscript'))
##     }
## }

#### Get missing metadata from data
## integers
if(any(is.na(varMins[integerVars]))){
    print('WARNING: missing min for some integer variates; using min from data')
    nanames <- integerVars[is.na(varMins[integerVars])]
    for(avar in nanames){
        varMins[avar] <- min(alldata[[avar]], na.rm=T)
    }
}
if(any(is.na(varMaxs[integerVars]))){
    print('WARNING: missing max for some integer variates; using max from data')
    nanames <- integerVars[is.na(varMaxs[integerVars])]
    for(avar in nanames){
        varMaxs[avar] <- max(alldata[[avar]], na.rm=T)
    }
}
## categories
if(any(is.na(varMins[categoryVars]))){
    print('WARNING: missing min for some category variates; using min from data')
    nanames <- categoryVars[is.na(varMins[categoryVars])]
    for(avar in nanames){
        varMins[avar] <- min(alldata[[avar]], na.rm=T)
    }
}
if(any(is.na(varMaxs[categoryVars]))){
    print('WARNING: missing max for some category variates; using max from data')
    nanames <- categoryVars[is.na(varMaxs[categoryVars])]
    for(avar in nanames){
        varMaxs[avar] <- max(alldata[[avar]], na.rm=T)
    }
}



## Normalization and standardization of real variates and calculation of hyperparams
sd2iqr <- 0.5/qnorm(0.75)
##
if(mcmcseed==1){ pdff(paste0(dirname,'densities_variances'),'a4') }
variateparameters <- NULL
for(avar in varNames){
    pars <- c(NA,NA)
    ##
    datum <- alldata[[avar]]
    datamin <- min(datum, na.rm=TRUE)
    datamax <- max(datum, na.rm=TRUE)
    datamedian <- quantile(datum, 0.5, type=8, na.rm=TRUE)
    dataQ1 <- quantile(datum, 0.25, type=8, na.rm=TRUE)
    dataQ2 <- quantile(datum, 0.75, type=8, na.rm=TRUE)
    ##
    if(avar %in% realVars){
        alocation <- median(datum)
        ascale <- IQR(datum) * sd2iqr
        if(ascale==0){ascale <- 1L}
        ## shift and rescale data
        datum <- (datum-alocation)/ascale
        rg <- diff(range(datum, na.rm=T))
        if(is.na(variateinfo[variate==avar, precision])){
            dmin <- min(diff(sort(unique(datum))))
        }else{
            dmin <- variateinfo[variate==avar, precision]/ascale
        }
        ##
        ## ## > ss <- 2; ss2 <- 1/2 ; test <- log10(rinvgamma(1e6, shape=ss, rate=rinvgamma(1e6, shape=ss2, scale=1)))/2; htest <- thist(test); testg <- seq(min(htest$breaks),max(htest$breaks),length.out=256) ; testo <- extraDistr::dbetapr(exp(log(10)*(-2*testg)),shape1=ss,shape2=ss2)*2*log(10)*exp(log(10)*(-2*testg)) ; tplot(list(htest$mids,testg),list(htest$density,testo))
        ## qts <- c(2^-14, 1-2^-14)
        ## fn <- function(p, target){
        ##     sum((log(extraDistr::qbetapr(qts,shape1=p[1],shape2=p[2]))/2 -log(target))^2)
        ##     }
        ## resu <- optim(par=c(2,1/2), fn=fn, target=c(dmin/2,rg*3))
        ## pars <- signif(resu$par, 3)
        ## vals <- sqrt(extraDistr::qbetapr(qts, shape1=pars[1], shape2=pars[2]))
        ## if(abs(vals[1] - dmin/2)/(vals[1] + dmin/2)*200 > 5 |
        ##    abs(vals[2] - rg*3)/(vals[2] + rg*3)*200 > 5
        ##    ){print(paste0('WARNING ', avar, ': bad parameters'))}
        ##
        ## Parameters and plot
        sgrid <- seq(log10(dmin/2), log10(rg*3), length.out=256)
        vv <- exp(log(10)*2*sgrid)
        ##
        if(compoundgamma){
            pars <- compoundgammapars
            ##pars <- c(1,1)
            test <- thist(log10(rinvgamma(1e6, shape=pars[2], rate=rinvgamma(1e6, shape=pars[1], scale=1)))/2)
            vdistr <- extraDistr::dbetapr(x=vv, shape1=pars[1], shape2=pars[2])*vv*2*log(10)
        }else{
            pars <- c(signif(2^4.5,3), 1/2)
            test <- thist(log10(rinvgamma(1e6, shape=pars[2], rate=pars[1]))/2)
            vdistr <- dinvgamma(x=vv, rate=pars[1], shape=pars[2])*vv*2*log(10)
        }
        ##
if(mcmcseed==1){
    tplot(x=list(test$mids,sgrid), y=list(test$density,vdistr),
          xlab=expression(lg~sigma), ylim=c(0,NA), ylab=NA,
          main=paste0(avar, ': shape1/rate = ', pars[1],', shape2/shape = ',pars[2]),
          xlim=log10(c(dmin/2/10, rg*3*10)))
    abline(v=log10(c(dmin/2, dmin, IQR(datum)*sd2iqr, rg, rg*3)), col=yellow, lwd=3)
}
        ##
        rm('test','sgrid','vdistr')
    }else if((avar %in% integerVars) | (avar %in% binaryVars)){
        alocation <- varMins[avar] # integer and binary start from 0
        ascale <- 1L
        dmin <- 1L
        ## shift and rescale data
        datum <- (datum-alocation)
    }else if(avar %in% categoryVars){
        alocation <- varMins[avar] - 1L # category start from 1
        ascale <- 1L
        dmin <- 1L
        ## shift and rescale data
        datum <- (datum-alocation)
    }
    ##
    variateparameters <- rbind(variateparameters,
                               cbind( location=alocation, scale=ascale,
                                 precision=dmin,
                                 min=min(datum, na.rm=T), max=max(datum, na.rm=T),
                                 shape1rate=pars[1], shape2shape=pars[2],
                                 datamin=datamin, datamax=datamax,
                                 datamedian=datamedian,
                                 dataQ1=dataQ1, dataQ2=dataQ2,
                                 thmin=varMins[avar], thmax=varMaxs[avar]
                                 ))
}
rownames(variateparameters) <- varNames
##
if(mcmcseed==1){
    dev.off()
    fwrite(cbind(data.table(variate=varNames),variateparameters), file=paste0(dirname,'variateparameters.csv'))
}
##
locvarMins <- varMins-variateparameters[,'location']
locvarMaxs <- varMaxs-variateparameters[,'location']


#################################
## Setup for Monte Carlo sampling
#################################

if(!exists('stagestart')){stagestart <- 0L}
if(stagestart>0){
    continue <- paste0('_finalstate-R',baseversion,'-V',length(varNames),'-D',ndata,'-K',nclusters,'-I',nsamples,'--',stagestart-1,'-',mcmcseed,'.rds')
}else{
    continue <- FALSE
}
##


for(obj in c('constants', 'dat', 'inits', 'bayesnet', 'model', 'Cmodel', 'confmodel', 'mcmcsampler', 'Cmcmcsampler')){if(exists(obj)){do.call(rm,list(obj))}}
gc()


## Data (standardized for real variates)
dat <- list()
if(nrvars>0){ dat$Real=t((t(data.matrix(alldata[, ..realVars])) - variateparameters[realVars,'location'])/variateparameters[realVars,'scale'])}
## if(nivars>0){ dat$Integer=data.matrix(alldata[, ..integerVars])}
if(nivars>0){ dat$Integer=t((t(data.matrix(alldata[, ..integerVars])) - variateparameters[integerVars,'location']))}
## if(ncvars>0){ dat$Category=data.matrix(alldata[, ..categoryVars])}
if(ncvars>0){ dat$Category=t((t(data.matrix(alldata[, ..categoryVars])) - variateparameters[categoryVars,'location']))}
## if(nbvars>0){ dat$Binary=data.matrix(alldata[, ..binaryVars])}
if(nbvars>0){ dat$Binary=t((t(data.matrix(alldata[, ..binaryVars])) - variateparameters[binaryVars,'location']))}



####  CONSTANTS, PRIOR PARAMETERS, INITIAL VALUES
source(paste0(origdir,'functions_mcmc.R')) # load functions for post-MCMC calculations
##
## In previous versions some statistics of the data were computed
## to decide on the hyperparameters.
## Now this is not done, because wrong in principle
## and because it can lead to silly hyperparameters
##
## Find max integer value in data
if(nivars > 0){
    ## maximum in data (for initial values)
    maxivars <- variateparameters[integerVars, 'max']
    thmaxivars <- locvarMaxs[integerVars] # theoretical maximum
    matrixprobivars <- matrix(0, nrow=nivars, ncol=max(thmaxivars), dimnames=list(integerVars))
    for(avar in integerVars){
        matrixprobivars[avar,1:thmaxivars[avar]] <- (1:thmaxivars[avar])/sum(1:thmaxivars[avar])
    }
}
##
## Find max number of categories in data
if(ncvars > 0){
    ncategories <- max(locvarMaxs[categoryVars]) # largest number of categories
    calphapad <- array(2^(-40), dim=c(ncvars, ncategories), dimnames=list(categoryVars,NULL))
    for(avar in categoryVars){
        if(categoryprior=='Haldane'){
            calphapad[avar,1:locvarMaxs[avar]] <- 1/locvarMaxs[avar]
        }else{
            calphapad[avar,1:locvarMaxs[avar]] <- categoryprior
        }
    }
}
## constants
constants <- list(nClusters=nclusters)
if(nrvars>0){constants$nRvars <- nrvars}
if(nivars>0){constants$nIvars <- nivars
    constants$maxIvars <- ncol(matrixprobivars)}
if(ncvars>0){constants$nCvars <- ncvars
    constants$nCategories <- ncategories}
if(nbvars>0){constants$nBvars <- nbvars}
if(posterior){constants$nData <- ndata}
##
initsFunction <- function(){
    c(list(
        qalpha0 = rep(alpha/nclusters, nclusters) # cluster probabilities
    ),
    if(nrvars > 0){# real variates
        list(
            meanRmean0 = variateparameters[realVars,'location']*0,
            meanRvar0 = (3 * (variateparameters[realVars,'scale']*0+1))^2,
            varRshape2shape = variateparameters[realVars,'shape2shape']
        )},
    if(compoundgamma & nrvars > 0){
        list(varRshape1 = variateparameters[realVars,'shape1rate'])
    }else{
        list(varRrate = variateparameters[realVars,'shape1rate'])
    },
    if(nivars > 0){# integer variates
        list(
            probIa0 = (log(thmaxivars)*15000-23000)^0.1,
            probIb0 = rep(1, nivars),
            sizeIprob0 = matrixprobivars,
            ## sizeI = matrix(c(maxivars,rep(NA,nivars*(nclusters-1))), nrow=nivars, ncol=nclusters)
            sizeI = cbind(maxivars, t(apply(matrixprobivars,1,function(prob){rcat(nclusters-1,prob=prob)})))
        )},
    if(ncvars > 0){# categorical variates
        list(
            calpha0 = calphapad
        )},
    if(nbvars > 0){# binary variates
        list(
            probBa0 = rep(1,nbvars),
            probBb0 = rep(1,nbvars)
        )},
    if((!casualinitvalues) & posterior){list(
                                            q = rep(1/nclusters, nclusters),
                                            C = rep(1, ndata)  # cluster occupations: all in one cluster at first
                                        )},
    if(casualinitvalues & posterior){
        list(q = rdirch(1, alpha=rep(1,nclusters)),
             ## C = rep(1, ndata))
             C = sample(1:nclusters, ndata, replace=TRUE))
        }
)}


##
#### Mathematical form of the long-run frequency distribution
bayesnet <- nimbleCode({
    q[1:nClusters] ~ ddirch(alpha=qalpha0[1:nClusters])
    ##
    if(nrvars>0){# real variates
        if(compoundgamma){
            for(avar in 1:nRvars){
                varRrate[avar] ~ dinvgamma(shape=varRshape1[avar], scale=1)
            }
        }
    }
    ##
    for(acluster in 1:nClusters){
        if(nrvars>0){# real variates
            for(avar in 1:nRvars){
                meanR[avar,acluster] ~ dnorm(mean=meanRmean0[avar], var=meanRvar0[avar])
                varR[avar,acluster] ~ dinvgamma(shape=varRshape2shape[avar], rate=varRrate[avar])
            }
        }
        if(nivars>0){# integer variates
            for(avar in 1:nIvars){
                probI[avar,acluster] ~ dbeta(shape1=probIa0[avar], shape2=probIb0[avar])
                sizeI[avar,acluster] ~ dcat(prob=sizeIprob0[avar,1:maxIvars])
            }
        }
        if(ncvars>0){# category variates
            for(avar in 1:nCvars){
                probC[avar,acluster,1:nCategories] ~ ddirch(alpha=calpha0[avar,1:nCategories])
            }
        }
        if(nbvars>0){# binary variates
            for(avar in 1:nBvars){
                probB[avar,acluster] ~ dbeta(shape1=probBa0[avar], shape2=probBb0[avar])
            }
        }
    }
    ##
    if(posterior){# cluster occupations
        for(adatum in 1:nData){
            C[adatum] ~ dcat(prob=q[1:nClusters])
            ##
            if(nrvars>0){# real variates
                for(avar in 1:nRvars){
                    Real[adatum,avar] ~ dnorm(mean=meanR[avar,C[adatum]], var=varR[avar,C[adatum]])
                }
            }
            if(nivars>0){# integer variates
                for(avar in 1:nIvars){
                    Integer[adatum,avar] ~ dbinom(prob=probI[avar,C[adatum]], size=sizeI[avar,C[adatum]])
                }
            }
            if(ncvars>0){# category variates
                for(avar in 1:nCvars){
                    Category[adatum,avar] ~ dcat(prob=probC[avar,C[adatum],1:nCategories])
                }
            }
            if(nbvars>0){# binary variates
                for(avar in 1:nBvars){
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
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants,
                         inits=initsFunction(), data=dat,
                         dimensions=c(
                             list(q=nclusters),
                             (if(nrvars>0){list(meanR=c(nrvars,nclusters), tauR=c(nrvars,nclusters))}),
                             (if(nivars>0){list(probI=c(nivars,nclusters))}),
                             (if(ncvars>0){list(probC=c(ncvars,nclusters,ncategories))}),
                             (if(nbvars>0){list(probB=c(nbvars,nclusters))}),
                             list(C=ndata),
                             if(compoundgamma){list(varRrate=nrvars)})
                         )
}else{
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants,
                         inits=initsFunction(), data=list(),
                         dimensions=c(
                             list(q=nclusters),
                             (if(nrvars>0){list(meanR=c(nrvars,nclusters), tauR=c(nrvars,nclusters))}),
                             (if(nivars>0){list(probI=c(nivars,nclusters))}),
                             (if(ncvars>0){list(probC=c(ncvars,nclusters,ncategories))}),
                             (if(nbvars>0){list(probB=c(nbvars,nclusters))}),
                             if(compoundgamma){list(varRrate=nrvars)})
                         )
}
Cmodel <- compileNimble(model, showCompilerOutput=FALSE)
gc()


##
if(posterior){# Samplers for posterior sampling
    confmodel <- configureMCMC(Cmodel, nodes=NULL,
                               monitors=c('q',
                                          if(nrvars > 0){c('meanR', 'varR')},
                                          if(nivars > 0){c('probI', 'sizeI')},
                                          if(ncvars > 0){c('probC')},
                                          if(nbvars > 0){c('probB')}
                                          ),
                               monitors2=c('C',
                                           if(compoundgamma & nrvars > 0){c('varRrate')}
                                           )
                               )
    ##
    for(adatum in 1:ndata){
        confmodel$addSampler(target=paste0('C[', adatum, ']'), type='categorical')
    }
    if(compoundgamma & nrvars>0){
        for(avar in 1:nrvars){
            confmodel$addSampler(target=paste0('varRrate[', avar, ']'), type='conjugate')
        }
    }
    for(acluster in 1:nclusters){
        if(nrvars>0){
            for(avar in 1:nrvars){
                confmodel$addSampler(target=paste0('varR[', avar, ', ', acluster, ']'), type='conjugate')
                confmodel$addSampler(target=paste0('meanR[', avar, ', ', acluster, ']'), type='conjugate')
            }
        }
        if(nivars>0){
            for(avar in 1:nivars){
                confmodel$addSampler(target=paste0('probI[', avar, ', ', acluster, ']'), type='conjugate')
                confmodel$addSampler(target=paste0('sizeI[', avar, ', ', acluster, ']'), type='categorical')
            }
        }
        if(ncvars>0){
            for(avar in 1:ncvars){
                confmodel$addSampler(target=paste0('probC[', avar, ', ', acluster, ', 1:', ncategories, ']'), type='conjugate')
            }
        }
        if(nbvars>0){
            for(avar in 1:nbvars){
                confmodel$addSampler(target=paste0('probB[', avar, ', ', acluster, ']'), type='conjugate')
            }
        }
    }
    confmodel$addSampler(target=paste0('q[1:', nclusters, ']'), type='conjugate')
##
}else{# sampler for prior sampling
    confmodel <- configureMCMC(Cmodel, 
                               monitors=c('q',
                                          if(nrvars>0){c('meanR', 'varR')},
                                          if(nivars>0){c('probI', 'sizeI')},
                                          if(ncvars>0){c('probC')},
                                          if(nbvars>0){c('probB')}
                                          ),
                               monitors2=c(if(compoundgamma & nrvars > 0){c('varRrate')}
                                           )
                               )
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
for(stage in is.character(continue):1){

    calctime <- Sys.time()
    print(paste0('==== STAGE ', stage, ' ===='))
    gc()
    if(stage==0){# burn-in stage
        inits0 <- initsFunction
        mcsamples <- runMCMC(Cmcmcsampler, nburnin=1, niter=niter0+1, thin=1, thin2=niter0, inits=inits0, setSeed=mcmcseed+stage+100)
    }else if(is.character(continue)){# continuing previous
        initsc <- readRDS(paste0(dirname,continue))
        inits0 <- initsFunction()
        for(aname in names(inits0)){inits0[[aname]] <- initsc[[aname]]}
        thin <- initsc[['thin']]
        mcsamples <- runMCMC(Cmcmcsampler, nburnin=0, niter=nsamples*thin, thin=thin, thin2=nsamples*thin, inits=inits0, setSeed=mcmcseed+stage+100)
    }else{# subsequent sampling stages
        Cmcmcsampler$run(niter=nsamples*thin, thin=thin, thin2=nsamples*thin, reset=FALSE, resetMV=TRUE)
    }
    ##
    mcsamples <- as.matrix(Cmcmcsampler$mvSamples)
    print('Time MCMC:')
    print(Sys.time() - calctime)
    ##
    if(any(is.na(mcsamples))){print('WARNING: SOME NA OUTPUTS')}
    if(any(!is.finite(mcsamples))){print('WARNING: SOME INFINITE OUTPUTS')}
    saveRDS(mcsamples,file=paste0(dirname,'_mcsamples-R',baseversion,'-V',length(varNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'--',stage,'-',mcmcseed,'.rds'))
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
    saveRDS(finalstate2list(finalstate, realVars=realVars, integerVars=integerVars, categoryVars=categoryVars, binaryVars=binaryVars, compoundgamma=compoundgamma), file=paste0(dirname,'_finalstate-R',baseversion,'-V',length(varNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'--',stage,'-',mcmcseed,'.rds'))
    ##
    ## SAVE THE PARAMETERS
    parmList <- mcsamples2parmlist(mcsamples, realVars, integerVars, categoryVars, binaryVars)
    saveRDS(parmList,file=paste0(dirname,'_frequencies-R',baseversion,'-V',length(varNames),'-D',ndata,'-K',nclusters,'-I',nrow(parmList$q),'--',stage,'-',mcmcseed,'.rds'))
    ##
    ## Diagnostics
    ## Log-likelihood
    if(posterior){
        ll <- llSamples(dat, parmList)
        flagll <- FALSE
        if(!posterior && !any(is.finite(ll))){
            flagll <- TRUE
            ll <- rep(0, length(ll))}
        condprobsd <- logsumsamplesF(Y=do.call(cbind,dat)[, mainvar, drop=F],
                                     X=do.call(cbind,dat)[, setdiff(varNames, mainvar),
                                                          drop=F],
                                     parmList=parmList, inorder=T)
        condprobsi <- logsumsamplesF(Y=do.call(cbind,dat)[, setdiff(varNames, mainvar),
                                                          drop=F],
                                     X=do.call(cbind,dat)[, mainvar, drop=F],
                                     parmList=parmList, inorder=T)
        ##
        traces <- cbind(loglikelihood=ll, 'mean of direct logprobabilities'=condprobsd, 'mean of inverse logprobabilities'=condprobsi)*10/log(10)/ndata #medians, iqrs, Q1s, Q3s
        badcols <- foreach(i=1:ncol(traces), .combine=c)%do%{if(all(is.na(traces[,i]))){i}else{NULL}}
        if(!is.null(badcols)){traces <- traces[,-badcols]}
        saveRDS(traces,file=paste0(dirname,'_probtraces-R',baseversion,'-V',length(varNames),'-D',ndata,'-K',nclusters,'-I',nrow(parmList$q),'--',stage,'-',mcmcseed,'.rds'))
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
        if(stage==0){ thin <- round(diagnIAT[1]) }
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
    }
    ##
    ## Plot various info and traces
    pdff(paste0(dirname,'mcmcsummary-R',baseversion,'-V',length(varNames),'-D',ndata,'-K',nclusters,'-I',nrow(parmList$q),'--',stage,'-',mcmcseed),'a4')
    ## Summary stats
    if(posterior){
        matplot(1:2, type='l', col='white', main=paste0('Stats stage ',stage), axes=FALSE, ann=FALSE)
        legendpositions <- c('topleft','topright','bottomleft','bottomright')
        for(alegend in 1:length(grouplegends)){
            legend(x=legendpositions[alegend], bty='n', cex=1.5,
                   legend=grouplegends[[alegend]] )
        }
        legend(x='center', bty='n', cex=1,
               legend=c(
                   paste0('Occupied clusters: ', usedclusters, ' of ', nclusters),
                   paste0('LL:  ( ', signif(mean(traces[,1]),3), ' +- ', signif(sd(traces[,1]),3),' ) dHart'),
                   'WARNINGS:',
                   if(any(is.na(mcsamples))){'some NA MC outputs'},
                   if(any(!is.finite(mcsamples))){'some infinite MC outputs'},
                   if(usedclusters > nclusters-5){'too many clusters occupied'},
                   if(flagll){'infinite values in likelihood'}
               ))
        ##
        ## Traces of likelihood and cond. probabilities
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
    }
    ##
    ## Samples of marginal frequency distributions
    if((stage > 0 & plotdistributions) |
       (stage==0 & plotburnindistributions)){
        if(plotmeans){nfsamples <- totsamples
        }else if(!posterior){nfsamples <- 256
        }else{nfsamples <- 63}
        subsample <- round(seq(1,nfsamples, length.out=63))
        for(avar in varNames){#print(avar)
            if(avar %in% realVars){
                rg <- signif((variateparameters[avar,c('thmin','thmax')] + 
                              7*variateparameters[avar,c('datamin','datamax')])/8, 2)
                if(!is.finite(rg[1])){rg[1] <- diff(variateparameters[avar,c('scale','datamin')])}
                if(!is.finite(rg[2])){rg[2] <- sum(variateparameters[avar,c('scale','datamax')])}
                Xgrid <- cbind(seq(rg[1], rg[2], length.out=256))
            }else{
                rg <- (variateparameters[avar,c('thmin','thmax')] + 
                       7*variateparameters[avar,c('datamin','datamax')])/8
                rg <- c(floor(rg[1]), ceiling(rg[2]))
                Xgrid <- cbind(rg[1]:rg[2])
            }
            colnames(Xgrid) <- avar
            plotsamples <- samplesF(Y=Xgrid, parmList=parmList, nfsamples=min(nfsamples,nrow(mcsamples)), inorder=FALSE, transform=variateparameters)
            fiven <- variateparameters[avar,c('datamin','dataQ1','datamedian','dataQ2','datamax')]
            ##
            if(posterior){
                par(mfrow=c(1,1))
                ymax <- quant(apply(plotsamples[,subsample],2,function(x){quant(x,99/100)}),99/100, na.rm=T)
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
            }else{
                par(mfrow=c(8,8),mar = c(0,0,0,0))
                tplot(x=list(Xgrid,Xgrid), y=list(rowMeans(plotsamples),rep(0,length(Xgrid))), type='l', col=c(paste0(palette()[3], 'FF'), '#bbbbbb80'), lty=1, lwd=c(2,1), xlab=NA, ylab=NA, ylim=c(0, NA), family=family,
                      xticks=NA, yticks=NA,
                      mar=c(1,1,1,1))
                abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'))
                text(sum(range(Xgrid))/2, par('usr')[4]*0.9, avar)
                ##
                for(aplot in 1:63){
                    tplot(x=list(Xgrid,Xgrid), y=list(plotsamples[,subsample[aplot]], rep(0,length(Xgrid))), type='l', col=c(paste0(palette()[1], 'FF'), '#bbbbbb80'), lty=1, lwd=c(1,1), xlab=NA, ylab=NA, ylim=c(0, NA), family=family,
                          xticks=NA, yticks=NA,
                          mar=c(1,1,1,1))
                    abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'))
                    if(aplot==1){ text(Xgrid[1], par('usr')[4]*0.9, variateinfo[variate==avar,type], pos=4)}
                    if(aplot==2){ text(Xgrid[1], par('usr')[4]*0.9, paste(signif(c(rg,diff(rg)),2),collapse=' -- '), pos=4)}
                }
            }
        }
    }
    dev.off()

    ##
    print('Time MCMC+diagnostics:')
    print(Sys.time() - calctime)
    ##
}

############################################################
## End MCMC
############################################################
plan(sequential)

stop('NONE. End of script')

alldata <- fread('data_ep.csv')
metadata <- fread('metadata.csv')
allnames <- names(alldata)
metadatanew <- metadata
metadatanew$precision <- as.integer(metadatanew$precision)
##
boundarycat <- 10
boundaryrea <- 100
for(i in allnames){
    if(metadata[variate==i,'type']=='integer'){
        print(i)
        if(is.na(metadata[variate==i,'max'])){
            print(paste0('(max)'))
            metadatanew[variate==i,'max'] <- max(alldata[[i]],na.rm=T)
        }
        if(is.na(metadata[variate==i,'min'])){
            print(paste0(': min'))
            metadatanew[variate==i,'min'] <- min(alldata[[i]],na.rm=T)
        }
        rg <- metadatanew[variate==i,'max']-metadatanew[variate==i,'min']+1
        if(rg < boundarycat){
            metadatanew[variate==i,'type'] <- 'categorical'
            print(paste0(': to cat, ',rg))
        }else if(rg>=boundaryrea){
            metadatanew[variate==i,'type'] <- 'real'
            print(paste0(': to real, ',rg))
            metadatanew[variate==i,'precision'] <- 1L
        }else{print(': no changes')}
    }
}

fwrite(x=metadatanew, file=paste0('metadata_noint',boundarycat,'_',boundaryrea,'.csv'))

nosw <- metadatanew$variate[grepl('_SW$',metadatanew$variate)]
metadatanew2 <- metadatanew[!(variate %in% nosw),]

fwrite(x=metadatanew2, file=paste0('metadata_noint',boundarycat,'_',boundaryrea,'_noSW.csv'))


for(i in allnames){if(metadata[variate==i,type]=='integer' & !(min(diff(sort(unique(alldata[[i]]))))==1)){print(c(i,min(diff(sort(unique(alldata[[i]]))))))}}



traces <- traces1
t(sapply(1:17,function(thin){c(
                                 thin,
                                 LaplacesDemon::IAT(traces[1:1024,1]),
                                 LaplacesDemon::IAT(traces[1025:2048,1]),
                                 rev(thisseq <- seq(from=2048+1,by=thin,length.out=2048))[1],
                                 thisiat <- LaplacesDemon::IAT(traces[thisseq,1]),
                                 thisess <- LaplacesDemon::ESS(traces[thisseq,1]),
                                 thisiat <- LaplacesDemon::IAT(traces[thisseq,1]),
                                 thisess <- LaplacesDemon::ESS(traces[thisseq,1]),
                                 floor(thisiat*thin),round(thisiat)*thin,
                                 NA,
                                 floor(length(thisseq)/thisess*thin),round(length(thisseq)/thisess)*thin
                             )}))


traces <- traces2
thin <-  round(1024/LaplacesDemon::ESS(traces[1:1024,1]))
thisseq <- seq(from=1024+1,by=thin,length.out=1024)
c(thin, range(thisseq), max(thisseq)/1024)
LaplacesDemon::IAT(traces[thisseq,1])
LaplacesDemon::ESS(traces[thisseq,1])
tplot(y=traces[thisseq,1])

traces <- traces2
leng <- 2048
thin <-  round(LaplacesDemon::IAT(traces[1:1024,1]))
thisseq <- seq(from=1024+1,by=thin,length.out=leng)
thisseq <- round(max(thisseq)/4)+thisseq-1024
x <- traces[thisseq,1]
c(thin, range(thisseq), max(thisseq)/1024)
LaplacesDemon::BMK.Diagnostic(cbind(x), batches=2)[,1]
LaplacesDemon::BMK.Diagnostic(cbind(x), batches=4)[,1]
LaplacesDemon::IAT(x)
LaplacesDemon::ESS(x)
LaplacesDemon::MCSE(x, method='batch.means')$se*100/sd(x) #6.27
LaplacesDemon::MCSE(x, method='batch.means')$se*100/sd(x) < 6.27
LaplacesDemon::is.stationary(cbind(x))
tplot(y=traces[thisseq,1])

##traces <- traces1
thin <-  round(LaplacesDemon::IAT(traces[1:1024,1]))
t(sapply(1:16,function(i){
    thisseq <- 1024+seq(from=128*(i-1)+1,by=thin,length.out=128)
    x <- traces[thisseq,1]
    c(thin, range(thisseq), max(thisseq)/128,NA,
    LaplacesDemon::BMK.Diagnostic(cbind(x), batches=2)[,1],
    LaplacesDemon::BMK.Diagnostic(cbind(x), batches=4)[,1],
    LaplacesDemon::IAT(x),
    LaplacesDemon::ESS(x),
    LaplacesDemon::MCSE(x, method='batch.means')$se*100/sd(x), #6.27
    LaplacesDemon::MCSE(x, method='batch.means')$se*100/sd(x) < 6.27,
    LaplacesDemon::is.stationary(cbind(x))
)}))

    
}
thisseq <- round(max(thisseq)/2)+thisseq-1024
tplot(y=traces[thisseq,1])





traces <- traces2
t(sapply(1:36,function(i){
    thisseq <- seq(from=1024*(i-1)+1,by=1,length.out=1024)
    x <- traces[thisseq,1]
    c(i,max(thisseq),
      LaplacesDemon::BMK.Diagnostic(cbind(x), batches=2)[,1],
      LaplacesDemon::BMK.Diagnostic(cbind(x), batches=4)[,1],
      LaplacesDemon::IAT(x),
      LaplacesDemon::ESS(x),
      LaplacesDemon::MCSE(x, method='batch.means')$se*100/sd(x),
      LaplacesDemon::is.stationary(cbind(x))
)}))
