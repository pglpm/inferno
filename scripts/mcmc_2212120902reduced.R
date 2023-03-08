## Author: PGL  Porta Mana
## Created: 2022-09-08T17:03:24+0200
## Last-Updated: 2022-12-13T18:56:34+0100
#########################################
## Inference of exchangeable variates (nonparametric density regression)
## using effectively-infinite mixture of product kernels
## Monte Carlo sampling
#########################################

#### USER INPUTS AND CHOICES ####
baseversion <- '_ingrid3' # *** ## Base name of output directory
## datafile <- 'testdataS1.csv'#'ingrid_data_nogds6.csv' #***
datafile <- 'ingrid_data_nogds6.csv' #***
predictorfile <- 'predictors.csv'
predictandfile <- NULL # 'predictors.csv'
varinfofile <- 'varinfo.rds'
requiredESS <- 1024*2/20 # required effective sample size
nsamples <- 8*ceiling((requiredESS*1.5)/8) # number of samples AFTER thinning
ndata <- NULL # set this if you want to use fewer data
shuffledata <- FALSE # useful if subsetting data
posterior <- TRUE # if set to FALSE it samples and plots prior samples
minstepincrease <- 8L
savetempsamples <- FALSE # save temporary MCMC samples
plottempdistributions <- FALSE # plot temporary sampled distributions
showdata <- TRUE # 'histogram' 'scatter' FALSE TRUE
plotmeans <- FALSE # plot frequency averages
totsamples <- 100 # number of samples if plotting frequency averages
##
niter0 <- 1024L * 1L # 3L # iterations burn-in
nclusters <- 64L
alpha0 <- 2^((-3):3)
casualinitvalues <- FALSE
## stagestart <- 3L # set this if continuing existing MC = last saved + 1
family <- 'Palatino'
####


#### Packages and setup ####
## load customized plot functions
if(!exists('tplot')){source('~/work/pglpm_plotfunctions.R')}
##
## Read MCMC seed from command line
mcmcseed <- as.integer(commandArgs(trailingOnly=TRUE))[1]
if(is.na(mcmcseed) || (!is.na(mcmcseed) && mcmcseed <= 0)){mcmcseed <- 1}
cat(paste0('\nMCMC seed = ',mcmcseed,'\n'))
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
cat('\navailableCores: ')
cat(availableCores())
cat('\navailableCores-multicore: ')
cat(availableCores('multicore'))
if(Sys.info()['nodename']=='luca-HP-Z2-G9'){
    ncores <- 1}else{
    ncores <- 6}
cat(paste0('\nusing ',ncores,' cores\n'))
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


###############################################
## READ DATA, META-INFORMATION, HYPERPARAMETERS
###############################################

if(Sys.info()['nodename']=='luca-HP-Z2-G9'){origdir <- '../'}else{origdir <- ''}
source(paste0(origdir,'functionsmcmc_2212120902.R')) # load functions for post-MCMC calculations
## varinfo <- data.matrix(read.csv(paste0(origdir,varinfofile), row.names=1))
varinfo <- readRDS(paste0(origdir,varinfofile))

variate <- lapply(variatetypes, function(x)names(varinfo[['type']])[varinfo[['type']]==x])
len <- lapply(variate,length)
names(variate) <- names(len) <- variatetypes

nalpha <- length(alpha0)

data0 <- fread(paste0(origdir,datafile), sep=',')
if(!all(unlist(variate) %in% colnames(data0))){cat('\nERROR: variates missing from datafile')}
data0 <- data0[, unlist(variate), with=F]
## shuffle data
if(exists('shuffledata') && shuffledata){data0 <- data0[sample(1:nrow(data0))]}
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(data0)}
data0 <- data0[1:ndata]

basename <- paste0(baseversion,'-V',sum(unlist(len)),'-D',ndata,'-K',nclusters,'-I',nsamples)
##
if(Sys.info()['nodename']=='luca-HP-Z2-G9'){
    dirname <- ''
}else{
    dirname <- paste0(basename,'/')
    dir.create(dirname)
}

if(mcmcseed == 1){saveRDS(varinfo,file=paste0(dirname,'_varinfo-R',basename,'.rds'))}

source(paste0(origdir,'functionsmcmc_2212120902.R')) # load functions for post-MCMC ca

if(!is.null(predictorfile)){predictorfile <- paste0(origdir,predictorfile) }
if(!is.null(predictandfile)){predictandfile <- paste0(origdir,predictandfile) }

#################################
## Setup for Monte Carlo sampling
#################################

## if(!exists('stagestart')){stagestart <- 0L}
## if(stagestart>0){
##     resume <- paste0('_finalstate-R',baseversion,'-V',sum(unlist(len)),'-D',ndata,'-K',nclusters,'-I',nsamples,'--',stagestart-1,'-',mcmcseed,'.rds')
## }else{
##     resume <- FALSE
## }
##

for(obj in c('constants', 'datapoints', 'inits', 'finitemix', 'finitemixnimble', 'Cfinitemixnimble', 'confnimble', 'mcsampler', 'Cmcsampler')){if(exists(obj)){do.call(rm,list(obj))}}
gc()

## data
datapoints = c(
    ## real
    if(len$R > 0){list( Rdata = transf(data0[,variate$R,with=F], varinfo) )},
    ## one-censored
    if(len$O > 0){list(
                   Odata = transf(data0[,variate$O,with=F], varinfo, Oout='data'),
                   Oleft = transf(data0[,variate$O,with=F], varinfo, Oout='left'),
                   Oaux = matrix(1L, nrow=ndata, ncol=len$O)
                  )},
    ## Two-bounded
    if(len$D > 0){list(
                   Ddata = transf(data0[,variate$D,with=F], varinfo, Dout='data'),
                   Dleft = transf(data0[,variate$D,with=F], varinfo, Dout='left'),
                   Dright = transf(data0[,variate$D,with=F], varinfo, Dout='right'),
                   Daux = matrix(1L, nrow=ndata, ncol=len$D)
                  )},
    ## integer
    if(len$I > 0){list(
                   Ileft = transf(data0[,variate$I,with=F], varinfo, Iout='left'),
                   Iright = transf(data0[,variate$I,with=F], varinfo, Iout='right'),
                   Iaux = matrix(1L, nrow=ndata, ncol=len$I)
    )},
    ## binary
    if(len$B > 0){list( Bdata = transf(data0[,variate$B,with=F], varinfo) )},
    ## categorical
    if(len$C > 0){list( Cdata = transf(data0[,variate$C,with=F], varinfo) )}
)

##
## constants
constants <- c(
    list(nclusters = nclusters),
    if(nalpha > 0){ list(nalpha = nalpha) },
    if(len$R > 0){ list(Rn = len$R) },
    if(len$O > 0){ list(On = len$O) },
    if(len$D > 0){ list(Dn = len$D) },
    if(len$I > 0){ list(In = len$I) },
    if(len$B > 0){ list(Bn = len$B) },
    if(len$C > 0){ list(Cn = len$C, Cmaxn = Cmaxn) },
    if(posterior){ list(ndata = ndata)}
)

##
## hyperparameters and some initial values
initsFunction <- function(){
    c(
        if(nalpha > 1){list( # distribution over concentration parameter
                           Alphaindex = length(alpha0),
                           probalpha0 = rep(1/nalpha, nalpha),
                           walpha0 = matrix(alpha0/nclusters,
                                            nrow=nalpha, ncol=nclusters)
                       )}else{list(
                                  walpha0 = rep(1/nclusters, nclusters)
                              )},
        if(len$R > 0){list( # real variate
                     Rmean0 = varinfo[[ 'hmean']][variate$R],
                     Rvar0 = varinfo[[ 'hsd']][variate$R]^2,
                     Rshapeout0 = varinfo[[ 'hshapeout']][variate$R],
                     Rshapein0 = varinfo[[ 'hshapein']][variate$R],
                     Rvarscale0 = varinfo[[ 'hvarscale']][variate$R]^2
                 )},
        if(len$O > 0){list( # doubly-bounded variate
                     Odata = transf(data0[,variate$O,with=F], varinfo, Oout='init'),
                     Omean0 = varinfo[[ 'hmean']][variate$O],
                     Ovar0 = varinfo[[ 'hsd']][variate$O]^2,
                     Oshapeout0 = varinfo[[ 'hshapeout']][variate$O],
                     Oshapein0 = varinfo[[ 'hshapein']][variate$O],
                     Ovarscale0 = varinfo[[ 'hvarscale']][variate$O]^2
                 )},
        if(len$D > 0){list( # doubly-bounded variate
                     Ddata = transf(data0[,variate$D,with=F], varinfo, Dout='init'),
                     Dmean0 = varinfo[[ 'hmean']][variate$D],
                     Dvar0 = varinfo[[ 'hsd']][variate$D]^2,
                     Dshapeout0 = varinfo[[ 'hshapeout']][variate$D],
                     Dshapein0 = varinfo[[ 'hshapein']][variate$D],
                     Dvarscale0 = varinfo[[ 'hvarscale']][variate$D]^2
                 )},
        if(len$I > 0){list( # integer ordinal variate
                     Icont = transf(data0[,variate$I,with=F], varinfo, Iout='init'),
                     Imean0 = varinfo[[ 'hmean']][variate$I],
                     Ivar0 = varinfo[[ 'hsd']][variate$I]^2,
                     Ishapeout0 = varinfo[[ 'hshapeout']][variate$I],
                     Ishapein0 = varinfo[[ 'hshapein']][variate$I],
                     Ivarscale0 = varinfo[[ 'hvarscale']][variate$I]^2
                     )},
        if(len$B > 0){list( # binay variate
                     Bshapeout0 = varinfo[[ 'hshapeout']][variate$B],
                     Bshapein0 = varinfo[[ 'hshapein']][variate$B]
                 )},
        if(len$C > 0){list( # categorical variate
                     Calpha0 = t(sapply(variate$B, function(v){
                         c( rep(varinfo[[ 'hshapeout']][v], varinfo[[ 'max']][v]),
                           rep(2^(-40), Cmaxn-varinfo[[ 'max']][v]) )
                     }))
                 )},
        if((!casualinitvalues) && posterior){list(
                                                 W = rep(1/nclusters, nclusters),
                                                 K = rep(1, ndata) # all in one cluster at first
                                            )},
        if(casualinitvalues && posterior){list(
                                             W = rdirch(1, alpha=rep(1,nclusters)),
                                             K = sample(1:nclusters, ndata, replace=TRUE)
                                         )}
    )}


##
#### Mathematical representation of long-run frequency distributions
finitemix <- nimbleCode({
    if(nalpha > 1){# distribution over concentration parameter
        Alphaindex ~ dcat(prob=probalpha0[1:nalpha])
        W[1:nclusters] ~ ddirch(alpha=walpha0[Alphaindex, 1:nclusters])
    }else{
        W[1:nclusters] ~ ddirch(alpha=walpha0[1:nclusters])
    }
    ##
    if(len$R > 0){# real variates
            for(v in 1:Rn){
                Rrate[v] ~ dinvgamma(shape=Rshapein0[v], scale=Rvarscale0[v])
            }
        }
    if(len$O > 0){# one-censored variates
            for(v in 1:On){
                Orate[v] ~ dinvgamma(shape=Oshapein0[v], scale=Ovarscale0[v])
            }
        }
    if(len$D > 0){# doubly-censored variates
            for(v in 1:Dn){
                Drate[v] ~ dinvgamma(shape=Dshapein0[v], scale=Dvarscale0[v])
            }
        }
    if(len$I > 0){# integer variates
            for(v in 1:In){
                Irate[v] ~ dinvgamma(shape=Ishapein0[v], scale=Ivarscale0[v])
            }
        }
    ##
    for(k in 1:nclusters){
        if(len$R > 0){# real variates
            for(v in 1:Rn){
                Rmean[v, k] ~ dnorm(mean=Rmean0[v], var=Rvar0[v])
                Rvar[v, k] ~ dinvgamma(shape=Rshapeout0[v], rate=Rrate[v])
            }
        }
        if(len$O > 0){# logarithmic censored variates
            for(v in 1:On){
                Omean[v, k] ~ dnorm(mean=Omean0[v], var=Ovar0[v])
                Ovar[v, k] ~ dinvgamma(shape=Oshapeout0[v], rate=Orate[v])
            }
        }
        if(len$D > 0){# bounded continuous variates
            for(v in 1:Dn){
                Dmean[v, k] ~ dnorm(mean=Dmean0[v], var=Dvar0[v])
                Dvar[v, k] ~ dinvgamma(shape=Dshapeout0[v], rate=Drate[v])
            }
        }
        if(len$I > 0){# bounded continuous variates
            for(v in 1:In){
                Imean[v, k] ~ dnorm(mean=Imean0[v], var=Ivar0[v])
                Ivar[v, k] ~ dinvgamma(shape=Ishapeout0[v], rate=Irate[v])
            }
        }
        if(len$B > 0){# binary variates
            for(v in 1:Bn){
                Bprob[v, k] ~ dbeta(shape1=Bshapein0[v], shape2=Bshapeout0[v])
            }
        }
        if(len$C > 0){# categorical variates
            for(v in 1:Cn){
                Cprob[v, k, 1:Cmaxn] ~ ddirch(alpha=Calpha0[v, 1:Cmaxn])
            }
        }
    }
    ##
    if(posterior){# cluster occupations
        for(d in 1:ndata){
            K[d] ~ dcat(prob=W[1:nclusters])
            ##
            if(len$R > 0){# real variates
                for(v in 1:Rn){
                    Rdata[d, v] ~ dnorm(mean=Rmean[v, K[d]], var=Rvar[v, K[d]])
                }
            }
            if(len$O > 0){# one-censored variates
                for(v in 1:On){
                    Oaux[d, v] ~ dconstraint(Odata[d, v] >= Oleft[d, v])
                    Odata[d, v] ~ dnorm(mean=Omean[v, K[d]], var=Ovar[v, K[d]])
                }
            }
            if(len$D > 0){# doubly-censored variates
                for(v in 1:Dn){
                    Daux[d, v] ~ dconstraint(Ddata[d, v] >= Dleft[d, v] & Ddata[d, v] <= Dright[d, v])
                    Ddata[d, v] ~ dnorm(mean=Dmean[v, K[d]], var=Dvar[v, K[d]])
                }
            }
            if(len$I > 0){# integer variates
                for(v in 1:In){
                    Iaux[d, v] ~ dconstraint(Icont[d, v] > Ileft[d, v] & Icont[d, v] < Iright[d, v])
                    Icont[d, v] ~ dnorm(mean=Imean[v, K[d]], var=Ivar[v, K[d]])
                }
            }
            if(len$B > 0){# binary variates
                for(v in 1:Bn){
                    Bdata[d, v] ~ dbern(prob=Bprob[v, K[d]])
                }
            }
            if(len$C > 0){# categorical variates
                for(v in 1:Cn){
                    Cdata[d, v] ~ dcat(prob=Cprob[v, K[d], 1:Cmaxn])
                }
            }
        }
    }
})

##
timecount <- Sys.time()
##
finitemixnimble <- nimbleModel(code=finitemix, name='finitemixnimble1',
                               constants=constants,
                               inits=initsFunction(),
                               data=(if(posterior){datapoints}else{list()}),
                               dimensions=c(
                                   list(W=nclusters),
                                   if(nalpha > 1){list(
                                                      walpha0=c(nalpha,nclusters)
                                                  )},
                                   if(nalpha ==1){list(
                                                      walpha0=nclusters
                                                  )},
                                   if(len$R > 0){list(
                                                     Rrate=len$R,
                                                     Rmean=c(len$R,nclusters),
                                                     Rvar=c(len$R,nclusters)
                                                 )},
                                   if(len$O > 0){list(
                                                     Orate=len$O,
                                                     Omean=c(len$O,nclusters),
                                                     Ovar=c(len$O,nclusters)
                                                 )},
                                   if(len$D > 0){list(
                                                     Drate=len$D,
                                                     Dmean=c(len$D,nclusters),
                                                     Dvar=c(len$D,nclusters)
                                                 )},
                                   if(len$I > 0){list(
                                                     Irate=len$I,
                                                     Imean=c(len$I,nclusters),
                                                     Ivar=c(len$I,nclusters)
                                                 )},
                                   if(len$B > 0){list(
                                                     Bprob=c(len$B,nclusters)
                                                 )},
                                   if(len$C > 0){list(
                                                     Cprob=c(len$C,nclusters,Cmaxn)
                                                 )},
                                   if(posterior){c(
                                       list(K=ndata),
                                       if(len$R > 0){list(Rdata=c(ndata,len$R))},
                                       if(len$O > 0){list(
                                                         Odata=c(ndata,len$O),
                                                         Oleft=c(ndata,len$O),
                                                         Oaux=c(ndata,len$O)
                                                     )},
                                       if(len$D > 0){list(
                                                         Ddata=c(ndata,len$D),
                                                         Dleft=c(ndata,len$D),
                                                         Dright=c(ndata,len$D),
                                                         Daux=c(ndata,len$D)
                                                     )},
                                       if(len$I > 0){list(
                                                         Icont=c(ndata,len$I),
                                                         Ileft=c(ndata,len$I),
                                                         Iright=c(ndata,len$I),
                                                         Iaux=c(ndata,len$I)
                                                     )},
                                       if(len$B > 0){list(Bdata=c(ndata,len$B))},
                                       if(len$C > 0){list(Cdata=c(ndata,len$C))}
                                   )}
                               )
                               )
##
Cfinitemixnimble <- compileNimble(finitemixnimble, showCompilerOutput=FALSE)
gc()

confnimble <- configureMCMC(Cfinitemixnimble, #nodes=NULL,
                            monitors=c('W',
                                       if(len$R > 0){c('Rmean', 'Rvar')},
                                       if(len$O > 0){c('Omean', 'Ovar')},
                                       if(len$D > 0){c('Dmean', 'Dvar')},
                                       if(len$I > 0){c('Imean', 'Ivar')},
                                       if(len$B > 0){c('Bprob')},
                                       if(len$C > 0){c('Cprob')}
                                       ),
                            monitors2=c(
                                ## if(posterior && len$R > 0){'Rdata'},
                                ## if(posterior && len$L > 0){'Ldata'},
                                ## if(posterior && len$T > 0){c('Tauxint', 'Tdatacont')},
                                ## if(posterior && len$I > 0){c('Iauxcont', 'Idataint')},
                                ## if(posterior && len$B > 0){'Bdata'},
                                ## if(posterior && len$C > 0){'Cdata'},
                                ## 'Icont',
                                if(nalpha > 1){'Alphaindex'},
                                if(posterior){'K'}
                            )
                            )
## confnimble$printSamplers(executionOrder=TRUE)
##
## takename <- function(x){sub('([^[]+)(.*)','\\1',confnimble$getSamplers(ind=x)[[1]]$target)}
## orde <- confnimble$getSamplerExecutionOrder()
## norde <- sapply(orde,function(i){sub('([^[]+)(.*)','\\1',confnimble$getSamplers(ind=i)[[1]]$target)})
## mysampleorder <- c( 'Rdata', 'Ldata', 'Tauxint', 'Tdatacont', 'Iauxcont', 'Idataint', 'Bdata', 'Cdata', 'K', 'Rrate', 'Lrate', 'Trate', 'Irate', 'Rvar', 'Rmean', 'Lvar', 'Lmean', 'Tvar', 'Tmean', 'Ivar', 'Imean', 'Bprob', 'Cprob', 'W', 'Alphaindex' )
## newsampleorder <- unlist(sapply(mysampleorder, function(v){which(norde == v)}))
## ##
## confnimble$setSamplerExecutionOrder(newsampleorder)
## ##
## ## confnimble$printSamplers(executionOrder=TRUE)
## if(!all(sort(orde)==sort(newsampleorder))){warning('sampler mismatch')}
confnimble$setSamplerExecutionOrder(rev(confnimble$getSamplerExecutionOrder()))

mcsampler <- buildMCMC(confnimble)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

cat('\nSetup time: ')
print(Sys.time() - timecount)

##################################################
## Monte Carlo sampler and plots of MC diagnostics
##################################################
traces <- mcsamples <- NULL
burnin <- 0
continue <- TRUE
stage <- -1
niter <- niter0
thin <- 1
totaliter <- 0
time0 <- Sys.time()
while(continue){
    calctime <- Sys.time()
    stage <- stage+1

    cat(paste0('\n\n==== STAGE ', stage, ' ===='))
    cat(paste0('\nIterations: ',niter,', thinning: ',thin,'\n'))
    gc()
    if(stage==0){# burn-in stage
        set.seed(mcmcseed+stage+100)
        Cfinitemixnimble$setInits(initsFunction())
        newmcsamples <- Cmcsampler$run(niter=niter+1, thin=thin, thin2=niter, nburnin=1, time=T)
        samplertimes <- Cmcsampler$getTimes()
        names(samplertimes) <- sapply(confnimble$getSamplers(),function(x)x$target)
        ## sum(samplertimes[grepl('^Rdata',names(samplertimes))])
        ##
        sprefixes <- unique(sub('^([^[]+)(\\[.*\\])', '\\1', names(samplertimes)))
        cat(paste0('\nSampler times:\n'))
        sort(sapply(sprefixes, function(x)sum(samplertimes[grepl(x,names(samplertimes))])),decreasing=T)
##        newmcsamples <- Cmcsampler$run(niter=1024, thin=1, thin2=1, nburnin=0, time=T)
    ## }else if(is.character(resume)){# continuing previous # must be fixed
    ##     initsc <- readRDS(paste0(dirname,resume))
    ##     inits0 <- initsFunction()
    ##     for(aname in names(inits0)){inits0[[aname]] <- initsc[[aname]]}
    ##     thin <- initsc[['thin']]
    ##     set.seed(mcmcseed+stage+100)
    ##     Cfinitemixnimble$setInits(initsc)
    ##     newmcsamples <- Cmcsampler$run(niter=niter*thin, thin=thin, thin2=niter*thin, nburnin=0)
    }else{# subsequent sampling stages
        cat('\nForecasted computation time: ')
        print(comptime*thin*niter)
        newmcsamples <- Cmcsampler$run(niter=niter*thin, thin=thin, thin2=niter*thin, reset=FALSE, resetMV=TRUE, time=F)
    }
    ##
    totaliter <- totaliter + niter*thin
    newmcsamples <- as.matrix(Cmcsampler$mvSamples)
    cat('\nTime MCMC: ')
    print(Sys.time() - calctime)
    ##
    if(any(is.na(newmcsamples))){cat('\nWARNING: SOME NA OUTPUTS')}
    if(any(!is.finite(newmcsamples))){cat('\nWARNING: SOME INFINITE OUTPUTS')}
    ##
    ## save final state of MCMC chain
    finalstate <- as.matrix(Cmcsampler$mvSamples2)
    finalstate <- c(newmcsamples[nrow(newmcsamples),], finalstate[nrow(finalstate),])
    ##
    ## Check how many "clusters" were occupied. Warns if too many
    occupations <- finalstate[grepl('^K\\[', names(finalstate))]
    usedclusters <- length(unique(occupations))
    if(usedclusters > nclusters-5){cat('\nWARNING: TOO MANY CLUSTERS OCCUPIED')}
    cat(paste0('\nOCCUPIED CLUSTERS: ', usedclusters, ' OF ', nclusters),'\n')
##    saveRDS(finalstate2list(finalstate, realVars=realVars, integerVars=integerVars, categoryVars=categoryVars, binaryVars=binaryVars, compoundgamma=compoundgamma), file=paste0(dirname,'_finalstate-R',basename,'--',mcmcseed,'-',stage,'.rds'))
    ##
    ## SAVE THE PARAMETERS
    ##    parmList <- mcsamples2parmlist(mcsamples, realVars, integerVars, categoryVars, binaryVars)
    ##  saveRDS(parmList,file=paste0(dirname,'_frequencies-R',baseversion,'-V',length(varNames),'-D',ndata,'-K',nclusters,'-I',nrow(parmList$q),'--',mcmcseed,'-',stage,'.rds'))
    ##
    ## Diagnostics
    ## Log-likelihood
    if(posterior){
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
        ll <- colSums(log(samplesFDistribution(Y=data.matrix(data0), X=NULL, mcsamples=newmcsamples, varinfo=varinfo, jacobian=FALSE)), na.rm=T) #- sum(log(invjacobian(data.matrix(data0), varinfo)), na.rm=T)
        if(!posterior && !any(is.finite(ll))){
            ll <- rep(0, length(ll))
        }
        lld <- colSums(log(samplesFDistribution(Y=data.matrix(data0[,..predictands]), X=data.matrix(data0[,..predictors]), mcsamples=newmcsamples, varinfo=varinfo, jacobian=FALSE)), na.rm=T) # - sum(log(invjacobian(data.matrix(data0[,..predictands]), varinfo)), na.rm=T)
        lli <- colSums(log(samplesFDistribution(Y=data.matrix(data0[,..predictors]), X=data.matrix(data0[,..predictands]), mcsamples=newmcsamples, varinfo=varinfo, jacobian=FALSE)), na.rm=T) #- sum(log(invjacobian(data.matrix(data0[,..predictors]), varinfo)), na.rm=T)
        ##
        traces <- rbind(traces,
                        10/log(10)/ndata *
                        cbind(loglikelihood=ll,
                              'mean of direct logprobabilities'=lld,
                              'mean of inverse logprobabilities'=lli)
                        )
        traces2 <- traces[apply(traces,1,function(x){all(is.finite(x))}),]
        saveRDS(traces,file=paste0(dirname,'_mctraces-R',basename,'--',mcmcseed,'-',stage,'.rds'))
        flagll <- nrow(traces) != nrow(traces2)
        ##
        if(nrow(traces2)>=1000){
            funMCSE <- function(x){LaplacesDemon::MCSE(x, method='batch.means')$se}
        }else{
            funMCSE <- function(x){LaplacesDemon::MCSE(x)}
        }
        diagnESS <- LaplacesDemon::ESS(traces2)
        diagnIAT <- apply(traces2, 2, function(x){LaplacesDemon::IAT(x)})
        diagnBMK <- LaplacesDemon::BMK.Diagnostic(traces2[1:(4*trunc(nrow(traces2)/4)),], batches=4)[,1]
        diagnMCSE <- 100*apply(traces2, 2, function(x){funMCSE(x)/sd(x)})
        diagnStat <- apply(traces2, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
        diagnBurn <- apply(traces2, 2, function(x){LaplacesDemon::burnin(matrix(x[1:(10*trunc(length(x)/10))], ncol=1))})
        diagnBurn2 <- proposeburnin(traces2, batches=10)
        diagnThin <- proposethinning(traces2)
        ##
        cat(paste0('\nESSs (',requiredESS,'): ',paste0(round(diagnESS), collapse=', ')))
        cat(paste0('\nIATs: ',paste0(round(diagnIAT), collapse=', ')))
        cat(paste0('\nBMKs: ',paste0(round(diagnBMK,3), collapse=', ')))
        cat(paste0('\nMCSEs: ',paste0(round(diagnMCSE,2), collapse=', ')))
        cat(paste0('\nStationary: ',paste0(diagnStat, collapse=', ')))
        cat(paste0('\nBurn-in I: ',paste0(diagnBurn, collapse=', ')))
        cat(paste0('\nBurn-in II: ',diagnBurn2))
        cat(paste0('\nProposed thinning: ',paste0(diagnThin, collapse=', ')),'\n')
        ##
        #########################################
        #### CHECK IF WE NEED TO SAMPLE MORE ####
        #########################################
        if(stage==0){
            thin <- round(max(diagnThin)*1.5)
            burnin <- niter0-1
            continue <- TRUE
        }else{
            if(min(diagnESS) >= ceiling(requiredESS) &
               #max(diagnMCSE) < 6.27 &
               sum(diagnStat) == 3 &
               diagnBurn2 == 0
               ){
                continue <- FALSE
                burnin <- 0
            }else{
                continue <- TRUE
                ## if(max(diagnThin) > 1){
                ##     thin <- thin*max(diagnThin)
                ##     burnin <- nsamples-1
                ## }else{
                burnin <- min(max(diagnBurn2, minstepincrease), nsamples-1)
            }
        }
        niter <- nsamples - nrow(traces) + burnin
        #########################################
        #### END CHECK                       ####
        #########################################
        mcsamples <- rbind(mcsamples, newmcsamples)
        rm(newmcsamples)
        if(savetempsamples | !continue){
            saveRDS(mcsamples,file=paste0(dirname,'_mcsamples-R',basename,'--',mcmcseed,'-',stage,'.rds'))
        }

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
              ## paste0('burn: ', signif(max(diagnBurn[tracegroups[[agroup]]]),6))
              paste0('burn: ', signif(diagnBurn2,6))
              )
        }
        colpalette <- c(7,2,1)
        names(colpalette) <- colnames(traces)
    ##
    ## Plot various info and traces
        cat('\nPlotting MCMC traces')
        graphics.off()
        pdff(paste0(dirname,'mcmcplottraces-R',basename,'--',mcmcseed,'-',stage),'a4')
    ## Summary stats
        matplot(1:2, type='l', col='white', main=paste0('Stats stage ',stage), axes=FALSE, ann=FALSE)
        legendpositions <- c('topleft','topright','bottomleft','bottomright')
        for(alegend in 1:length(grouplegends)){
            legend(x=legendpositions[alegend], bty='n', cex=1.5,
                   legend=grouplegends[[alegend]] )
        }
        legend(x='center', bty='n', cex=1,
               legend=c(
                   paste0('STAGE ',stage),
                   paste0('Occupied clusters: ', usedclusters, ' of ', nclusters),
                   paste0('LL:  ( ', signif(mean(traces[,1]),3), ' +- ', signif(sd(traces[,1]),3),' ) dHart'),
                   'NOTES:',
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
                              ' | MCSE = ', signif(diagnMCSE[avar], 3),
                              ' | stat: ', diagnStat[avar],
                              ' | burn I: ', diagnBurn[avar],
                              ' | burn II: ', diagnBurn2
                              ),
                  ylab=paste0(avar,'/dHart'), xlab='sample', family=family
                  )
        }
dev.off()
    }
    ##
    ## Samples of marginal frequency distributions
    if(plottempdistributions | !continue){
        if(plotmeans){nfsamples <- totsamples
        }else if(!posterior){nfsamples <- 256
        }else{nfsamples <- 100}
        nfsamples <- min(nfsamples, nrow(mcsamples))
        subsample <- round(seq(1,nfsamples, length.out=64))
        ##
        cat('\nPlotting samples of frequency distributions')
        graphics.off()
        pdff(paste0(dirname,'mcmcdistributions-R',basename,'--',mcmcseed,'-',stage),'a4')
        for(v in unlist(variate)){#cat(avar)
            contvar <- varinfo[['type']][v] %in% c('R','O','D')
            rg <- c(varinfo[['plotmin']][v], varinfo[['plotmax']][v])
            if(contvar){
                Xgrid <- cbind(seq(rg[1], rg[2], length.out=256))
            }else{
                Xgrid <- seq(varinfo[['min']][v], varinfo[['max']][v], length.out=varinfo[['n']][v])
                Xgrid <- cbind(Xgrid[Xgrid >= rg[1] & Xgrid <= rg[2]])
            }
            colnames(Xgrid) <- v
            plotsamples <- samplesFDistribution(Y=Xgrid, mcsamples=mcsamples, varinfo=varinfo, subsamples=round(seq(1,nrow(mcsamples),length.out=nfsamples)), jacobian=TRUE)
            ##
            if(posterior){
                par(mfrow=c(1,1))
                ymax <- tquant(apply(plotsamples[,subsample],2,function(x){tquant(x,31/32)}),31/32, na.rm=T)
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
                    fiven <- fivenum(datum)
                    if(!(varinfo[['type']][v] %in% c('O','D'))){
                        if(contvar){
                            nh <- max(10,round(length(datum)/64))
                        }else{
                            nh <- (varinfo[['max']][v]-varinfo[['min']][v])/(varinfo[['n']][v]-1)
                            nh <- seq(varinfo[['min']][v]-nh/2, varinfo[['max']][v]+nh/2, length.out=varinfo[['n']][v]+1)
                        }
                        histo <- thist(datum, n=nh)
                        histomax <- max(rowMeans(plotsamples))/max(histo$density)
                        tplot(x=histo$breaks, y=histo$density*histomax, col=7, alpha=3/4, border=darkgrey, border.alpha=3/4, lty=1, lwd=1, family=family, ylim=c(0,NA), add=TRUE)
                    }else{ # histogram for censored variate
                        interior <- which(datum > varinfo[['min']][v] & datum < varinfo[['max']][v])
                        histo <- thist(datum[interior], n=max(10,round(length(interior)/64)))
                        interiorgrid <- which(Xgrid > varinfo[['min']][v] & Xgrid < varinfo[['max']][v])
                        histomax <- max(rowMeans(plotsamples)[interiorgrid])/max(histo$density)
                        tplot(x=histo$breaks, y=histo$density*histomax, col=7, alpha=3/4, border=darkgrey, border.alpha=3/4, lty=1, lwd=1, family=family, ylim=c(0,NA), add=TRUE)
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
                    abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'),lwd=4)
                }else if((showdata=='scatter')|(showdata==TRUE & contvar)){
                    datum <- data0[[v]]
                    datum <- datum[!is.na(datum)]
                    diffdatum <- c(apply(cbind(c(0,diff(datum)),c(diff(datum),0)),1,min))/2
                    scatteraxis(side=1, n=NA, alpha='88', ext=8, x=rnorm(length(datum),mean=datum,sd=diffdatum),col=yellow)
                }
            }else{
                par(mfrow=c(8,8),mar = c(0,0,0,0))
                tplot(x=list(Xgrid,Xgrid), y=list(rowMeans(plotsamples),rep(0,length(Xgrid))), type='l', col=c(paste0(palette()[3], 'FF'), '#bbbbbb80'), lty=1, lwd=c(2,1), xlab=NA, ylab=NA, ylim=c(0, NA), family=family,
                      xticks=NA, yticks=NA,
                      mar=c(1,1,1,1))
                abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'))
                text(sum(range(Xgrid))/2, par('usr')[4]*0.9, v)
                ##
                for(aplot in 1:63){
                    tplot(x=list(Xgrid,Xgrid), y=list(plotsamples[,subsample[aplot]], rep(0,length(Xgrid))), type='l', col=c(paste0(palette()[1], 'FF'), '#bbbbbb80'), lty=1, lwd=c(1,1), xlab=NA, ylab=NA, ylim=c(0, NA), family=family,
                          xticks=NA, yticks=NA,
                          mar=c(1,1,1,1))
                    abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'))
                    ## if(aplot==1){ text(Xgrid[1], par('usr')[4]*0.9, variateinfo[variate==avar,type], pos=4)}
                    if(aplot==2){ text(Xgrid[1], par('usr')[4]*0.9, paste(signif(c(rg,diff(rg)),2),collapse=' -- '), pos=4)}
                }
            }
            }
    dev.off()
    }
    ##
    cat('\nTime MCMC+diagnostics: ')
    print(Sys.time() - calctime)
    comptime <- (Sys.time() - time0)/totaliter
    ##
    ## mcsamples <- mcsamples[(burnin+1):nrow(mcsamples),,drop=FALSE]
    ## traces <- traces[(burnin+1):nrow(traces),,drop=FALSE]
    ##
    ## if(savetempsamples | !continue){
    ##         saveRDS(mcsamples,file=paste0(dirname,'_mcsamples-R',basename,'--',mcmcseed,'-',stage,'.rds'))
    ##     }

    if(!continue){
        ## Change name of rds files with final results
        file.rename(from=paste0(dirname,'_mcsamples-R',basename,'--',mcmcseed,'-',stage,'.rds'), to=paste0(dirname,'_mcsamples-R',basename,'--',mcmcseed,'-','F','.rds') )
        file.rename(from=paste0(dirname,'_mctraces-R',basename,'--',mcmcseed,'-',stage,'.rds'), to=paste0(dirname,'_mctraces-R',basename,'--',mcmcseed,'-','F','.rds') )
        ## Change name of pdf files with final plots
        file.rename(from=paste0(dirname,'mcmcplottraces-R',basename,'--',mcmcseed,'-',stage,'.pdf'), to=paste0(dirname,'mcmcplottraces-R',basename,'--',mcmcseed,'-','F','.pdf') )
        file.rename(from=paste0(dirname,'mcmcdistributions-R',basename,'--',mcmcseed,'-',stage,'.pdf'), to=paste0(dirname,'mcmcdistributions-R',basename,'--',mcmcseed,'-','F','.pdf'))
        ##
        cat('\n==== RESULTS SEEM STATIONARY. END ====\n\n')
    }else{
        mcsamples <- mcsamples[(burnin+1):nrow(traces),,drop=F]
        traces <- traces[(burnin+1):nrow(traces),,drop=F]
    }
    ##
    
}

############################################################
## End MCMC
############################################################
plan(sequential)

stop('NONE. End of script')
