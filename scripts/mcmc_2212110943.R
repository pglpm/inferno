## Author: PGL  Porta Mana
## Created: 2022-09-08T17:03:24+0200
## Last-Updated: 2022-12-11T13:44:47+0100
#########################################
e## Inference of exchangeable variates (nonparametric density regression)
## using effectively-infinite mixture of product kernels
## Monte Carlo sampling
#########################################

#### USER INPUTS AND CHOICES ####
baseversion <- '_testnewS1' # *** ## Base name of output directory
datafile <- 'testdataS1.csv'#'ingrid_data_nogds6.csv' #***
## datafile <- 'ingrid_data_nogds6.csv' #***
predictorfile <- 'predictors.csv'
predictandfile <- NULL # 'predictors.csv'
varinfofile <- 'varinfoS1.rds'
requiredESS <- 1024*2/20 # required effective sample size
nsamples <- 8*ceiling((requiredESS*1.5)/8) # number of samples AFTER thinning
ndata <- NULL # set this if you want to use fewer data
shuffledata <- FALSE # useful if subsetting data
posterior <- TRUE # if set to FALSE it samples and plots prior samples
minstepincrease <- 8L
savetempsamples <- TRUE # save temporary MCMC samples
plottempdistributions <- TRUE # plot temporary sampled distributions
showdata <- 'histogram' # 'histogram' 'scatter' FALSE TRUE
plotmeans <- TRUE # plot frequency averages
totsamples <- 100 # number of samples if plotting frequency averages
##
niter0 <- 1024L * 1L # 3L # iterations burn-in
nclusters <- 64L
alpha0 <- c(0.5, 1, 2)
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
source(paste0(origdir,'functionsmcmc_2212110943.R')) # load functions for post-MCMC calculations
## varinfo <- data.matrix(read.csv(paste0(origdir,varinfofile), row.names=1))
varinfo <- readRDS(paste0(origdir,varinfofile))

variate <- lapply(variatetypes, function(x)rownames(varinfo)[varinfo[,'type']==x])
len <- lapply(variate,length)

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
    ## logarithmic
    if(len$L > 0){list( Ldata = transf(data0[,variate$L,with=F], varinfo) )},
    ## logarithmic censored
    if(len$S > 0){list(
                   Sdata = transf(data0[,variate$S,with=F], varinfo, Sout='data'),
                   Sleft = transf(data0[,variate$S,with=F], varinfo, Sout='left'),
                   Saux = matrix(1L, nrow=ndata, ncol=len$S)
                  )},
    ## Two-bounded
    if(len$T > 0){list(
                   Tdata = transf(data0[,variate$T,with=F], varinfo, Tout='data'),
                   Tleft = transf(data0[,variate$T,with=F], varinfo, Tout='left'),
                   Tright = transf(data0[,variate$T,with=F], varinfo, Tout='right'),
                   Taux = matrix(1L, nrow=ndata, ncol=len$T)
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
    if(len$S > 0){ list(Sn = len$S) },
    if(len$L > 0){ list(Ln = len$L) },
    if(len$T > 0){ list(Tn = len$T) },
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
                           Alpha = 2,
                           probalpha0 = rep(1/nalpha, nalpha),
                           walpha0 = matrix(alpha0/nclusters,
                                            nrow=nalpha, ncol=nclusters)
                       )}else{list(
                                  walpha0 = rep(1/nclusters, nclusters)
                              )},
        if(len$R > 0){list( # real variate
                     Rmean0 = varinfo[variate$R, 'hmean'],
                     Rvar0 = varinfo[variate$R, 'hsd']^2,
                     Rshapeout0 = varinfo[variate$R, 'hshapeout'],
                     Rshapein0 = varinfo[variate$R, 'hshapein'],
                     Rvarscale0 = varinfo[variate$R, 'hvarscale']^2
                 )},
        if(len$L > 0){list( # logarithmic variate
                     Lmean0 = varinfo[variate$L, 'hmean'],
                     Lvar0 = varinfo[variate$L, 'hsd']^2,
                     Lshapeout0 = varinfo[variate$L, 'hshapeout'],
                     Lshapein0 = varinfo[variate$L, 'hshapein'],
                     Lvarscale0 = varinfo[variate$L, 'hvarscale']^2
                 )},
        if(len$S > 0){list( # doubly-bounded variate
                     Sdata = transf(data0[,variate$S,with=F], varinfo, Sout='init'),
                     Smean0 = varinfo[variate$S, 'hmean'],
                     Svar0 = varinfo[variate$S, 'hsd']^2,
                     Sshapeout0 = varinfo[variate$S, 'hshapeout'],
                     Sshapein0 = varinfo[variate$S, 'hshapein'],
                     Svarscale0 = varinfo[variate$S, 'hvarscale']^2
                 )},
        if(len$T > 0){list( # doubly-bounded variate
                     Tdata = transf(data0[,variate$T,with=F], varinfo, Tout='init'),
                     Tmean0 = varinfo[variate$T, 'hmean'],
                     Tvar0 = varinfo[variate$T, 'hsd']^2,
                     Tshapeout0 = varinfo[variate$T, 'hshapeout'],
                     Tshapein0 = varinfo[variate$T, 'hshapein'],
                     Tvarscale0 = varinfo[variate$T, 'hvarscale']^2
                 )},
        if(len$I > 0){list( # integer ordinal variate
                     Icont = transf(data0[,variate$I,with=F], varinfo, Iout='init'),
                     Imean0 = varinfo[variate$I, 'hmean'],
                     Ivar0 = varinfo[variate$I, 'hsd']^2,
                     Ishapeout0 = varinfo[variate$I, 'hshapeout'],
                     Ishapein0 = varinfo[variate$I, 'hshapein'],
                     Ivarscale0 = varinfo[variate$I, 'hvarscale']^2
                     )},
        if(len$B > 0){list( # binay variate
                     Bshapeout0 = varinfo[variate$B, 'hshapeout'],
                     Bshapein0 = varinfo[variate$B, 'hshapein']
                 )},
        if(len$C > 0){list( # categorical variate
                     Calpha0 = t(sapply(variate$B, function(v){
                         c( rep(varinfo[v, 'hshapeout'], varinfo[v, 'max']),
                           rep(2^(-40), Cmaxn-varinfo[v, 'max']) )
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
        Alpha ~ dcat(prob=probalpha0[1:nalpha])
        W[1:nclusters] ~ ddirch(alpha=walpha0[Alpha, 1:nclusters])
    }else{
        W[1:nclusters] ~ ddirch(alpha=walpha0[1:nclusters])
    }
    ##
    if(len$R > 0){# real variates
            for(v in 1:Rn){
                Rrate[v] ~ dinvgamma(shape=Rshapein0[v], scale=Rvarscale0[v])
            }
        }
    if(len$L > 0){# logarithmic variates
            for(v in 1:Ln){
                Lrate[v] ~ dinvgamma(shape=Lshapein0[v], scale=Lvarscale0[v])
            }
        }
    if(len$S > 0){# logarithmic censored variates
            for(v in 1:Sn){
                Srate[v] ~ dinvgamma(shape=Sshapein0[v], scale=Svarscale0[v])
            }
        }
    if(len$T > 0){# bounded continuous variates
            for(v in 1:Tn){
                Trate[v] ~ dinvgamma(shape=Tshapein0[v], scale=Tvarscale0[v])
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
        if(len$L > 0){# logarithmic variates
            for(v in 1:Ln){
                Lmean[v, k] ~ dnorm(mean=Lmean0[v], var=Lvar0[v])
                Lvar[v, k] ~ dinvgamma(shape=Lshapeout0[v], rate=Lrate[v])
            }
        }
        if(len$S > 0){# logarithmic censored variates
            for(v in 1:Sn){
                Smean[v, k] ~ dnorm(mean=Smean0[v], var=Svar0[v])
                Svar[v, k] ~ dinvgamma(shape=Sshapeout0[v], rate=Srate[v])
            }
        }
        if(len$T > 0){# bounded continuous variates
            for(v in 1:Tn){
                Tmean[v, k] ~ dnorm(mean=Tmean0[v], var=Tvar0[v])
                Tvar[v, k] ~ dinvgamma(shape=Tshapeout0[v], rate=Trate[v])
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
            if(len$L > 0){# logarithmic variates
                for(v in 1:Ln){
                    Ldata[d, v] ~ dnorm(mean=Lmean[v, K[d]], var=Lvar[v, K[d]])
                }
            }
            if(len$S > 0){# logarithmic censored variates
                for(v in 1:Sn){
                    Saux[d, v] ~ dconstraint(Sdata[d, v] >= Sleft[d, v])
                    Sdata[d, v] ~ dnorm(mean=Smean[v, K[d]], var=Svar[v, K[d]])
                }
            }
            if(len$T > 0){# bounded continuous variates
                for(v in 1:Tn){
                    Taux[d, v] ~ dconstraint(Tdata[d, v] >= Tleft[d, v] & Tdata[d, v] <= Tright[d, v])
                    Tdata[d, v] ~ dnorm(mean=Tmean[v, K[d]], var=Tvar[v, K[d]])
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
                                   if(len$L > 0){list(
                                                     Lrate=len$L,
                                                     Lmean=c(len$L,nclusters),
                                                     Lvar=c(len$L,nclusters)
                                                 )},
                                   if(len$S > 0){list(
                                                     Srate=len$S,
                                                     Smean=c(len$S,nclusters),
                                                     Svar=c(len$S,nclusters)
                                                 )},
                                   if(len$T > 0){list(
                                                     Trate=len$T,
                                                     Tmean=c(len$T,nclusters),
                                                     Tvar=c(len$T,nclusters)
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
                                       if(len$L > 0){list(Ldata=c(ndata,len$L))},
                                       if(len$S > 0){list(
                                                         Sdata=c(ndata,len$S),
                                                         Sleft=c(ndata,len$S),
                                                         Saux=c(ndata,len$S)
                                                     )},
                                       if(len$T > 0){list(
                                                         Tdata=c(ndata,len$T),
                                                         Tleft=c(ndata,len$T),
                                                         Tright=c(ndata,len$T),
                                                         Taux=c(ndata,len$T)
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
                                       if(len$L > 0){c('Lmean', 'Lvar')},
                                       if(len$S > 0){c('Smean', 'Svar')},
                                       if(len$T > 0){c('Tmean', 'Tvar')},
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
                                'Sdata',
                                if(nalpha > 1){'Alpha'},
                                if(posterior){'K'}
                            )
                            )
## confnimble$printSamplers(executionOrder=TRUE)
##
## takename <- function(x){sub('([^[]+)(.*)','\\1',confnimble$getSamplers(ind=x)[[1]]$target)}
## orde <- confnimble$getSamplerExecutionOrder()
## norde <- sapply(orde,function(i){sub('([^[]+)(.*)','\\1',confnimble$getSamplers(ind=i)[[1]]$target)})
## mysampleorder <- c( 'Rdata', 'Ldata', 'Tauxint', 'Tdatacont', 'Iauxcont', 'Idataint', 'Bdata', 'Cdata', 'K', 'Rrate', 'Lrate', 'Trate', 'Irate', 'Rvar', 'Rmean', 'Lvar', 'Lmean', 'Tvar', 'Tmean', 'Ivar', 'Imean', 'Bprob', 'Cprob', 'W', 'Alpha' )
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
        ## newmcsamples <- Cmcsampler$run(niter=niter+1, thin=thin, thin2=niter, nburnin=1, time=T)
        newmcsamples <- Cmcsampler$run(niter=1024, thin=1, thin2=1, nburnin=0, time=T)
    }else if(is.character(resume)){# continuing previous # must be fixed
        initsc <- readRDS(paste0(dirname,resume))
        inits0 <- initsFunction()
        for(aname in names(inits0)){inits0[[aname]] <- initsc[[aname]]}
        thin <- initsc[['thin']]
        set.seed(mcmcseed+stage+100)
        Cfinitemixnimble$setInits(initsc)
        newmcsamples <- Cmcsampler$run(niter=niter*thin, thin=thin, thin2=niter*thin, nburnin=0)
    }else{# subsequent sampling stages
        cat('\nForecasted computation time: ')
        print(comptime*thin*niter)
        newmcsamples <- Cmcsampler$run(niter=niter*thin, thin=thin, thin2=niter*thin, reset=FALSE, resetMV=TRUE)
    }
    ##
    totaliter <- totaliter + niter*thin

    newmcsamples <- as.matrix(Cmcsampler$mvSamples)
    newmcsamples2 <- as.matrix(Cmcsampler$mvSamples2)

    pdff('debugplotsS1')
        v <- colnames(newmcsamples)[grepl('Smean',colnames(newmcsamples))]
        tplot(y=newmcsamples[,v],ylab='Smeans',xlab=NA,lty=1,alpha=0.25)
        v <- colnames(newmcsamples)[grepl('Svar',colnames(newmcsamples))]
        tplot(y=log10(newmcsamples[,v])/2,ylab='Ssds',xlab=NA,lty=1,alpha=0.25)
    dev.off()




        ## log-censored
    v <- 'TRAASCOR_neuro'
    vn <- which(variate$S == v)
    xgrid <- cbind(seq(varinfo[v,'datamin'],varinfo[v,'max'],length.out=128))
    colnames(xgrid) <- v
    extremes <- c(length(xgrid))
    txgrid <- transf(xgrid,varinfo,Sout='')
    subsamples <- 1:nrow(newmcsamples)
testt <- sapply(subsamples, function(sam){
    rowSums(sapply(1:nclusters, function(clu){
        dnorm(x=txgrid[-extremes],
              mean=newmcsamples[sam,paste0('Smean[',vn,', ',clu,']')],
              sd=sqrt(newmcsamples[sam,paste0('Svar[',vn,', ',clu,']')])
              ) *
            newmcsamples[sam,paste0('W[',clu,']')]
    }))
})
    testte <- sapply(subsamples, function(sam){
    sum(sapply(1:nclusters, function(clu){
        pnorm(q=txgrid[extremes],
              mean=newmcsamples[sam,paste0('Smean[',vn,', ',clu,']')],
              sd=sqrt(newmcsamples[sam,paste0('Svar[',vn,', ',clu,']')]),
              lower.tail=F
              ) * newmcsamples[sam,paste0('W[',clu,']')]
    }))
})

    testrt <- samplesFDistribution(Y=xgrid,X=NULL,mcsamples=newmcsamples,varinfo=varinfo, subsamples=subsamples, jacobian=F)
    
    tplot(x=txgrid[extremes],y=rowMeans(testrt)[extremes]*max(rowMeans(testrt)[-extremes]),ylim=c(0, max(rowMeans(testrt)[-extremes])),lty=1,col=1,type='p')
    tplot(x=txgrid[-extremes],y=rowMeans(testrt)[-extremes],ylim=c(0,NA),lty=1,col=1,add=T)
    
        tplot(x=txgrid[extremes],y=mean(testte)*max(rowMeans(testrt)[-extremes]),ylim=c(0,NA),lty=1,col=2,type='p',pch=2,add=T)
tplot(x=txgrid[-extremes],y=rowMeans(testt),ylim=c(0,NA),lty=2,col=2,add=T)

    histoT <- thist(transf(data.matrix(data0[c(data0[,..v]<varinfo[v,'max']),v,with=F]),varinfo,Sout=''),n=100)
    histoTe <- sum(c(data0[,..v] >= varinfo[v,'max']))/nrow(data0)
    ##
    tplot(x=histoT$breaks,y=histoT$density/max(histoT$density)*max(rowMeans(testrt)[-extremes]),add=T,col=3)
    tplot(x=txgrid[extremes],y=histoTe,type='p',pch=3,col=3, add=T)

        ## log-censored, with Jacobian
    v <- 'TRAASCOR_neuro'
    vn <- which(variate$S == v)
    xgrid <- cbind(seq(varinfo[v,'datamin'],varinfo[v,'max'],length.out=128))
    colnames(xgrid) <- v
    extremes <- c(length(xgrid))
    txgrid <- transf(xgrid,varinfo,Sout='')
    subsamples <- 1:nrow(newmcsamples)
testt <- sapply(subsamples, function(sam){
    rowSums(sapply(1:nclusters, function(clu){
        dnorm(x=txgrid[-extremes],
              mean=newmcsamples[sam,paste0('Smean[',vn,', ',clu,']')],
              sd=sqrt(newmcsamples[sam,paste0('Svar[',vn,', ',clu,']')])
              ) *
            newmcsamples[sam,paste0('W[',clu,']')]
    }))
})/(xgrid[-extremes]*varinfo[v,'scale'])
    testte <- sapply(subsamples, function(sam){
    sum(sapply(1:nclusters, function(clu){
        pnorm(q=txgrid[extremes],
              mean=newmcsamples[sam,paste0('Smean[',vn,', ',clu,']')],
              sd=sqrt(newmcsamples[sam,paste0('Svar[',vn,', ',clu,']')]),
              lower.tail=F
              ) * newmcsamples[sam,paste0('W[',clu,']')]
    }))
})

    testrt <- samplesFDistribution(Y=xgrid,X=NULL,mcsamples=newmcsamples,varinfo=varinfo, subsamples=subsamples, jacobian=T)
    
    tplot(x=xgrid[-extremes],y=rowMeans(testrt)[-extremes],ylim=c(0,NA),lty=1,col=1,xlim=range(xgrid))
    tplot(x=xgrid[extremes],y=rowMeans(testrt)[extremes]*max(rowMeans(testrt)[-extremes]),ylim=c(0, max(rowMeans(testrt)[-extremes])),lty=1,col=1,type='p',add=T)
    
    tplot(x=xgrid[-extremes],y=rowMeans(testt),ylim=c(0,NA),lty=2,col=2,add=T)
    tplot(x=xgrid[extremes],y=mean(testte)*max(rowMeans(testrt)[-extremes]),ylim=c(0,NA),lty=1,col=2,type='p',pch=2,add=T)

    histoT <- thist(data.matrix(data0[c(data0[,..v]<varinfo[v,'max']),v,with=F]),n=100)
    histoTe <- sum(c(data0[,..v] >= varinfo[v,'max']))/nrow(data0)
    ##
    tplot(x=histoT$breaks,y=histoT$density/max(histoT$density)*max(rowMeans(testrt)[-extremes]),add=T,col=3)
    tplot(x=xgrid[extremes],y=histoTe*max(rowMeans(testrt)[-extremes]),type='p',pch=3,col=3, add=T)




    ## real
    v <- 'AGE'
    vn <- which(variate$L == v)
    xgrid <- cbind(seq(varinfo[v,'plotmin'],varinfo[v,'plotmax'],length.out=128))
    colnames(xgrid) <- v
    txgrid <- transf(xgrid,varinfo,Tout='')
    subsamples <- 1:nrow(newmcsamples)
testt <- sapply(subsamples, function(sam){
    rowSums(sapply(1:nclusters, function(clu){
        dnorm(x=txgrid,
              mean=newmcsamples[sam,paste0('Lmean[',vn,', ',clu,']')],
              sd=sqrt(newmcsamples[sam,paste0('Lvar[',vn,', ',clu,']')])
              ) *
            newmcsamples[sam,paste0('W[',clu,']')]
    }))
})

    testrt <- testsamplesFDistribution(Y=xgrid,X=NULL,mcsamples=newmcsamples,varinfo=varinfo, subsamples=subsamples, jacobian=F)
    
    tplot(x=txgrid,y=rowMeans(testrt),ylim=c(0,NA),lty=1,col=1)
    
tplot(x=txgrid,y=rowMeans(testt),ylim=c(0,NA),lty=2,col=2,add=T)

    histoT <- thist(transf(data.matrix(data0[,v,with=F]),varinfo,Tout=''),n=50)
tplot(x=histoT$breaks,y=histoT$density/max(histoT$density)*max(rowMeans(testrt)[-extremes]),add=T)

    ## doubly-bounded
    v <- 'TRABSCOR_neuro'
    vn <- which(variate$T == v)
    xgrid <- cbind(seq(varinfo[v,'min'],varinfo[v,'max'],length.out=128))
    colnames(xgrid) <- v
    extremes <- c(1,length(xgrid))
    txgrid <- transf(xgrid,varinfo,Tout='')
    subsamples <- 1:nrow(newmcsamples)
testt <- sapply(subsamples, function(sam){
    rowSums(sapply(1:nclusters, function(clu){
        dnorm(x=txgrid[-extremes],
              mean=newmcsamples[sam,paste0('Tmean[',vn,', ',clu,']')],
              sd=sqrt(newmcsamples[sam,paste0('Tvar[',vn,', ',clu,']')])
              ) *
            newmcsamples[sam,paste0('W[',clu,']')]
    }))
})
testte <- sapply(subsamples, function(sam){
    rowSums(sapply(1:nclusters, function(clu){
        c(pnorm(q=txgrid[extremes[1]],
              mean=newmcsamples[sam,paste0('Tmean[',vn,', ',clu,']')],
              sd=sqrt(newmcsamples[sam,paste0('Tvar[',vn,', ',clu,']')])
              ),
          pnorm(q=txgrid[extremes[2]],
              mean=newmcsamples[sam,paste0('Tmean[',vn,', ',clu,']')],
              sd=sqrt(newmcsamples[sam,paste0('Tvar[',vn,', ',clu,']')]),
              lower.tail=F
              )
          ) * newmcsamples[sam,paste0('W[',clu,']')]
    }))
})

    testrt <- testsamplesFDistribution(Y=xgrid,X=NULL,mcsamples=newmcsamples,varinfo=varinfo, subsamples=subsamples, jacobian=F)
    
    tplot(x=txgrid[extremes],y=rowMeans(testrt)[extremes]*max(rowMeans(testrt)[-extremes]),ylim=c(0, max(rowMeans(testrt)[-extremes])),lty=1,col=1,type='p')
    tplot(x=txgrid[-extremes],y=rowMeans(testrt)[-extremes],ylim=c(0,NA),lty=1,col=1,add=T)
    
        tplot(x=txgrid[extremes],y=rowMeans(testte)*max(rowMeans(testrt)[-extremes]),ylim=c(0,NA),lty=1,col=2,type='p',pch=2,add=T)
tplot(x=txgrid[-extremes],y=rowMeans(testt),ylim=c(0,NA),lty=2,col=2,add=T)

    histoT <- thist(transf(data.matrix(data0[,v,with=F]),varinfo,Tout=''),n=100)
tplot(x=histoT$breaks,y=histoT$density/max(histoT$density)*max(rowMeans(testrt)[-extremes]),add=T)



    ## integer
    v <- 'RAVLT_immediate'
    vn <- which(variate$I == v)
    xgrid <- cbind(seq(varinfo[v,'min'],varinfo[v,'max'],length.out=varinfo[v,'n']))#[120:122,,drop=F]
    colnames(xgrid) <- v
    txgrid <- transf(xgrid,varinfo,Iout='init')
    lxgrid <- transf(xgrid,varinfo,Iout='left')
    rxgrid <- transf(xgrid,varinfo,Iout='right')
    subsamples <- 1:nrow(newmcsamples)
testt <- sapply(subsamples, function(sam){
    rowSums(sapply(1:nclusters, function(clu){
        ( pnorm(q=rxgrid,
              mean=newmcsamples[sam,paste0('Imean[',vn,', ',clu,']')],
              sd=sqrt(newmcsamples[sam,paste0('Ivar[',vn,', ',clu,']')])
              ) -
         pnorm(q=lxgrid,
              mean=newmcsamples[sam,paste0('Imean[',vn,', ',clu,']')],
              sd=sqrt(newmcsamples[sam,paste0('Ivar[',vn,', ',clu,']')])
              )) *
            newmcsamples[sam,paste0('W[',clu,']')]
    }))
})

    testrt <- testsamplesFDistribution(Y=xgrid,X=NULL,mcsamples=newmcsamples,varinfo=varinfo, subsamples=subsamples, jacobian=F)
    
tplot(x=txgrid,y=rowMeans(testrt),ylim=c(0,NA),lty=1,col=1)
    
tplot(x=txgrid,y=rowMeans(testt),ylim=c(0,NA),lty=2,col=2,add=T)

    histoT <- thist(transf(data.matrix(data0[,v,with=F]),varinfo,Iout='init'),n=c(txgrid[1],lxgrid[-1],txgrid[length(txgrid)]+0.01))
tplot(x=histoT$breaks,y=histoT$counts/sum(histoT$counts),add=T)






    pdff('debugplotsIS2')
    for(v in colnames(newmcsamples2)){
        tplot(y=newmcsamples2[,v],ylab=v,xlab=paste0(signif(range(newmcsamples2[,v]),3),collapse=' , '))
    }
    dev.off()
    

    pdff('debugplotsIS2')
    for(v in colnames(newmcsamples2)){
        tplot(y=newmcsamples2[,v],ylab=v,xlab=paste0(signif(range(newmcsamples2[,v]),3),collapse=' , '))
    }
    for(v in colnames(newmcsamples)){
        if(grepl('mean',v)){
        tplot(y=newmcsamples[,v],ylab=v,xlab=paste0(signif(range(newmcsamples[,v]),3),collapse=' , '))
    }}
        v <- colnames(newmcsamples)[grepl('Tmean',colnames(newmcsamples))]
        tplot(y=newmcsamples[,v],ylab='Tmeans',xlab=NA,lty=1,alpha=0.25)
        v <- colnames(newmcsamples)[grepl('Imean',colnames(newmcsamples))]
        tplot(y=newmcsamples[,v],ylab='Imeans',xlab=NA,lty=1,alpha=0.25)
        v <- colnames(newmcsamples)[grepl('Tvar',colnames(newmcsamples))]
        tplot(y=log10(newmcsamples[,v])/2,ylab='Tsds',xlab=NA,lty=1,alpha=0.25)
        v <- colnames(newmcsamples)[grepl('Ivar',colnames(newmcsamples))]
        tplot(y=log10(newmcsamples[,v])/2,ylab='Isds',xlab=NA,lty=1,alpha=0.25)
    dev.off()
    
    cat('\nTime MCMC: ')
    print(Sys.time() - calctime)
    ##
    ## ## Check sample-time partition
    ## times <- Cmcsampler$getTimes()
    ## names(times) <- sapply(confnimble$getSamplers(),function(x)x$target)
    ## ##
    ## cbind(sort(times[c('C[1]','q[1:64]','meanR[1, 1]', 'varR[1, 1]', 'probB[1, 1]', 'probC[1, 1, 1:21]', 'varRrate[1]')]))
    ## ##
    ## test <- sapply(c('C','q','meanR', 'varR', 'probB', 'probC', 'varRrate'),
    ##        function(x){
    ##            sum(times[grep(paste0('^',x),names(times))])
    ##        })
    ## names(test) <- c('C','q','meanR', 'varR', 'probB', 'probC', 'varRrate')
    ## cbind(sort(test))
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
        ll <- colSums(log(samplesFDistribution(Y=data.matrix(data0), X=NULL, mcsamples=newmcsamples, varinfo=varinfo, jacobian=FALSE))) + sum(log(invjacobian(data.matrix(data0), varinfo)), na.rm=T)
        flagll <- FALSE
        if(!posterior && !any(is.finite(ll))){
            flagll <- TRUE
            ll <- rep(0, length(ll))
        }
        lld <- colSums(log(samplesFDistribution(Y=data.matrix(data0[,..predictands]), X=data.matrix(data0[,..predictors]), mcsamples=newmcsamples, varinfo=varinfo, jacobian=FALSE))) + sum(log(invjacobian(data.matrix(data0[,..predictands]), varinfo)), na.rm=T)
        lli <- colSums(log(samplesFDistribution(Y=data.matrix(data0[,..predictors]), X=data.matrix(data0[,..predictands]), mcsamples=newmcsamples, varinfo=varinfo, jacobian=FALSE))) + sum(log(invjacobian(data.matrix(data0[,..predictors]), varinfo)), na.rm=T)
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
        }else{nfsamples <- 63}
        nfsamples <- min(nfsamples, nrow(mcsamples))
        subsample <- round(seq(1,nfsamples, length.out=63))
        ##
        cat('\nPlotting samples of frequency distributions')
        graphics.off()
        pdff(paste0(dirname,'mcmcdistributions-R',basename,'--',mcmcseed,'-',stage),'a4')
        for(v in unlist(variate)){#cat(avar)
            contvar <- varinfo[v,'type'] %in% variatetypes[c('R','L','T')]
            rg <- varinfo[v,c('plotmin','plotmax')]
            if(contvar){
                Xgrid <- cbind(seq(rg[1], rg[2], length.out=256))
            }else{
                Xgrid <- seq(varinfo[v,'min'], varinfo[v,'max'], length.out=varinfo[v,'n'])
                Xgrid <- cbind(Xgrid[Xgrid >= rg[1] & Xgrid <= rg[2]])
            }
            colnames(Xgrid) <- v
            plotsamples <- samplesFDistribution(Y=Xgrid, mcsamples=mcsamples, varinfo=varinfo, subsamples=round(seq(1,nrow(mcsamples),length.out=nfsamples)), jacobian=TRUE)
            ##
            if(posterior){
                par(mfrow=c(1,1))
                ymax <- tquant(apply(plotsamples[,subsample],2,function(x){tquant(x,99/100)}),99/100, na.rm=T)
                if(varinfo[v,'type'] != variatetypes['T']){
                    tplot(x=Xgrid, y=plotsamples[,subsample], type='l', col=paste0(palette()[5], '44'), lty=1, lwd=2, xlab=paste0(v), ylab=paste0('frequency'), ylim=c(0, ymax), family=family)
                    if(plotmeans){
                        tplot(x=Xgrid, y=rowMeans(plotsamples, na.rm=T), type='l', col=paste0(palette()[1], '88'), lty=1, lwd=4, add=T)
                    }
                }else{ # plot of a continuous doubly-bounded variate
                    interior <- which(Xgrid > varinfo[v,'min'] & Xgrid < varinfo[v,'max'])
                    tplot(x=Xgrid[interior], y=plotsamples[interior,subsample], type='l', col=paste0(palette()[5], '44'), lty=1, lwd=2, xlab=paste0(v), ylab=paste0('frequency'), ylim=c(0, ymax), family=family)
                    if(length(interior) < length(Xgrid)){
                        tplot(x=Xgrid[-interior], y=plotsamples[-interior,subsample,drop=F]*ymax, type='p', pch=1, cex=1, col=paste0(palette()[5], '44'), lty=1, lwd=2, xlab=paste0(v), ylab=paste0('frequency'), ylim=c(0, ymax), family=family,add=T)
                        }
                    if(plotmeans){
                        tplot(x=Xgrid[interior], y=rowMeans(plotsamples, na.rm=T)[interior], type='l', col=paste0(palette()[1], '88'), lty=1, lwd=4, add=T)
                    if(length(interior) < length(Xgrid)){
                        tplot(x=Xgrid[-interior], y=rowMeans(plotsamples, na.rm=T)[-interior]*ymax, type='p', pch=2, cex=1, col=paste0(palette()[1], '88'), lty=1, lwd=3, add=T)
                        }
                    }
                }
                ##
                if((showdata=='histogram')||(showdata==TRUE && !contvar)){
                    datum <- data0[[v]]
                    datum <- datum[!is.na(datum)]
                    ## fiven <- varinfo[v,c('min','Q1','Q2','Q3','max')]
                    fiven <- fivenum(datum)
                    if(varinfo[v,'type'] != variatetypes['T']){
                        ## histo <- thist(datum, n=(if(contvar){min(max(10,sqrt(ndata)),100)}else{'i'}))#-exp(mean(log(c(round(sqrt(length(datum))), length(Xgrid))))))
                        histo <- thist(datum, n=round(ndata/64))#-exp(mean(log(c(round(sqrt(length(datum))), length(Xgrid))))))
                        histomax <- max(rowMeans(plotsamples))/max(histo$density)
                        tplot(x=histo$breaks, y=histo$density*histomax, col=grey, alpha=0.75, border=NA, lty=1, lwd=1, family=family, ylim=c(0,NA), add=TRUE)
                    }else{ # histogram for a continuous doubly-bounded variate
                        interior <- which(datum > varinfo[v,'min'] & datum < varinfo[v,'max'])
                        histo <- thist(datum[interior], n=round(ndata/64)) #(if(contvar){min(max(10,sqrt(ndata)),100)}else{'i'}))#-exp(mean(log(c(round(sqrt(length(datum))), length(Xgrid))))))
                        interior2 <- which(Xgrid > varinfo[v,'min'] & Xgrid < varinfo[v,'max'])

                        histomax <- max(rowMeans(plotsamples)[interior2])/max(histo$density) # (if(contvar){max(rowMeans(plotsamples)[interior])/max(histo$density)}else{1L})
                        tplot(x=histo$breaks, y=histo$density*histomax, col=grey, alpha=0.75, border=NA, lty=1, lwd=1, family=family, ylim=c(0,NA), add=TRUE)
                        ##
                        pborder <- sum(datum <= varinfo[v,'min'])/length(datum)
                        if(pborder > 0){
                            tplot(x=rep(varinfo[v,'min'],2), y=c(0,pborder*ymax), type='l', col=grey, alpha=0.75, lty=1, lwd=1.5, family=family, ylim=c(0,NA), add=TRUE)
                        }
                        ##
                        pborder <- sum(datum >= varinfo[v,'max'])/length(datum)
                        if(pborder > 0){
                            tplot(x=rep(varinfo[v,'max'],2), y=c(0,pborder*ymax), type='l', col=grey, alpha=0.75, lty=1, lwd=1.5, family=family, ylim=c(0,NA), add=TRUE)
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

alldata <- fread('data_ep.csv')
metadata <- fread('metadata.csv')
allnames <- names(alldata)
metadatanew <- metadata
metadatanew$precision <- as.integer(metadatanew$precision)
##
boundarycat <- 0
boundaryrea <- 0
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


for(i in allnames){if(metadata[variate==i,type]=='integer' & !(min(diff(sort(unique(alldata[[i]]))))==1)){cat(c(i,min(diff(sort(unique(alldata[[i]]))))))}}



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

## q         0.029227
## varRrate  0.707669
## probB     1.131390
## C         6.549711
## meanR    11.425076
## varR     15.414792
## probC    94.969083




xgrid <- seq(-5,-3,length.out=512)
fn <- function(x,dx){pnorm(x)/(2*dx)-1}
testy <- fn(xgrid,1e-4)
tplot(x=xgrid,y=((testy)))

qnorm(2*(1e-4))
## qnorm of twice precision of variate = boundary


mm <- c(0,2,4)
ss <- c(1,2,0.2)
qq <- c(7,4,0.4)
##
f <- function(x){sapply(x,function(y)sum(qq*dnorm(y,mean=mm,sd=ss))/sum(qq))}
##
xgrid <- seq(-3,6,length.out=256)
graphics.off()
pdff('testcog1')
tplot(xgrid,f(xgrid),ylim=c(0,NA),xlabels=NA,xlab='cog score #4',ylab='probability',ylabels=NA,lwd=5,ly=2,yticks=NA)
plotquantiles(xgrid,rbind(
                        f(xgrid)*
                        (1+3*(dnorm(xgrid-min(xgrid))+0.2*dnorm(xgrid-2.5,sd=1)+dnorm(max(xgrid)-xgrid))),
                        f(xgrid)*
                        (1-3*(dnorm(xgrid-min(xgrid))+0.2*dnorm(xgrid-2.5,sd=1)+dnorm(max(xgrid)-xgrid)))
                    ))
dev.off()
me <- sum(xgrid*f(xgrid))/sum(f(xgrid))
sqrt(sum((xgrid-me)^2*f(xgrid))/sum(f(xgrid)))
##
##
mm <- c(0,2,3)
ss <- c(1,2,0.2)*0.5
qq <- c(7,2,0.1)
##
f <- function(x){sapply(x,function(y)sum(qq*dnorm(y,mean=mm,sd=ss))/sum(qq))}
##
xgrid <- seq(-3,6,length.out=256)
graphics.off()
pdff('testcog2')
tplot(xgrid,f(xgrid),ylim=c(0,NA),xlabels=NA,xlab='cog score #4',ylab='probability',ylabels=NA,lwd=5,ly=2,yticks=NA)
plotquantiles(xgrid,rbind(
                        f(xgrid)*
                        (1+7*(dnorm(xgrid-min(xgrid))+0.1*dnorm(xgrid-2.5,sd=1)+dnorm(max(xgrid)-xgrid))),
                        f(xgrid)*
                        (1-7*(dnorm(xgrid-min(xgrid))+0.1*dnorm(xgrid-2.5,sd=1)+dnorm(max(xgrid)-xgrid)))
                    ))
dev.off()
me <- sum(xgrid*f(xgrid))/sum(f(xgrid))
sqrt(sum((xgrid-me)^2*f(xgrid))/sum(f(xgrid)))




tplot(xgrid,dbeta((xgrid-min(xgrid))/diff(range(xgrid)), 0.2,0.2))





bounds <- c(-1, 1)
tra1 <- function(x){
    x[x<=bounds[1]] <- bounds[1]
    x[x>=bounds[2]] <- bounds[2]
    x
}
## slower:
tra2 <- function(x){
    lx <- length(x)
    ind <- rinterval(lx, x, bounds)
    bounds[1]*(ind==0)+bounds[2]*(ind==2)+x*(ind==1)
}


bounds <- ((-2):(2-1))+0.5
##
itra2 <- function(x){
    rinterval(length(x), x, bounds)
}
## slower:
itra0 <- function(x){
    rowSums(sapply(bounds, function(i){x>=i}))
}
## slower:
itra1 <- function(x){
    boundse <- c(-Inf, bounds, Inf)
    rowSums(sapply(2:length(boundse), function(i){
        (x>boundse[i-1] & x<=boundse[i])*i
    }))-2
}
