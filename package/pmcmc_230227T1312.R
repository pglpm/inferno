## Author: PGL  Porta Mana
## Created: 2022-09-08T17:03:24+0200
## Last-Updated: 2023-05-13T20:22:57+0200
#########################################
## Inference of exchangeable variates (nonparametric density regression)
## using effectively-infinite mixture of product kernels
## Monte Carlo sampling
#########################################

#### USER INPUTS AND CHOICES ####
baseversion <- '_testnewthresh2' # Base name of output directory
nsamples <- 64L # 256 gives 4096 samples with 16 parallel runs
ncores <- 1
datafile <- 'traindataset.csv'
predictorfile <- 'predictors.csv'
predictandfile <- NULL # 'predictors.csv'
varinfofile <- 'varinfo.rds'
ndata <- 20 # set this if you want to use fewer data
shuffledata <- FALSE # useful if subsetting data
posterior <- TRUE # if set to FALSE it samples and plots prior samples
savetempsamples <- FALSE # save temporary MCMC samples
plottempdistributions <- FALSE # plot temporary sampled distributions
showdata <- TRUE # 'histogram' 'scatter' FALSE TRUE
plotmeans <- TRUE # plot frequency averages
totsamples <- 'all' # 'all' number of samples if plotting frequency averages
showsamples <- 100 # number of samples to show. Shown separately for posterior=F
##
niter0 <- 1024L # 3L # iterations to try
testLength <- TRUE
nthreshold <- 2 # multiple of threshold for acceptable number of burn-in samples
casualinitvalues <- FALSE
showhyperparametertraces <- FALSE ##
showsamplertimes <- FALSE ##
family <- 'Palatino'

#### Hyperparameters
nclusters <- 64L
alpha0 <- 2^((-3):3)
##
hwidth <- 2 # number of powers of 2 to consider in either direction
##
rmean0 <- 0
rvar0 <- 2^2
rshapein0 <- 1
rshapeout0 <- 1
rvarscales <- (1 * 2^((-hwidth):hwidth))^2
##
dmean0 <- 0
dvar0 <- 2^2
dshapein0 <- 1
dshapeout0 <- 1
dvarscales <- (1 * 2^((-hwidth):hwidth))^2
##
## imean0 <- 0
## ivar0 <- 2^2
## ishapein0 <- 1
## ishapeout0 <- 1
## ivarscales <- (1 * 2^((-hwidth):hwidth))^2
## ifunction <- 'qt'
### alternative
imean0 <- 0
ivar0 <- (7/8)^2
ishapein0 <- 1
ishapeout0 <- 1
ivarscales <- ((1/4) * 2^((-hwidth):hwidth))^2
ifunction <- 'qnorm'
##
bshapein0 <- 1
bshapeout0 <- 1


#### Packages and setup ####
## load customized plot functions
if(!exists('tplot')){source('~/work/tplotfunctions.R')}
##
## Read MCMC seed from command line
arguments <- as.integer(commandArgs(trailingOnly=TRUE))[1]
mcmcseed <- arguments
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
if(!is.na(arguments)){## one of many chains
    if(is.na(ncores) || is.null(ncores)){
        ncores <- 1
    }else if(ncores > 1){
        warning('MORE THAN 1 CHAIN IN WHAT LOOKS LIKE A PARALLEL JOB')
    }
}else if(is.na(ncores) || is.null(ncores)){
    ncores <- 6
}
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

## if(Sys.info()['nodename']=='luca-HP-Z2-G9'){origdir <- '../'}else{origdir <- ''}
if(is.na(arguments)){origdir <- ''}else{origdir <- '../'}
source(paste0(origdir,'functionsmcmc_2212120902.R')) # load functions for post-MCMC calculations
## varinfo <- data.matrix(read.csv(paste0(origdir,varinfofile), row.names=1))
varinfo <- readRDS(paste0(origdir,varinfofile))

variate <- lapply(variatetypes, function(x)names(varinfo[['type']])[varinfo[['type']]==x])
len <- lapply(variate,length)
names(variate) <- names(len) <- variatetypes

nalpha <- length(alpha0)
nrvarscales <- length(rvarscales)
ndvarscales <- length(dvarscales)
novarscales <- length(ovarscales)
nivarscales <- length(ivarscales)


data0 <- fread(paste0(origdir,datafile), sep=',')
if(!all(unlist(variate) %in% colnames(data0))){cat('\nERROR: variates missing from datafile')}
data0 <- data0[, unlist(variate), with=F]
## shuffle data
if(exists('shuffledata') && shuffledata){data0 <- data0[sample(1:nrow(data0))]}
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(data0)}
data0 <- data0[1:ndata]

basename <- paste0(baseversion,'-V',sum(unlist(len)),'-D',ndata,'-K',nclusters,'-I',nsamples)
##
if(!is.na(arguments)){
    dirname <- ''
}else{
    dirname <- paste0(basename,'/')
    dir.create(dirname)
}

if(mcmcseed == 1){saveRDS(varinfo,file=paste0(dirname,'_varinfo-R',basename,'.rds'))}

source(paste0(origdir,'functionsmcmc_2212120902.R')) # load functions for post-MCMC ca

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


#################################
## Setup for Monte Carlo sampling
#################################

for(obj in c('constants', 'datapoints', 'inits', 'finitemix', 'finitemixnimble', 'Cfinitemixnimble', 'confnimble', 'mcsampler', 'Cmcsampler')){if(exists(obj)){do.call(rm,list(obj))}}
gc()

## data
datapoints = c(
    ## real
    if(len$R > 0){list( Rdata = transf(data0[,variate$R,with=F], varinfo) )},
    ## ## one-censored
    ## if(len$O > 0){list(
    ##                Odata = transf(data0[,variate$O,with=F], varinfo, Oout='data'),
    ##                Oleft = transf(data0[,variate$O,with=F], varinfo, Oout='left'),
    ##                Oaux = matrix(1L, nrow=ndata, ncol=len$O)
    ##               )},
    ## Censored
    if(len$D > 0){list(
                   Ddata = transf(data0[,variate$D,with=F], varinfo, Dout='data'),
                   Dleft = transf(data0[,variate$D,with=F], varinfo, Dout='left'),
                   Dright = transf(data0[,variate$D,with=F], varinfo, Dout='right'),
                   Daux = matrix(1L, nrow=ndata, ncol=len$D)
                  )},
    ## integer
    if(len$I > 0){list(
                   Ileft = transf(data0[,variate$I,with=F], varinfo, Iout='left', Ifunction=ifunction),
                   Iright = transf(data0[,variate$I,with=F], varinfo, Iout='right', Ifunction=ifunction),
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
    list(nalpha = nalpha),
    if(len$R > 0){ list(Rn = len$R,
                        nrvarscales = nrvarscales
                        ) },
    ## if(len$O > 0){ list(On = len$O,
    ##                     novarscales = novarscales
    ##                     ) },
    if(len$D > 0){ list(Dn = len$D,
                        ndvarscales = ndvarscales
                        ) },
    if(len$I > 0){ list(In = len$I,
                        nivarscales = nivarscales
                        ) },
    if(len$B > 0){ list(Bn = len$B) },
    if(len$C > 0){ list(Cn = len$C, Cmaxn = Cmaxn) },
    if(posterior){ list(ndata = ndata)}
)

##
## hyperparameters and some initial values
initsFunction <- function(){
    probalpha0 <- rep(1/nalpha, nalpha)
    walpha0 <- matrix(alpha0/nclusters, nrow=nalpha, ncol=nclusters)
    Alphaindex <- extraDistr::rcat(n=1, prob=probalpha0)
    W <- rdirch(n=1, alpha=walpha0[Alphaindex,])
    ##
    if(H$Rlen > 0){
        Rmean0 <- rep(H$Rmean0, H$Rlen)
        Rmeanv0 <- rep(H$Rmeanv0, H$Rlen)
        Rshapein0 <- rep(H$Rshapein0, H$Rlen)
        Rshapeout0 <- rep(H$Rshapeout0, H$Rlen)
        Rratev0 <- rep(H$Rratev0, H$Rlen)
        ##
        Rratev1 <- rinvgamma(n=H$Rlen, shape=Rshapein0, rate=Rratev0)

        Rvarscale <- rvarscales
        rprobvarscale0 <-rep(1/nrvarscales, nrvarscales)
        sRvarscaleindex <- extraDistr::rcat(n=H$Rlen, prob=rprobvarscale0)
        sRrate <- rinvgamma(n=H$Rlen, shape=Rshapein0, rate=Rvarscale[sRvarscaleindex])
        sRmean <- matrix(rnorm(n=nclusters*H$Rlen, mean=Rmean0, sd=sqrt(Rvar0)),
                         nrow=H$Rlen, ncol=nclusters)
        sRvar <- matrix(rinvgamma(n=nclusters*H$Rlen, shape=Rshapeout0, rate=sRrate),
                         nrow=H$Rlen, ncol=nclusters)
        }
    ##
    ## if(H$Olen > 0){
    ##     Omean0 <- rep(omean0, H$Olen)
    ##     Ovar0 <- rep(ovar0, H$Olen)
    ##     Oshapein0 <- rep(oshapein0, H$Olen)
    ##     Oshapeout0 <- rep(oshapeout0, H$Olen)
    ##     Ovarscale <- ovarscales
    ##     oprobvarscale0 <-rep(1/novarscales, novarscales)
    ##     sOvarscaleindex <- extraDistr::rcat(n=H$Olen, prob=oprobvarscale0)
    ##     sOrate <- rinvgamma(n=H$Olen, shape=Oshapein0, rate=Ovarscale[sOvarscaleindex])
    ##     sOmean <- matrix(rnorm(n=nclusters*H$Olen, mean=Omean0, sd=sqrt(Ovar0)),
    ##                      nrow=H$Olen, ncol=nclusters)
    ##     sOvar <- matrix(rinvgamma(n=nclusters*H$Olen, shape=Oshapeout0, rate=sOrate),
    ##                      nrow=H$Olen, ncol=nclusters)
    ##     }
    ##
    if(H$Dlen > 0){
        Dmean0 <- rep(dmean0, H$Dlen)
        Dvar0 <- rep(dvar0, H$Dlen)
        Dshapein0 <- rep(dshapein0, H$Dlen)
        Dshapeout0 <- rep(dshapeout0, H$Dlen)
        Dvarscale <- dvarscales
        dprobvarscale0 <-rep(1/ndvarscales, ndvarscales)
        sDvarscaleindex <- extraDistr::rcat(n=H$Dlen, prob=dprobvarscale0)
        sDrate <- rinvgamma(n=H$Dlen, shape=Dshapein0, rate=Dvarscale[sDvarscaleindex])
        sDmean <- matrix(rnorm(n=nclusters*H$Dlen, mean=Dmean0, sd=sqrt(Dvar0)),
                         nrow=H$Dlen, ncol=nclusters)
        sDvar <- matrix(rinvgamma(n=nclusters*H$Dlen, shape=Dshapeout0, rate=sDrate),
                         nrow=H$Dlen, ncol=nclusters)
        }
    ##
    if(H$Ilen > 0){
        Imean0 <- rep(imean0, H$Ilen)
        Ivar0 <- rep(ivar0, H$Ilen)
        Ishapein0 <- rep(ishapein0, H$Ilen)
        Ishapeout0 <- rep(ishapeout0, H$Ilen)
        Ivarscale <- ivarscales
        iprobvarscale0 <-rep(1/nivarscales, nivarscales)
        sIvarscaleindex <- extraDistr::rcat(n=H$Ilen, prob=iprobvarscale0)
        sIrate <- rinvgamma(n=H$Ilen, shape=Ishapein0, rate=Ivarscale[sIvarscaleindex])
        sImean <- matrix(rnorm(n=nclusters*H$Ilen, mean=Imean0, sd=sqrt(Ivar0)),
                         nrow=H$Ilen, ncol=nclusters)
        sIvar <- matrix(rinvgamma(n=nclusters*H$Ilen, shape=Ishapeout0, rate=sIrate),
                         nrow=H$Ilen, ncol=nclusters)
        }
    ##
    if(H$Blen > 0){
        Bshapeout0 <- rep(bshapeout0, H$Blen)
        Bshapein0 <- rep(bshapein0, H$Blen)
        sBprob <- matrix(rbeta(n=nclusters*H$Blen, shape1=Bshapein0, shape2=Bshapeout0),
                         nrow=H$Blen, ncol=nclusters)
    }
    ##
    ##
    c(
        list( # distribution over concentration parameter
            probalpha0 = probalpha0,
            walpha0 = walpha0,
            Alphaindex = Alphaindex,
            W = W
        ),
        if(H$Rlen > 0){list( # real variate
                     Rmean0 = Rmean0,
                     Rvar0 = Rvar0,
                     Rshapein0 = Rshapein0,
                     Rshapeout0 = Rshapeout0,
                     Rvarscale = Rvarscale,
                     rprobvarscale0 = rprobvarscale0,
                     Rvarscaleindex = sRvarscaleindex,
                     Rrate = sRrate,
                     Rmean = sRmean,
                     Rvar = sRvar
                 )},
        ## if(H$Olen > 0){c(list( # one-bounded variate
        ##              Omean0 = Omean0,
        ##              Ovar0 = Ovar0,
        ##              Oshapein0 = Oshapein0,
        ##              Oshapeout0 = Oshapeout0,
        ##              Ovarscale = Ovarscale,
        ##              oprobvarscale0 = oprobvarscale0,
        ##              Ovarscaleindex = sOvarscaleindex,
        ##              Orate = sOrate,
        ##              Omean = sOmean,
        ##              Ovar = sOvar
        ##              ),
        ##              if(posterior){list(
        ##              Odata = transf(data0[,variate$O,with=F], varinfo, Oout='init')
        ##                            )}
        ##              )},
        if(H$Dlen > 0){c(list( # censored variate
                     Dmean0 = Dmean0,
                     Dvar0 = Dvar0,
                     Dshapein0 = Dshapein0,
                     Dshapeout0 = Dshapeout0,
                     Dvarscale = Dvarscale,
                     dprobvarscale0 = dprobvarscale0,
                     Dvarscaleindex = sDvarscaleindex,
                     Drate = sDrate,
                     Dmean = sDmean,
                     Dvar = sDvar
                     ),
                     if(posterior){list(
                     Ddata = transf(data0[,variate$D,with=F], varinfo, Dout='init')
                                   )}
                     )},
        if(H$Ilen > 0){c(list( # discrete ordinal variate
                     Imean0 = Imean0,
                     Ivar0 = Ivar0,
                     Ishapein0 = Ishapein0,
                     Ishapeout0 = Ishapeout0,
                     Ivarscale = Ivarscale,
                     iprobvarscale0 = iprobvarscale0,
                     Ivarscaleindex = sIvarscaleindex,
                     Irate = sIrate,
                     Imean = sImean,
                     Ivar = sIvar
                     ),
                     if(posterior){list(
                     Icont = transf(data0[,variate$I,with=F], varinfo, Iout='init', Ifunction=ifunction)
                                   )}
                     )},
        if(H$Blen > 0){list( # binay variate
                     Bshapeout0 = Bshapeout0,
                     Bshapein0 = Bshapein0,
                     Bprob = sBprob
                 )},
        if(H$Clen > 0){list( # categorical variate ***must be fixed***
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
    Alphaindex ~ dcat(prob=probalpha0[1:nalpha])
    W[1:nclusters] ~ ddirch(alpha=walpha0[Alphaindex, 1:nclusters])
    ##
    if(H$Rlen > 0){# real variates
        for(v in 1:Rn){
            Rmean1[v] ~ dnorm(mean=Rmean0, var=Rmeanv0)
            Rratev1[v] ~ dinvgamma(shape=Rshapein0, rate=Rratev0)
            Rvarv1[v] ~ dinvgamma(shape=Rshapeout0, rate=Rratev1[v])
        }
    }
    if(H$Llen > 0){# censored variates
        for(v in 1:Ln){
            Lmean1[v] ~ dnorm(mean=Lmean0, var=Lmeanv0)
            Lratev1[v] ~ dinvgamma(shape=Lshapein0, rate=Lratev0)
            Lvarv1[v] ~ dinvgamma(shape=Lshapeout0, rate=Lratev1[v])
        }
    }
    ##
    for(k in 1:nclusters){
        if(H$Rlen > 0){# continuous variates
            for(v in 1:Rn){
                Rmean[v, k] ~ dnorm(mean=Rmean1[v], var=Rmeanv0)
                Rrate[v, k] ~ dinvgamma(shape=Rshapein0, rate=Rvarv1[v])
                Rvar[v, k] ~ dinvgamma(shape=Rshapeout0, rate=Rrate[v, k])
            }
        }
        if(H$Llen > 0){# latent-based variates
            for(v in 1:Ln){
                Lmean[v, k] ~ dnorm(mean=Lmean1[v], var=Lmeanv0)
                Lrate[v, k] ~ dinvgamma(shape=Lshapein0, rate=Lvarv1[v])
                Lvar[v, k] ~ dinvgamma(shape=Lshapeout0[v], rate=Lrate[v, k])
            }
        }
        if(H$Blen > 0){# binary variates
            for(v in 1:Bn){
                Bprob[v, k] ~ dbeta(shape1=Bshapein0[v], shape2=Bshapeout0[v])
            }
        }
        if(H$Clen > 0){# nominal variates
            for(v in 1:Cn){
                Cprob[v, k, 1:Cnmax] ~ ddirch(alpha=Calpha0[v, 1:Cnmax])
            }
        }
    }
    ##
    if(posterior){# cluster occupations
        for(d in 1:ndata){
            K[d] ~ dcat(prob=W[1:nclusters])
            ##
            if(H$Rlen > 0){# continuous variates
                for(v in 1:Rn){
                    Rdata[d, v] ~ dnorm(mean=Rmean[v, K[d]], var=Rvar[v, K[d]])
                }
            }
            if(H$Llen > 0){# latent-based variates
                for(v in 1:Ln){
                    Laux[d, v] ~ dconstraint(Ldata[d, v] >= Lleft[d, v] & Ldata[d, v] <= Lright[d, v])
                    Ldata[d, v] ~ dnorm(mean=Lmean[v, K[d]], var=Lvar[v, K[d]])
                }
            }
            if(H$Blen > 0){# binary variates
                for(v in 1:Bn){
                    Bdata[d, v] ~ dbern(prob=Bprob[v, K[d]])
                }
            }
            if(H$Clen > 0){# categorical variates
                for(v in 1:Cn){
                    Cdata[d, v] ~ dcat(prob=Cprob[v, K[d], 1:Cnmax])
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
                                   list(walpha0=c(nalpha,nclusters)),
                                   if(H$Rlen > 0){list(
                                                     Rmean=c(H$Rlen,nclusters),
                                                     Rvar=c(H$Rlen,nclusters),
                                                     Rrate=H$Rlen,
                                                     Rvarscaleindex=H$Rlen
                                                 )},
                                   ## if(H$Olen > 0){list(
                                   ##                   Omean=c(H$Olen,nclusters),
                                   ##                   Ovar=c(H$Olen,nclusters),
                                   ##                   Orate=H$Olen,
                                   ##                   Ovarscaleindex=H$Olen
                                   ##               )},
                                   if(H$Dlen > 0){list(
                                                     Dmean=c(H$Dlen,nclusters),
                                                     Dvar=c(H$Dlen,nclusters),
                                                     Drate=H$Dlen,
                                                     Dvarscaleindex=H$Dlen
                                                 )},
                                   if(H$Ilen > 0){list(
                                                     Imean=c(H$Ilen,nclusters),
                                                     Ivar=c(H$Ilen,nclusters),
                                                     Irate=H$Ilen,
                                                     Ivarscaleindex=H$Ilen
                                                 )},
                                   if(H$Blen > 0){list(
                                                     Bprob=c(H$Blen,nclusters)
                                                 )},
                                   if(H$Clen > 0){list(
                                                     Cprob=c(H$Clen,nclusters,Cmaxn)
                                                 )},
                                   if(posterior){c(
                                       list(K=ndata),
                                       if(H$Rlen > 0){list(Rdata=c(ndata,H$Rlen))},
                                       ## if(H$Olen > 0){list(
                                       ##                   Odata=c(ndata,H$Olen),
                                       ##                   Oleft=c(ndata,H$Olen),
                                       ##                   Oaux=c(ndata,H$Olen)
                                       ##               )},
                                       if(H$Dlen > 0){list(
                                                         Ddata=c(ndata,H$Dlen),
                                                         Dleft=c(ndata,H$Dlen),
                                                         Dright=c(ndata,H$Dlen),
                                                         Daux=c(ndata,H$Dlen)
                                                     )},
                                       if(H$Ilen > 0){list(
                                                         Icont=c(ndata,H$Ilen),
                                                         Ileft=c(ndata,H$Ilen),
                                                         Iright=c(ndata,H$Ilen),
                                                         Iaux=c(ndata,H$Ilen)
                                                     )},
                                       if(H$Blen > 0){list(Bdata=c(ndata,H$Blen))},
                                       if(H$Clen > 0){list(Cdata=c(ndata,H$Clen))}
                                   )}
                               )
                               )
##
Cfinitemixnimble <- compileNimble(finitemixnimble, showCompilerOutput=FALSE)
gc()

confnimble <- configureMCMC(Cfinitemixnimble, #nodes=NULL,
                            monitors=c('W',
                                       if(H$Rlen > 0){c('Rmean', 'Rvar')},
                                       ## if(H$Olen > 0){c('Omean', 'Ovar')},
                                       if(H$Dlen > 0){c('Dmean', 'Dvar')},
                                       if(H$Ilen > 0){c('Imean', 'Ivar')},
                                       if(H$Blen > 0){c('Bprob')},
                                       if(H$Clen > 0){c('Cprob')}
                                       ),
                            monitors2=c( 'Alphaindex',
                                        if(posterior){'K'},
                                        if(showhyperparametertraces){c(
                                            if(H$Rlen > 0){ c( 'Rvarscaleindex') },
                                            ## if(H$Olen > 0){ c( 'Ovarscaleindex') },
                                            if(H$Dlen > 0){ c( 'Dvarscaleindex') },
                                            if(H$Ilen > 0){ c( 'Ivarscaleindex') }
                                        )}
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
time0 <- Sys.time()
traces <- NULL
niterfinal <- niter0
nitertot <- 0L
stage <- 0L
continue <- TRUE

set.seed(mcmcseed+1000)
Cfinitemixnimble$setInits(initsFunction())
while(continue){
    showsamplertimes0 <- showsamplertimes && (stage==0)
    showhyperparametertraces0 <- showhyperparametertraces && (stage==0)
    if(testLength){
        nsamples0 <- 1L
        nburnin <- 0L
        reset <- FALSE
        niter <- niterfinal - nitertot
    }else{
        continue <- FALSE
        nsamples0 <- nsamples
        nburnin <- niterfinal-1L
        niter <- niterfinal
        reset <- TRUE
        traces <- NULL
        set.seed(mcmcseed+1000)
    }

    cat(paste0('Iterations: ', niter),'\n')
    mcsamples <- finalstate <- NULL
    calctime <- Sys.time()
    for(achain in 1:nsamples0){
        cat(paste0('chain: ', achain,', est. remaining time: ')); print((Sys.time()-calctime)/(achain-1)*(nsamples0-achain+1))
        if(!continue){
        ## set.seed(mcmcseed+achain+999)
            Cfinitemixnimble$setInits(initsFunction())
        }
        todelete <- Cmcsampler$run(niter=niter, thin=1, thin2=1, nburnin=nburnin, time=showsamplertimes0, reset=reset, resetMV=TRUE)
        mcsamples <- rbind(mcsamples, as.matrix(Cmcsampler$mvSamples))
        finalstate <- rbind(finalstate, as.matrix(Cmcsampler$mvSamples2))
    }
    cat('\nTime MCMC: ')
    print(Sys.time() - calctime)

    if(any(!is.finite(mcsamples))){cat('\nWARNING: SOME NON-FINITE OUTPUTS')}
    
    if(showsamplertimes){
        samplertimes <- Cmcsampler$getTimes()
        names(samplertimes) <- sapply(confnimble$getSamplers(),function(x)x$target)
        sprefixes <- unique(sub('^([^[]+)(\\[.*\\])', '\\1', names(samplertimes)))
        cat(paste0('\nSampler times:\n'))
        print(sort(sapply(sprefixes, function(x)sum(samplertimes[grepl(x,names(samplertimes))])),decreasing=T))
    }
    ##
    if(showhyperparametertraces){
        occupations <- apply(finalstate[,grepl('^K\\[', colnames(finalstate)),drop=F], 1, function(xxx){length(unique(xxx))})
        cat(paste0('\nSTATS OCCUPIED CLUSTERS:\n'))
        print(summary(occupations))
        ##
        pdff(paste0(dirname,'traces_hyperparameters-',mcmcseed,'-',stage))
        tplot(y=occupations, ylab='occupied clusters',xlab=NA,ylim=c(0,nclusters))
        histo <- thist(occupations,n='i')
        tplot(x=histo$mids,y=histo$density,xlab='occupied clusters',ylab=NA)
        tplot(y=log2(alpha0[finalstate[,'Alphaindex']]), ylab='log2(alpha)',xlab=NA)
        histo <- thist(log2(alpha0[finalstate[,'Alphaindex']]))
        tplot(x=histo$mids,y=histo$density,xlab='log2(alpha)',ylab='')
        for(vtype in c('R','I',#'O',
                       'D')){
            if(len[[vtype]] > 0){
                for(togrep in c('varscaleindex')){
                    for(v in colnames(finalstate)[grepl(paste0('^',vtype,togrep,'\\['), colnames(finalstate))]){
                        tplot(y=finalstate[,v],ylab=v,xlab=NA,ylim=c(1,(2*hwidth+1)))
                    }
                }
            }
        }
        dev.off()
    }
    ##
    finalstate <- c(mcsamples[nrow(mcsamples),], finalstate[nrow(finalstate),])
    ## Check how many "clusters" were occupied. Warns if too many
    occupations <- finalstate[grepl('^K\\[', names(finalstate))]
    usedclusters <- length(unique(occupations))
    if(usedclusters > nclusters-5){cat('\nWARNING: TOO MANY CLUSTERS OCCUPIED')}
    cat(paste0('\nOCCUPIED CLUSTERS: ', usedclusters, ' OF ', nclusters),'\n')
    ##
    ## Diagnostics
    ## Log-likelihood
    diagntime <- Sys.time()
    ll <- colSums(log(samplesFDistribution(Y=data.matrix(data0), X=NULL, mcsamples=mcsamples, varinfo=varinfo, jacobian=FALSE)), na.rm=T) #- sum(log(invjacobian(data.matrix(data0), varinfo)), na.rm=T)
    lld <- colSums(log(samplesFDistribution(Y=data.matrix(data0[,..predictands]), X=data.matrix(data0[,..predictors]), mcsamples=mcsamples, varinfo=varinfo, jacobian=FALSE)), na.rm=T) # - sum(log(invjacobian(data.matrix(data0[,..predictands]), varinfo)), na.rm=T)
    lli <- colSums(log(samplesFDistribution(Y=data.matrix(data0[,..predictors]), X=data.matrix(data0[,..predictands]), mcsamples=mcsamples, varinfo=varinfo, jacobian=FALSE)), na.rm=T) #- sum(log(invjacobian(data.matrix(data0[,..predictors]), varinfo)), na.rm=T)
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
        cat(paste0('\nESSs: ',paste0(round(diagnESS), collapse=', ')))
        cat(paste0('\nIATs: ',paste0(round(diagnIAT), collapse=', ')))
        cat(paste0('\nBMKs: ',paste0(round(diagnBMK,3), collapse=', ')))
        cat(paste0('\nMCSEs: ',paste0(round(diagnMCSE,2), collapse=', ')))
        cat(paste0('\nStationary: ',paste0(diagnStat, collapse=', ')))
        cat(paste0('\nBurn-in I: ',paste0(diagnBurn, collapse=', ')))
        cat(paste0('\nBurn-in II: ',diagnBurn2))
        cat(paste0('\nProposed thinning: ',paste0(diagnThin, collapse=', ')),'\n')

        cat('\nTime diagnostics: ')
        print(Sys.time() - diagntime)
        ##
        #########################################
        #### CHECK IF WE NEED TO SAMPLE MORE ####
        #########################################
    lengthmeasure <- ceiling(nthreshold * max(diagnBurn2,diagnIAT))
    ## lengthmeasure <- ceiling(max(diagnBurn2))
    if(testLength){
            if(niterfinal < lengthmeasure){
                cat(paste0('Number of iterations/threshold is too small: ', signif(niterfinal/lengthmeasure,2)), '. ')
                niterfinal <- lengthmeasure
                cat(paste0('Increasing to ', niterfinal), '\n')
                nitertot <- nitertot + niter
                stage <- stage + 1L
            }else{
                cat(paste0('Number of iterations/threshold: ', signif(niterfinal/lengthmeasure,2)), '\n')
                testLength <- FALSE
            }
        }else{
            saveRDS(mcsamples, file=paste0(dirname,'_mcsamples-R',basename,'--',mcmcseed,'-F.rds'))
            saveRDS(traces, file=paste0(dirname,'_mctraces-R',basename,'--',mcmcseed,'-F.rds'))
            stage <- 'F'
        }
        #########################################
        #### END CHECK                       ####
        #########################################

        ##
        ###############
        #### PLOTS ####
        ###############
        tracegroups <- list(loglikelihood=1,
                            'predictand given predictor'=2,
                            'predictor given predictand'=3
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
    ##
    ## Samples of marginal frequency distributions
    if(plottempdistributions | !continue){
        subsamples <- (if(totsamples=='all'){1:nrow(mcsamples)}else{round(seq(1, nrow(mcsamples), length.out=totsamples))})
        showsubsample <- round(seq(1, length(subsamples), length.out=showsamples))
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
            plotsamples <- samplesFDistribution(Y=Xgrid, mcsamples=mcsamples, varinfo=varinfo, subsamples=subsamples, jacobian=TRUE)
            ##
            if(posterior){
                par(mfrow=c(1,1))
                ymax <- tquant(apply(plotsamples[,showsubsample],2,function(x){tquant(x,31/32)}),31/32, na.rm=T)
                if((showdata=='histogram' || showdata==TRUE) && !contvar){
                    datum <- data0[[v]]
                    datum <- datum[!is.na(datum)]
                    nh <- (varinfo[['max']][v]-varinfo[['min']][v])/(varinfo[['n']][v]-1)
                    nh <- seq(varinfo[['min']][v]-nh/2, varinfo[['max']][v]+nh/2, length.out=varinfo[['n']][v]+1)
                    histo <- thist(datum, n=nh)
                    ymax <- max(ymax,histo$counts/sum(histo$counts))
                }
                if(!(varinfo[['type']][v] %in% c('O','D'))){
                    ##
                    tplot(x=Xgrid, y=plotsamples[,showsubsample], type='l', col=5, alpha=7/8, lty=1, lwd=2,
                          xlab=paste0(v, (if(varinfo[['type']][v] %in% c('I','B','C')){' (discrete)'}else{' (continuous)'})),
                          ylab=paste0('frequency', (if(varinfo[['type']][v] %in% c('R','O','D')){' density'}else{''})),
                          ylim=c(0, ymax), family=family)
                    ##
                    if(plotmeans){
                        tplot(x=Xgrid, y=rowMeans(plotsamples, na.rm=T), type='l', col=1, alpha=0.25, lty=1, lwd=4, add=T)
                    }
                }else{ # plot of a continuous doubly-bounded variate
                    interior <- which(Xgrid > varinfo[['tmin']][v] & Xgrid < varinfo[['tmax']][v])
                    tplot(x=Xgrid[interior], y=plotsamples[interior,showsubsample], type='l', col=5, alpha=7/8, lty=1, lwd=2,
                          xlab=paste0(v, ' (continuous with deltas)'),
                          ylab=paste0('frequency (density)'),
                          ylim=c(0, ymax), family=family)
                    if(length(interior) < length(Xgrid)){
                        tplot(x=Xgrid[-interior], y=plotsamples[-interior,showsubsample,drop=F]*ymax, type='p', pch=2, cex=2, col=5, alpha=7/8, lty=1, lwd=2, xlab=paste0(v), ylab=paste0('frequency'), ylim=c(0, ymax), family=family,add=T)
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
                    ##
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
                            tplot(x=histo$mids, y=histo$density*histomax, col=yellow, alpha=2/4, border=darkgrey, border.alpha=3/4, lty=1, lwd=4, family=family, ylim=c(0,NA), add=TRUE)
                        }else{
                            tplot(x=histo$mids, y=histo$counts/sum(histo$counts), col=yellow, alpha=2/4, border=darkgrey, border.alpha=3/4, lty=1, lwd=4, family=family, ylim=c(0,NA), add=TRUE)
                        }
                    }else{ # histogram for censored variate
                        interior <- which(datum > varinfo[['tmin']][v] & datum < varinfo[['tmax']][v])
                        histo <- thist(datum[interior], n=max(10,round(length(interior)/64)))
                        interiorgrid <- which(Xgrid > varinfo[['tmin']][v] & Xgrid < varinfo[['tmax']][v])
                        histomax <- 1#max(rowMeans(plotsamples)[interiorgrid])/max(histo$density)
                        tplot(x=histo$mids, y=histo$density*histomax, col=yellow, alpha=2/4, border=darkgrey, border.alpha=3/4, lty=1, lwd=4, family=family, ylim=c(0,NA), add=TRUE)
                        ##
                        pborder <- sum(datum <= varinfo[['tmin']][v])/length(datum)
                        if(pborder > 0){
                            tplot(x=varinfo[['tmin']][v], y=pborder*ymax, type='p', pch=0, cex=2, col=yellow, alpha=0, lty=1, lwd=5, family=family, ylim=c(0,NA), add=TRUE)
                        }
                        ##
                        pborder <- sum(datum >= varinfo[['tmax']][v])/length(datum)
                        if(pborder > 0){
                            tplot(x=varinfo[['tmax']][v], y=pborder*ymax, type='p', pch=0, cex=2, col=yellow, alpha=0, lty=1, lwd=5, family=family, ylim=c(0,NA), add=TRUE)
                        }
                    }
                }else if((showdata=='scatter')|(showdata==TRUE & contvar)){
                    datum <- data0[[v]]
                    datum <- datum[!is.na(datum)]
                    diffdatum <- c(apply(cbind(c(0,diff(datum)),c(diff(datum),0)),1,min))/2
                    scatteraxis(side=1, n=NA, alpha='88',
                                ext=5, x=datum+runif(length(datum),
                                                     min=-min(diff(sort(unique(datum))))/4,
                                                     max=min(diff(sort(unique(datum))))/4),
                                col=yellow)
                }
                ## fiven <- sapply(c('datamin','Q1','Q2','Q3','datamax'),function(xxx){varinfo[[xxx]][v]})
                fiven <- fivenum(datum)
                abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'),lwd=4)
            }else{
                fiven <- sapply(c('datamin','Q1','Q2','Q3','datamax'),function(xxx){varinfo[[xxx]][v]})
                par(mfrow=c(floor(sqrt(showsamples)),floor(sqrt(showsamples))),mar = c(0,0,0,0))
                ##
                for(aplot in showsubsample[-1]){
                    tplot(x=Xgrid, y=plotsamples[,aplot], type='l', col=1, lty=1, lwd=c(1,1), xlab=NA, ylab=NA, ylim=c(0, NA), family=family,
                          xticks=NA, yticks=NA,
                          mar=c(1,1,1,1))
                    abline(h=0, col=7, lwd=1)
                    abline(v=fiven, col=paste0(palette()[c(2,4,5,4,2)], '44'))
                    ## if(aplot==1){ text(Xgrid[1], par('usr')[4]*0.9, variateinfo[variate==avar,type], pos=4)}
                    ## if(aplot==2){ text(Xgrid[1], par('usr')[4]*0.9, paste(signif(c(rg,diff(rg)),2),collapse=' -- '), pos=4)}
                }
                ##
                tplot(x=Xgrid, y=rowMeans(plotsamples), type='l', col=3, lty=1, lwd=2, xlab=NA, ylab=NA, ylim=c(0, NA), family=family,
                      xticks=NA, yticks=NA,
                      mar=c(1,1,1,1))
                abline(h=0, col=7, lwd=1)
                abline(v=fiven, col=paste0(palette()[c(2,4,5,4,2)], '44'),lwd=4)
                text(sum(range(Xgrid))/2, par('usr')[4]*0.9, v)
            }
            }
    dev.off()
    }
    ##
    cat('\nTime MCMC+diagnostics: ')
    print(Sys.time() - calctime)
    ##

}

##
cat('\nTotal time: ')
print(Sys.time() - time0)

############################################################
## End MCMC
############################################################
plan(sequential)

stop('NONE. End of script')
