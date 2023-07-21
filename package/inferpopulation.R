inferpopulation <- function(data, auxmetadata, outputdir, nsamples=1200, nchains=120, nsamplesperchain, parallel=TRUE, niterini=1024, plottraces=TRUE, showclusterstraces=TRUE, seed=701, subsampledata, miniter=0, useOquantiles=TRUE, output=FALSE, cleanup=TRUE){

    cat('\n')
#### Determine the status of parallel processing
    if(!missing(parallel) && is.logical(parallel) && parallel){
        if(getDoParRegistered()){
            cat('Using already registered', getDoParName(), 'with', getDoParWorkers(), 'workers\n')
            ncores <- getDoParWorkers()
        }else{
            cat('No parallel backend registered.\n')
            ncores <- 1
        }
    }else if(!missing(parallel) && is.numeric(parallel) && parallel >= 2){
        if(getDoParRegistered()){
            ncores <- min(getDoParWorkers(), parallel)
            cat('Using already registered', getDoParName(), 'with', getDoParWorkers(), 'workers\n')
        }else{
            ## registerDoSEQ()
            ## cl <- makePSOCKcluster(ncores)
            cl <- makeCluster(parallel)
            registerDoParallel(cl)
            cat('Registered', getDoParName(), 'with', getDoParWorkers(), 'workers\n')
            ncores <- parallel
        }
    }else{
        cat('No parallel backend registered.\n')
        ncores <- 1
    }

    
#### Consistency checks for numbers of samples, chains, cores
    if(!missing(nsamples) && !missing(nchains) && missing(nsamplesperchain)){
        ## nsamples and nchains given
        if(nchains < ncores){ ncores <- nchains }
        nchainspercore <- ceiling(nchains/ncores)
        if(nchainspercore*ncores > nchains){
            nchains <- nchainspercore*ncores
            cat('Increasing number of chains to',nchains,'\n')
        }
        nsamplesperchain <- ceiling(nsamples/nchains)
        if(nsamplesperchain*nchains > nsamples){
            nsamples <- nchains*nsamplesperchain
            cat('Increasing number of samples to',nsamples,'\n')
        }
        ##
    }else if(!missing(nsamples) && missing(nchains) && !missing(nsamplesperchain)){
        ## nsamples and nsamplesperchain given
        nchains <- ceiling(nsamples/nsamplesperchain)
        if(nchains < ncores){ ncores <- nchains }
        nchainspercore <- ceiling(nchains/ncores)
        if(nchainspercore*ncores > nchains){
            nchains <- nchainspercore*ncores
        }
        if(nsamplesperchain*nchains > nsamples){
            nsamples <- nchains*nsamplesperchain
            cat('Increasing number of samples to',nsamples,'\n')
        }
        ##
    }else if(missing(nsamples) && !missing(nchains) && !missing(nsamplesperchain)){
        ## nchains and nsamplesperchain given
        if(nchains < ncores){ ncores <- nchains }
        nchainspercore <- ceiling(nchains/ncores)
        if(nchainspercore*ncores > nchains){
            nchains <- nchainspercore*ncores
            cat('Increasing number of chains to',nchains,'\n')
        }
        nsamples <- nchains * nsamplesperchain
        ##
    }else{
        stop('please specify exactly two among "nsamples", "nchains", "nsamplesperchain"')
    }

    ## ## set number of cores for parallel computation
    ## if(missing(ncores) || !is.numeric(ncores)){
    ##     warning('The number of cores has not been given.\nIt is much preferable that it be set by the user.')
    ##     ncores <- round(parallel::detectCores()/2)
    ##     cat('\nTrying to use ',ncores,' cores\n')
    ## }

    if(ncores < 2){ `%dochains%` <- `%do%` }else{ `%dochains%` <- `%dorng%` }
    
    timestart0 <- Sys.time()

    ## cat('\n')
    ## library('data.table')
    ## library('LaplacesDemon', include.only=NULL)
    ## library('foreach')
    ## library('doParallel')
    ## library('doRNG')


    ## auxmetadata
    if(is.character(auxmetadata) && file.exists(auxmetadata)){
        auxmetadata <- paste0(sub('.rds$', '', auxmetadata), '.rds')
        auxmetadata <- readRDS(auxmetadata)
    }

    ## read dataset
    datafile <- NULL
    if(missing(data) || is.null(data) || (is.logical(data) && data==FALSE)){
        message('Missing data: calculating prior distribution')
        data <- as.data.table(matrix(NA,nrow=1,ncol=nrow(auxmetadata),dimnames=list(NULL,auxmetadata[['name']])))
    }
    if(is.character(data) && file.exists(data)){
        datafile <- paste0(sub('.csv$', '', data), '.csv')
        data <- fread(datafile, na.strings='')
    }
    data <- as.data.table(data)
    if(!missing(subsampledata) && is.numeric(subsampledata)){
        ##@@TODO: use faster and memory-saving subsetting
        data <- data[sample(1:nrow(data), min(subsampledata,nrow(data)), replace=F),]
    }

    if(missing(outputdir) || outputdir==TRUE){
        outputdir <- paste0('_output_', datafile)
        outputdir <- paste0(sub('.csv$', '', outputdir))
    }

    ## Check consistency of variate names
    if(!all(auxmetadata[['name']] %in% colnames(data))){
        stop('Missing variates in data file.')
    }
    if(!all(colnames(data) %in% auxmetadata[['name']])){
        message('Data file has additional variates. Dropping them.')
        subvar <- intersect(colnames(data), auxmetadata[['name']])
        ## cat(subvar,'\n')
        data <- data[,..subvar]
        ## auxmetadata <- auxmetadata[name %in% subvar]
    }

    ## ## file to save output
    ## if(is.character(file) || (is.logical(file) && file)){ # must save to file
    ##     if(is.character(file)){
    ##         file <- paste0(sub('.rds$', '', file), '.rds')
    ##     }else{
    ##         if(file.exists('longrunsamples.rds')){
    ##             file <- paste0('longrunsamples_',format(Sys.time(), '%y%m%dT%H%M%S'),'.rds')
    ##         }else{
    ##             file <- 'longrunsamples.rds'
    ##         }
    ##     }
    ## }

##################################################
#### various internal parameters
    ## source('hyperparameters.R')
    ## niterini # initial iterations to try
#### Hyperparameters
    nclusters <- 64L # ****
    minalpha <- -3L
    maxalpha <- 3L
    Rshapelo <- 0.5
    Rshapehi <- 0.5
    Cshapelo <- 0.5
    Cshapehi <- 0.5
    Dshapelo <- 0.5
    Dshapehi <- 0.5
    Oshapelo <- 0.5
    Oshapehi <- 0.5
    Bshapelo <- 1
    Bshapehi <- 1
    ##
    nalpha <- length(minalpha:maxalpha)
    npoints <- nrow(data)

#### other options
    Alphatoslice <- TRUE
    Ktoslice <- FALSE
    RWtoslice <- FALSE
    ##
    ## showdata <- TRUE # 'histogram' 'scatter' FALSE TRUE
    plotmeans <- TRUE # plot frequency averages
    totsamples <- 'all' # 'all' number of samples if plotting frequency averages
    showsamples <- 100 # number of samples to show.
    showquantiles <- c(1,31)/32 # quantiles to show
    nclustersamples <- 128 ## number of samples of Alpha and occupations
    showsamplertimes <- FALSE ##
    family <- 'Palatino'
##################################################


    nameroot <- paste0(outputdir,'-V',nrow(auxmetadata),'-D',(if(npoints==1 && all(is.na(data))){0}else{npoints}),'-K',nclusters,'-S',nsamples)
    ##
    dirname <- paste0(nameroot,'/')
    dir.create(dirname)
    cat('\n',paste0(rep('*',max(nchar(dirname),26)),collapse=''),
        '\n Saving output in directory\n',dirname,'\n',
        paste0(rep('*',max(nchar(dirname),26)),collapse=''),'\n')
    nameroot <- basename(nameroot)

    saveRDS(auxmetadata,file=paste0(dirname,'_auxmetadata_copy.rds'))

    source('tplotfunctions.R')
    source('vtransform.R')
    source('samplesFDistribution.R')
    source('proposeburnin.R')
    source('proposethinning.R')
    source('plotFsamples.R')
    printtime <- function(tim){paste0(signif(tim,2),' ',attr(tim,'units'))}
    printnull <- function(message, outcon){
        sink(NULL,type='message')
        message(message, appendLF=FALSE)
        flush.console()
        sink(outcon,type='message')
    }

    ## printtime <- function(tim){sub('^Time difference of (.*)', '\\1', capture.output(print(tim)))}
    vn <- list()
    vnames <- list()
    for(atype in c('R','C','D','O','N','B')){
        vn[[atype]] <- length(auxmetadata[mcmctype == atype, name])
        vnames[[atype]] <- auxmetadata[mcmctype == atype, name]
    }
    ## ## choose predictands if not chosen
    ## if(is.null(predictands) || !exists('predictands')){
    ##     predictands <- sapply(vnames,function(xx)xx[1])
    ##     cat('\nSelf-choosing predictand variates:\n', predictands, '\n')
    ## }

    if(vn$N > 0){
        Nmaxn <- max(auxmetadata[mcmctype == 'N', Nvalues])
        Nalpha0 <- matrix(1e-100, nrow=vn$N, ncol=Nmaxn)
        for(avar in 1:length(vnames$N)){
            Nalpha0[avar, 1:auxmetadata[name == vnames$N[avar], Nvalues]] <- 1
        }
    }

#### CONSTANTS OF MONTE-CARLO SAMPLER
    constants <- c( list(
        nclusters = nclusters,
        npoints = npoints,
        nalpha = nalpha,
        probalpha0 = rep(1/nalpha, nalpha),
        basealphas = rep((2^(minalpha-1L))/nclusters, nclusters)
    ),
    if(vn$R > 0){# continuous
        list(Rn = vn$R,
             Rmean1 = rep(0, 1),
             Rvarm1 = rep(1, 1),
             Rvar1 = rep(1, 1),
             Rshapelo = rep(Rshapelo, 1),
             Rshapehi = rep(Rshapehi, 1)
             ) },
    if(vn$C > 0){# censored
        list(Cn = vn$C,
             Cmean1 = rep(0, 1),
             Cvarm1 = rep(1, 1),
             Cvar1 = rep(1, 1),
             Cshapelo = rep(Cshapelo, 1),
             Cshapehi = rep(Cshapehi, 1),
             Cleft = vtransform(data[,vnames$C, with=F], auxmetadata, Cout='left'),
             Cright = vtransform(data[,vnames$C, with=F], auxmetadata, Cout='right')
             ) },
    if(vn$D > 0){# discretized
        list(Dn = vn$D,
             Dmean1 = rep(0, 1),
             Dvarm1 = rep(1, 1),
             Dvar1 = rep(1, 1),
             Dshapelo = rep(Dshapelo, 1),
             Dshapehi = rep(Dshapehi, 1),
             Dleft = vtransform(data[,vnames$D, with=F], auxmetadata, Dout='left'),
             Dright = vtransform(data[,vnames$D, with=F], auxmetadata, Dout='right')
             ) },
    if(vn$O > 0){# ordinal
        list(On = vn$O,
             Omean1 = rep(0, 1),
             Ovarm1 = rep(1, 1),
             Ovar1 = rep(1, 1),
             Oshapelo = rep(Oshapelo, 1),
             Oshapehi = rep(Oshapehi, 1),
             Oleft = vtransform(data[,vnames$O, with=F], auxmetadata, Oout='left', useOquantiles=useOquantiles),
             Oright = vtransform(data[,vnames$O, with=F], auxmetadata, Oout='right', useOquantiles=useOquantiles)
             ) },
    if(vn$N > 0){# nominal
        list(Nn = vn$N,
             Nmaxn = Nmaxn,
             Nalpha0 = Nalpha0
             ) },
    if(vn$B > 0){# binary
        list(Bn = vn$B,
             Bshapelo = rep(Bshapelo, 1),
             Bshapehi = rep(Bshapehi, 1)
             ) }
    )

#### DATAPOINTS
    datapoints <- c(
        if(vn$R > 0){# continuous
            list(
                Rdata = vtransform(data[,vnames$R, with=F], auxmetadata)
            ) },
        if(vn$C > 0){# censored
            list(
                Caux = vtransform(data[,vnames$C, with=F], auxmetadata, Cout='aux'),
                Clat = vtransform(data[,vnames$C, with=F], auxmetadata, Cout='lat')
            ) },
        if(vn$D > 0){# discretized
            list(
                Daux = vtransform(data[,vnames$D, with=F], auxmetadata, Dout='aux')
            ) },
        if(vn$O > 0){# ordinal
            list(
                Oaux = vtransform(data[,vnames$O, with=F], auxmetadata, Oout='aux', useOquantiles=useOquantiles)
            ) },
        if(vn$N > 0){# nominal
            list(
                Ndata = vtransform(data[,vnames$N,with=F], auxmetadata, Nout='numeric')                
            ) },
        if(vn$B > 0){# binary
            list(
                Bdata = vtransform(data[,vnames$B,with=F], auxmetadata, Bout='numeric')                
            ) }
    )

#### Output info
    cat('Starting Monte Carlo sampling of',nsamples,'samples by',nchains,'chains across', ncores, 'cores.\n')
    cat(nsamplesperchain,'samples per chain,',nchainspercore,'chains per core.\n')
    cat('Core logs are being saved in individual files.\n')
    cat('Setting up samplers (this can take tens of minutes with many data or variates).\n...\r')
    ##cat('Estimating remaining time...\r')
    ## stopCluster(cluster)
    ## stopImplicitCluster()
    ## registerDoSEQ()
    ## ## cl <- makePSOCKcluster(ncores)
    ## if(ncores > 1){
    ##     cl <- makeCluster(ncores)
    ##     registerDoParallel(cl)
    ## }else{
    ##     registerDoSEQ()
    ## }
    ## toexport <- c('constants', 'datapoints', 'vn', 'vnames', 'nalpha', 'nclusters')
    ## toexport <- c('vtransform','samplesFDistribution','proposeburnin','proposethinning','plotFsamples')

#### BEGINNING OF FOREACH LOOP OVER CORES
    chaininfo <- foreach(acore=1:ncores, .combine=rbind, .inorder=FALSE)%dochains%{

        outcon <- file(paste0(dirname,'_log-',nameroot,'-',acore,'.log'), open = "w")
        sink(outcon)
        sink(outcon, type = "message")
        suppressPackageStartupMessages(library('data.table'))
        suppressPackageStartupMessages(library('nimble'))

        
#### Load and define various functions
        source('tplotfunctions.R')
        source('vtransform.R')
        source('samplesFDistribution.R')
        source('proposeburnin.R')
        source('proposethinning.R')
        source('plotFsamples.R')
        source('mcsubset.R')

        ## Function for diagnostics
        funMCSE <- function(x){
            if(length(x) >= 1000){
                LaplacesDemon::MCSE(x, method='batch.means')$se
            }else{
                LaplacesDemon::MCSE(x)
            }
        }

        ## ## To join MC samples in list form
        ## joinmc <- function(mc1, mc2){
        ##     if(is.null(mc1)){
        ##         mc2
        ##     }else{
        ##         mapply(function(xx,yy){
        ##             temp <- c(xx,yy)
        ##             dx <- dim(xx)[-length(dim(xx))]
        ##             dim(temp) <- c(dx, length(temp)/dx)
        ##             temp
        ##         },
        ##         mc1, mc2)
        ##     }
        ## }
        ## ## To remove iterations with non-finite values
        ## cleanmc <- function(mcx,toremove){
        ##     lapply(mcx,function(xx){
        ##         do.call('[',c(list(xx),rep(TRUE,length(dim(xx))-1), list(-toremove)) )
        ##     })
        ## }



        ## Parameter and function to test MCMC convergence
        if(miniter == 0){
            multcorr <- 2L
            thresholdfn <- function(diagnESS, diagnIAT, diagnBMK, diagnMCSE, diagnStat, diagnBurn, diagnBurn2, diagnThin){
                ceiling(2* max(diagnBurn2) + (nsamplesperchain-1L) * multcorr * ceiling(max(diagnIAT, diagnThin)))
            }
        }else if(miniter > 0){
            multcorr <- 2L
            thresholdfn <- function(diagnESS, diagnIAT, diagnBMK, diagnMCSE, diagnStat, diagnBurn, diagnBurn2, diagnThin){
                max( ceiling(2* max(diagnBurn2) + (nsamplesperchain-1L) * multcorr * ceiling(max(diagnIAT, diagnThin))),
                    miniter )
            }
        }else if(miniter < 0){
            multcorr <- 0L
            thresholdfn <- function(diagnESS, diagnIAT, diagnBMK, diagnMCSE, diagnStat, diagnBurn, diagnBurn2, diagnThin){
                -miniter
            }
            niterini <- min(-miniter, niterini)
        }else{
            stop('Invalid "miniter" argument.')
        }

        ## printtime <- function(tim){sub('^Time difference of (.*)', '\\1', capture.output(print(tim)))}
        printtime <- function(tim){paste0(signif(tim,2),' ',attr(tim,'units'))}
        printnull <- function(message, outcon){
            sink(NULL,type='message')
            message(message, appendLF=FALSE)
            flush.console()
            sink(outcon,type='message')
        }

        ## predictors <- setdiff(unlist(vnames), predictands)

#### CLUSTER REPRESENTATION OF FREQUENCY SPACE
        ## hierarchical probability structure
        finitemix <- nimbleCode({
            ## Component weights
            Alpha ~ dcat(prob=probalpha0[1:nalpha])
            alphas[1:nclusters] <- basealphas[1:nclusters] * 2^Alpha
            W[1:nclusters] ~ ddirch(alpha=alphas[1:nclusters])
            ##
            ## if(FALSE){
            ##     for(v in 1:nvar){
            ##         ## Rmean1[v] ~ dnorm(mean=0, var=1)
            ##         ## Rratem1[v] ~ dinvgamma(shape=shapehi0, rate=1)
            ##         Rvarm1[v] ~ dinvgamma(shape=shapelo0, rate=Rratem1[v])
            ##         ##
            ##         Shapelo[v] ~ dinvgamma(shape=minshape, scale=maxshape)
            ##         Shapehi[v] ~ dinvgamma(shape=minshape, scale=maxshape)
            ##         ## Rrate1[v] ~ dinvgamma(shape=shapehi0, rate=1)
            ##         ## Rvar1[v] ~ dinvgamma(shape=shapelo0, rate=Rrate1[v])
            ##     }
            ## }
            ## Probability density for the parameters of the components
            for(k in 1:nclusters){
                if(vn$R > 0){# continuous
                    for(v in 1:Rn){
                        Rmean[v, k] ~ dnorm(mean=Rmean1, var=Rvarm1)
                        Rrate[v, k] ~ dinvgamma(shape=Rshapehi, rate=Rvar1)
                        Rvar[v, k] ~ dinvgamma(shape=Rshapelo, rate=Rrate[v, k])
                    }
                }
                if(vn$C > 0){# censored
                    for(v in 1:Cn){
                        Cmean[v, k] ~ dnorm(mean=Cmean1, var=Cvarm1)
                        Crate[v, k] ~ dinvgamma(shape=Cshapehi, rate=Cvar1)
                        Cvar[v, k] ~ dinvgamma(shape=Cshapelo, rate=Crate[v, k])
                    }
                }
                if(vn$D > 0){# discretized
                    for(v in 1:Dn){
                        Dmean[v, k] ~ dnorm(mean=Dmean1, var=Dvarm1)
                        Drate[v, k] ~ dinvgamma(shape=Dshapehi, rate=Dvar1)
                        Dvar[v, k] ~ dinvgamma(shape=Dshapelo, rate=Drate[v, k])
                    }
                }
                if(vn$O > 0){# ordinal
                    for(v in 1:On){
                        Omean[v, k] ~ dnorm(mean=Omean1, var=Ovarm1)
                        Orate[v, k] ~ dinvgamma(shape=Oshapehi, rate=Ovar1)
                        Ovar[v, k] ~ dinvgamma(shape=Oshapelo, rate=Orate[v, k])
                    }
                }
                if(vn$N > 0){# nominal
                    for(v in 1:Nn){
                        Nprob[v, k, 1:Nmaxn] ~ ddirch(alpha=Nalpha0[v, 1:Nmaxn])
                    }
                }
                if(vn$B > 0){# binary
                    for(v in 1:Bn){
                        Bprob[v, k] ~ dbeta(shape1=Bshapelo, shape2=Bshapehi)
                    }
                }
            }
            ## Probability of data
            for(d in 1:npoints){
                K[d] ~ dcat(prob=W[1:nclusters])
                ##
                if(vn$R > 0){# continuous
                    for(v in 1:Rn){
                        Rdata[d, v] ~ dnorm(mean=Rmean[v, K[d]], var=Rvar[v, K[d]])
                    }
                }
                if(vn$C > 0){# censored
                    for(v in 1:Cn){
                        Caux[d, v] ~ dconstraint(Clat[d, v] >= Cleft[d, v] & Clat[d, v] <= Cright[d, v])
                        Clat[d, v] ~ dnorm(mean=Cmean[v, K[d]], var=Cvar[v, K[d]])
                    }
                }
                if(vn$D > 0){# discretized
                    for(v in 1:Dn){
                        Daux[d, v] ~ dconstraint(Dlat[d, v] >= Dleft[d, v] & Dlat[d, v] < Dright[d, v])
                        Dlat[d, v] ~ dnorm(mean=Dmean[v, K[d]], var=Dvar[v, K[d]])
                    }
                }
                if(vn$O > 0){# ordinal
                    for(v in 1:On){
                        Oaux[d, v] ~ dconstraint(Olat[d, v] >= Oleft[d, v] & Olat[d, v] < Oright[d, v])
                        Olat[d, v] ~ dnorm(mean=Omean[v, K[d]], var=Ovar[v, K[d]])
                    }
                }
                if(vn$N > 0){# nominal
                    for(v in 1:Nn){
                        Ndata[d, v] ~ dcat(prob=Nprob[v, K[d], 1:Nmaxn])
                    }
                }
                if(vn$B > 0){# binary
                    for(v in 1:Bn){
                        Bdata[d, v] ~ dbern(prob=Bprob[v, K[d]])
                    }
                }
            }
        })


#### INITIAL-VALUE FUNCTION
        initsfn <- function(){
            Alpha <- sample(1:nalpha, 1, prob=constants$probalpha0, replace=T)
            W <- 1/nclusters + 0*nimble::rdirch(n=1, alpha=constants$basealphas*2^Alpha)
            outlist <- list(
                Alpha = Alpha,
                W = W,
                ## K = rep(which(W>0)[1], npoints)
                ## K = sample(rep(which(W>0),2), npoints, replace=T)
                K = rep(sample(rep(which(W>0), 2), 1, replace=T), npoints)
            )
            ##
            if(vn$R > 0){# continuous
                Rrate <- matrix(nimble::rinvgamma(n=vn$R*nclusters, shape=constants$Rshapehi, rate=constants$Rvar1), nrow=vn$R, ncol=nclusters)
                outlist <- c(outlist,
                             list(
                                 Rmean = matrix(rnorm(n=vn$R*nclusters, mean=constants$Rmean1, sd=sqrt(constants$Rvarm1)), nrow=vn$R, ncol=nclusters),
                                 Rrate = Rrate,
                                 Rvar = matrix(nimble::rinvgamma(n=vn$R*nclusters, shape=constants$Rshapelo, rate=Rrate), nrow=vn$R, ncol=nclusters)
                             ))
            }
            if(vn$C > 0){# censored
                Crate <- matrix(nimble::rinvgamma(n=vn$C*nclusters, shape=constants$Cshapehi, rate=constants$Cvar1), nrow=vn$C, ncol=nclusters)
                outlist <- c(outlist,
                             list(
                                 Cmean = matrix(rnorm(n=vn$C*nclusters, mean=constants$Cmean1, sd=sqrt(constants$Cvarm1)), nrow=vn$C, ncol=nclusters),
                                 Crate = Crate,
                                 Cvar = matrix(nimble::rinvgamma(n=vn$C*nclusters, shape=constants$Cshapelo, rate=Crate), nrow=vn$C, ncol=nclusters),
                                 Clat = vtransform(data[,vnames$C, with=F], auxmetadata, Cout='init') ## for data with boundary values
                             ))
            }
            if(vn$D > 0){# discretized
                Drate <- matrix(nimble::rinvgamma(n=vn$D*nclusters, shape=constants$Dshapehi, rate=constants$Dvar1), nrow=vn$D, ncol=nclusters)
                outlist <- c(outlist,
                             list(
                                 Dmean = matrix(rnorm(n=vn$D*nclusters, mean=constants$Dmean1, sd=sqrt(constants$Dvarm1)), nrow=vn$D, ncol=nclusters),
                                 Drate = Drate,
                                 Dvar = matrix(nimble::rinvgamma(n=vn$D*nclusters, shape=constants$Dshapelo, rate=Drate), nrow=vn$D, ncol=nclusters),
                                 Dlat = vtransform(data[,vnames$D, with=F], auxmetadata, Dout='init') ## for data with boundary values
                             ))
            }
            if(vn$O > 0){# ordinal
                Orate <- matrix(nimble::rinvgamma(n=vn$O*nclusters, shape=constants$Oshapehi, rate=constants$Ovar1), nrow=vn$O, ncol=nclusters)
                outlist <- c(outlist,
                             list(
                                 Omean = matrix(rnorm(n=vn$O*nclusters, mean=constants$Omean1, sd=sqrt(constants$Ovarm1)), nrow=vn$O, ncol=nclusters),
                                 Orate = Orate,
                                 Ovar = matrix(nimble::rinvgamma(n=vn$O*nclusters, shape=constants$Oshapelo, rate=Orate), nrow=vn$O, ncol=nclusters),
                                 Olat = vtransform(data[,vnames$O, with=F], auxmetadata, Oout='init', useOquantiles=useOquantiles) ## for data with boundary values
                             ))
            }
            if(vn$N > 0){# nominal
                outlist <- c(outlist,
                             list(
                                 Nprob = aperm(array(sapply(1:vn$N, function(avar){sapply(1:nclusters, function(aclus){rdirch(n=1, alpha=Nalpha0[avar,])})}), dim=c(Nmaxn,nclusters,vn$N)))
                             ))
            }
            if(vn$B > 0){# binary
                outlist <- c(outlist,
                             list(
                                 Bprob = matrix(rbeta(n=vn$B*nclusters, shape1=constants$Bshapelo, shape2=constants$Bshapehi), nrow=vn$B, ncol=nclusters)
                             ))
            }
            ##
            outlist
        }

        timecount <- Sys.time()

##################################################
#### NIMBLE SETUP
##################################################

        finitemixnimble <- nimbleModel(
            code=finitemix, name='finitemixnimble1',
            constants=constants,
            data=datapoints,
            inits=initsfn()
        )

        Cfinitemixnimble <- compileNimble(finitemixnimble, showCompilerOutput=FALSE)
        gc()

        confnimble <- configureMCMC(
            Cfinitemixnimble, #nodes=NULL,
            monitors=c('W',
                       if(vn$R > 0){c('Rmean', 'Rvar')},
                       if(vn$C > 0){c('Cmean', 'Cvar')},
                       if(vn$D > 0){c('Dmean', 'Dvar')},
                       if(vn$O > 0){c('Omean', 'Ovar')},
                       if(vn$N > 0){c('Nprob')},
                       if(vn$B > 0){c('Bprob')}
                       ),
            monitors2=c(if(showclusterstraces){'Alpha'}, 'K')
        )
        ## print(confnimble$getUnsampledNodes())
        ## confnimble$printSamplers(executionOrder=TRUE)
        
        targetslist <- sapply(confnimble$getSamplers(), function(xx)xx$target)
        nameslist <- sapply(confnimble$getSamplers(), function(xx)xx$name)
        ## cat('\n******** NAMESLIST',nameslist,'\n')

        ## replace Alpha's cat-sampler with slice
        if(Alphatoslice && !('Alpha' %in% targetslist[nameslist == 'posterior_predictive'])){
            confnimble$removeSamplers('Alpha')
            confnimble$addSampler(target='Alpha', type='slice')
        }

        ## replace K's cat-sampler with slice
        if(Ktoslice){
            for(asampler in grep('^K\\[', targetslist,value=T)){
                if(!(asampler %in% targetslist[nameslist == 'posterior_predictive'])){
                    confnimble$removeSamplers(asampler)
                    confnimble$addSampler(target=asampler, type='slice')
                }
            }
        }

        ## replace all RW samplers with slice
        testreptime <- Sys.time()
        if(RWtoslice){
            for(asampler in targetslist[nameslist == 'RW']){
                confnimble$removeSamplers(asampler)
                confnimble$addSampler(target=asampler, type='slice')
                ## confnimble$replaceSampler(target=asampler, type='slice')
            }
        }

        ## print(confnimble$getUnsampledNodes())
        ## call this to do a first reordering
        mcsampler <- buildMCMC(confnimble)
        ## print(confnimble$getUnsampledNodes())

#### change execution order for some variates
        samplerorder <- c(
            'K',
            if(vn$R > 0){c('Rmean','Rrate','Rvar')},
            if(vn$C > 0){c('Cmean','Crate','Cvar')},
            if(vn$D > 0){c('Dmean','Drate','Dvar')},
            if(vn$O > 0){c('Omean','Orate','Ovar')},
            if(vn$N > 0){c('Nprob')},
            if(vn$B > 0){c('Bprob')},
            'W','Alpha'
        )
        ##
        neworder <- foreach(var=samplerorder, .combine=c)%do%{grep(paste0('^',var,'(\\[.+\\])*$'), sapply(confnimble$getSamplers(), function(x){
            if(!(x$name == 'posterior_predictive')){ x$target }else{NULL}
        }))}
        ##
        ## cat('\n********NEW ORDER',neworder,'\n')
        confnimble$setSamplerExecutionOrder(c(setdiff(confnimble$getSamplerExecutionOrder(), neworder), neworder))
        print(confnimble)
        ## print(confnimble$getUnsampledNodes())

#### Compile Monte Carlo sampler
        mcsampler <- buildMCMC(confnimble)
        Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

        
        cat('\nSetup time', printtime(Sys.time() - timecount), '\n')

        if(acore == 1){
            printnull('Done. Estimating remaining time, please be patient...', outcon)
        }

##################################################
#### LOOP OVER CHAINS (WITHIN ONE CORE)
##################################################
        maxusedclusters <- 0
        ## build test data for assessing stationarity
        testdata <- data.table(sapply(1:nrow(auxmetadata),function(xx){
            xx <- auxmetadata[xx,]
            toadd <- xx[,paste0('mctest',1:3),with=F]
            if(xx[['mcmctype']] %in% c('B','N')){
                toadd <- xx[1,paste0('V',toadd),with=F]
            }
            toadd
        }))
        colnames(testdata) <- auxmetadata[['name']]
        ## testdata <- t(auxmetadata[,paste0('mctest',1:3)])
        ## colnames(testdata) <- auxmetadata[,name]
        ## testdata <- rbind(auxmetadata[['Q2']], varinfo[['Q1']], varinfo[['Q3']]) #, varinfo[['plotmin']], varinfo[['plotmax']], varinfo[['datamin']], varinfo[['datamax']])

        starttime <- Sys.time()
#### LOOP
        allflagmc <- FALSE

        for(achain in 1:nchainspercore){
            showsamplertimes0 <- showsamplertimes && (achain==1)
            ## showclusterstraces0 <- showclusterstraces && (achain==1)
            niter <- niterini
            ## ## Experimental: decrease number of iterations based on previous chain
            ## if(achain == 1){
            ##     niter <- niterini
            ## }else{
            ##     niter <- max(min(niterini,lengthmeasure*2), 128)
            ##     }
            nitertot <- 0
            lengthmeasure <- +Inf
            reset <- TRUE
            traces <- NULL
            ## allmcsamples <- NULL
            allmcsamples <- NULL
            allclusterhypar <- list(Alpha=NULL, K=NULL)
            if(showclusterstraces){
                clusterhypar <- NULL
            }
            flagll <- FALSE
            flagmc <- FALSE
            gc()
            chainnumber <- (acore-1L)*nchainspercore + achain
            cat('\nChain #', chainnumber,
                '(chain', achain,'of',nchainspercore,'for this core)\n')
            cat('Seed:', chainnumber+seed, '\n')
            set.seed(chainnumber+seed)
            Cfinitemixnimble$setInits(initsfn())

#### WHILE LOOP CONTINUING UNTIL CONVERGENCE
            while(nitertot < lengthmeasure){
                cat('Iterations:', niter,'\n')
                
#### MONTE-CARLO CALL
                Cmcsampler$run(niter=niter, thin=1, thin2=(if(showclusterstraces){max(1,round(niter/nclustersamples))}else{niter}), nburnin=0, time=showsamplertimes0, reset=reset, resetMV=TRUE)
                
                ## mcsamples <- as.matrix(Cmcsampler$mvSamples)
                mcsamples <- as.list(Cmcsampler$mvSamples, iterationAsLastIndex=T)

                clusterhypar <- as.list(Cmcsampler$mvSamples2, iterationAsLastIndex=F)
                clusterhypar$K <- apply(clusterhypar$K, 1, function(xx)length(unique(xx)))
                clusterhypar$Alpha <- c(clusterhypar$Alpha)

                cat('\nMCMC time', printtime(Sys.time() - starttime), '\n')

                ## #### Remove iterations with non-finite values
                ##                 toremove <- which(!is.finite(mcsamples), arr.ind=T)
                ##                 ##
                ##                 if(length(toremove) > 0){
                ##                     printnull('\nWARNING: SOME NON-FINITE OUTPUTS\n', outcon)
                ##                     ##
                ##                     flagmc <- TRUE
                ##                     allflagmc <- TRUE
                ##                     if(length(unique(toremove[,1])) == nrow(mcsamples)){
                ##                         suppressWarnings(sink())
                ##                         suppressWarnings(sink(NULL,type='message'))
                ##                         registerDoSEQ()
                ##                         if(ncores > 1){ stopCluster(cl) }
                ##                         stop('...TOO MANY NON-FINITE OUTPUTS. ABORTING')
                ##                     }else{
                ##                         mcsamples <- mcsamples[-unique(toremove[,1]),,drop=F]
                ##                     }
                ##                 }
#### Remove iterations with non-finite values
                toremove <- sort(unique(unlist(lapply(mcsamples,function(xx){temp <- which(is.na(xx),arr.ind=T);temp[,ncol(temp)]}))))
                if(length(toremove) > 0){
                    cat('\nWARNING:',length(toremove),'NON-FINITE SAMPLES\n')
                    ##
                    flagmc <- TRUE
                    allflagmc <- TRUE
                    if(length(unique(toremove[,1])) == ncol(mcsamples$W)){
                        suppressWarnings(sink())
                        suppressWarnings(sink(NULL,type='message'))
                        ## registerDoSEQ()
                        if(exists('cl')){ stopCluster(cl) }
                        stop('...TOO MANY NON-FINITE OUTPUTS. ABORTING')
                    }else{
                        mcsamples <- mcsubset(mcsamples, -toremove)
                    }
                }

                
                ##
                if(showsamplertimes0){
                    samplertimes <- Cmcsampler$getTimes()
                    names(samplertimes) <- sapply(confnimble$getSamplers(),function(x)x$target)
                    sprefixes <- unique(sub('^([^[]+)(\\[.*\\])', '\\1', names(samplertimes)))
                    cat('\nSampler times:\n')
                    print(sort(sapply(sprefixes, function(x)sum(samplertimes[grepl(x,names(samplertimes))])),decreasing=T))
                }
                ##
                
                ## Check how many clusters were occupied. Warns if too many
                usedclusters <- clusterhypar$K[length(clusterhypar$K)]
                ## maxusedclusters <- max(usedclusters, maxusedclusters)
                ## if(usedclusters > nclusters-5){
                ##     cat('\nWARNING: TOO MANY CLUSTERS OCCUPIED')
                ##     printnull('\nWARNING: TOO MANY CLUSTERS OCCUPIED\n', outcon)
                ## }
                cat('\nOCCUPIED CLUSTERS:', usedclusters, 'OF', nclusters,'\n')
                
                ##
                ## Diagnostics
                ## Log-likelihood
                diagntime <- Sys.time()
                ##
                ll <- t(
                    log(samplesFDistribution(Y=testdata, X=NULL, mcsamples=mcsamples, auxmetadata=auxmetadata, jacobian=FALSE, useOquantiles=useOquantiles, parallel=FALSE, silent=TRUE)) #- sum(log(invjacobian(data.matrix(data0), varinfo)), na.rm=T)
                )
                colnames(ll) <- paste0('log-',c('mid','lo','hi')) #,'pm','pM','dm','dM'))
                ##
                traces <- rbind(traces, 10/log(10) * ll)
                ## pdff('testdifferenttraces2')
                ## tplot(y=traces, lwd=1, lty=1)
                ## for(i in 1:ncol(traces)){
                ##     tplot(y=traces[,i], main=colnames(traces)[i])
                ## }
                ## dev.off()
                toremove <- which(!is.finite(traces), arr.ind=T)
                if(length(toremove) > 0){
                    flagll <- TRUE
                    traces <- traces[-unique(toremove[,1]),,drop=F]
                }
                ##
                diagnESS <- LaplacesDemon::ESS(traces)
                cat('\nESSs:',paste0(round(diagnESS), collapse=', '))
                diagnIAT <- apply(traces, 2, function(x){LaplacesDemon::IAT(x)})
                cat('\nIATs:',paste0(round(diagnIAT), collapse=', '))
                diagnBMK <- LaplacesDemon::BMK.Diagnostic(traces[1:(4*trunc(nrow(traces)/4)),], batches=4)[,1]
                cat('\nBMKs:',paste0(round(diagnBMK,3), collapse=', '))
                diagnMCSE <- 100*apply(traces, 2, function(x){funMCSE(x)/sd(x)})
                cat('\nMCSEs:',paste0(round(diagnMCSE,2), collapse=', '))
                diagnStat <- apply(traces, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
                cat('\nStationary:',paste0(diagnStat, collapse=', '))
                diagnBurn <- apply(traces, 2, function(x){LaplacesDemon::burnin(matrix(x[1:(10*trunc(length(x)/10))], ncol=1))})
                cat('\nBurn-in I:',paste0(diagnBurn, collapse=', '))
                diagnBurn2 <- proposeburnin(traces, batches=10)
                cat('\nBurn-in II:',diagnBurn2)
                diagnThin <- proposethinning(traces)
                cat('\nProposed thinning:',paste0(diagnThin, collapse=', '),'\n')

                cat('\nDiagnostics time', printtime(Sys.time() - diagntime), '\n')


                ## allmcsamples <- rbind(allmcsamples, mcsamples)
                ## rm(mcsamples)
                ## gc()
                
                if(is.null(allmcsamples)){
                    allmcsamples <- mcsamples
                }else{
                    allmcsamples <- mapply(function(xx,yy){
                        temp <- c(xx,yy)
                        dx <- dim(xx)[-length(dim(xx))]
                        dim(temp) <- c(dx, length(temp)/prod(dx))
                        temp
                    },
                    allmcsamples, mcsamples)
                }
                rm(mcsamples)
                gc()

                if(showclusterstraces){
                    allclusterhypar <- mapply(function(xx,yy){c(xx,yy)}, allclusterhypar, clusterhypar, SIMPLIFY=F)
                }
                nitertot <- ncol(allmcsamples$W)
                
#########################################
#### CHECK IF CHAIN MUST BE CONTINUED ####
#########################################
                lengthmeasure <- thresholdfn(diagnESS=diagnESS, diagnIAT=diagnIAT, diagnBMK=diagnBMK, diagnMCSE=diagnMCSE, diagnStat=diagnStat, diagnBurn=diagnBurn, diagnBurn2=diagnBurn2, diagnThin=diagnThin)
                cat('\nNumber of iterations', nitertot, ', required', lengthmeasure,'\n')
                ##
                if(nitertot < lengthmeasure){
                    ## limit number of iterations per loop, to save memory
                    niter <- min(lengthmeasure - nitertot + 1L, niterini)
                    cat('Increasing by', niter, '\n')
                }
                reset <- FALSE
            }

            
#########################################
#### SAVE CHAIN ####
#########################################
            ## tokeep <- seq(to=nrow(allmcsamples), length.out=nsamplesperchain, by=max(1,multcorr*ceiling(max(diagnIAT,diagnThin)), na.rm=T))
            ## allmcsamples <- allmcsamples[tokeep,,drop=F]
            ## ##
            ## saveRDS(allmcsamples, file=paste0(dirname,'_mcsamples-',nameroot,'--',chainnumber,'.rds'))
            ## ## rm(allmcsamples)

            tokeep <- seq(to=ncol(allmcsamples$W), length.out=nsamplesperchain, by=max(1,multcorr*ceiling(max(diagnIAT,diagnThin)), na.rm=T))
            ##
            saveRDS(mcsubset(allmcsamples, tokeep), file=paste0(dirname,'_mcsamples-',nameroot,'--',chainnumber,'.rds'))
            rm(allmcsamples)
            ## nitertot <- ncol(allmcsamples$W)

            gc()

            saveRDS(traces[tokeep,],file=paste0(dirname,'_mctraces-',nameroot,'--',chainnumber,'.rds'))

            for(i in 1:length(allclusterhypar))
                ## Check how many clusters were occupied. Warns if too many
                usedclusters <- allclusterhypar$K[length(allclusterhypar$K)]
            if(usedclusters > nclusters-5){
                cat('\nWARNING: TOO MANY CLUSTERS OCCUPIED')
                ## printnull('\nWARNING: TOO MANY CLUSTERS OCCUPIED\n', outcon)
                ## printnull(paste0('\nOCCUPIED CLUSTERS:', usedclusters, 'OF', nclusters,'\n'), outcon)
            }
            cat('\nOCCUPIED CLUSTERS:', usedclusters, 'OF', nclusters,'\n')

###############
#### PLOTS ####
###############
#### Plot Alpha and cluster occupation, if required
            if(showclusterstraces){
                cat('Plotting Alpha and cluster information.\n')
                cat('\nSTATS OCCUPIED CLUSTERS:\n')
                print(summary(allclusterhypar$K))
                ##
                pdff(paste0(dirname,'_hyperparams_traces-',nameroot,'--',chainnumber,'_',achain,'-',acore), apaper=4)
                tplot(y=allclusterhypar$K, ylab='occupied clusters',xlab='iteration',ylim=c(0,nclusters))
                tplot(x=((-1):nclusters)+0.5,y=tabulate(allclusterhypar$K + 1, nbins=nclusters+1), type='h', xlab='occupied clusters', ylab=NA, ylim=c(0,NA))
                ##
                cat('\nSTATS log2(alpha):\n')
                print(summary(allclusterhypar$Alpha+minalpha-1L, na.rm=T))
                tplot(y=allclusterhypar$Alpha+minalpha-1L, ylab=bquote(log2(alpha)),xlab='iteration',ylim=c(minalpha,maxalpha))
                tplot(x=(minalpha:(maxalpha+1))-0.5,y=tabulate(allclusterhypar$Alpha,nbin=nalpha), type='h', xlab=bquote(log2(alpha)), ylab='', ylim=c(0,NA))
                dev.off()
            }
            
#### Plot diagnostic traces of current chain
            if(plottraces){
                cat('\nPlotting traces and samples.\n')

                tracegroups <- list(1,2,3)
                names(tracegroups) <- colnames(traces)
                grouplegends <- foreach(agroup=1:length(tracegroups))%do%{
                    c( paste0('-- STATS ', names(tracegroups)[agroup], ' --'),
                      paste0('min ESS = ', signif(min(diagnESS[tracegroups[[agroup]]]),6)),
                      paste0('max IAT = ', signif(max(diagnIAT[tracegroups[[agroup]]]),6)),
                      paste0('max BMK = ', signif(max(diagnBMK[tracegroups[[agroup]]]),6)),
                      paste0('max MCSE = ', signif(max(diagnMCSE[tracegroups[[agroup]]]),6)),
                      paste0('stationary: ', sum(diagnStat[tracegroups[[agroup]]]),'/',length(diagnStat[tracegroups[[agroup]]])),
                      ## paste0('burn: ', signif(max(diagnBurn[tracegroups[[agroup]]]),6))
                      paste0('burn: ', signif(diagnBurn2,6)),
                      paste0('max thin = ', signif(max(diagnThin[tracegroups[[agroup]]]),6))
                      )
                }
                colpalette <- c(7,2,1)
                names(colpalette) <- colnames(traces)
                ##
                ## Plot various info and traces
                cat('\nPlotting MCMC traces')
                graphics.off()
                pdff(paste0(dirname,'_mcmcpartialtraces-',nameroot,'--',chainnumber,'_',achain,'-',acore), apaper=4)
                ## Summary stats
                matplot(1:2, type='l', col='white', main=paste0('Stats chain ',achain), axes=FALSE, ann=FALSE)
                legendpositions <- c('topleft','topright','bottomleft','bottomright')
                for(alegend in 1:length(grouplegends)){
                    legend(x=legendpositions[alegend], bty='n', cex=1.5,
                           legend=grouplegends[[alegend]] )
                }
                legend(x='center', bty='n', cex=1,
                       legend=c(
                           paste0('Chain ', chainnumber,'_',achain,'-',acore),
                           paste0('Occupied clusters: ', usedclusters, ' of ', nclusters),
                           paste0('LL:  ( ', signif(mean(traces[,1]),3), ' +- ', signif(sd(traces[,1]),3),' ) dHart'),
                           'NOTES:',
                           if(flagmc){'some non-finite MC outputs'},
                           if(usedclusters > nclusters-5){'too many clusters occupied'},
                           if(flagll){'non-finite values in diagnostics'}
                       ))
                ##
                ## Traces of likelihood and cond. probabilities
                par(mfrow=c(1,1))
                for(avar in colnames(traces)){
                    tplot(y=traces[,avar], type='l', lty=1, col=colpalette[avar],
                          main=paste0( 'ESS = ', signif(diagnESS[avar], 3),
                                      ' | IAT = ', signif(diagnIAT[avar], 3),
                                      ' | BMK = ', signif(diagnBMK[avar], 3),
                                      ' | MCSE = ', signif(diagnMCSE[avar], 3),
                                      ' | stat: ', diagnStat[avar]*1L,
                                      ' | burnI: ', diagnBurn[avar],
                                      ' | burnII: ', diagnBurn2,
                                      ' | thin: ', diagnThin[avar]
                                      ), 
                          ylab=paste0(avar,'/dHart'),
                          xlab='Monte Carlo sample',
                          family=family, mar=c(NA,6,NA,NA) )
                }
                dev.off()
            }
            
            ##
## #### Plot samples from current chain
##             if(plotpartialsamples){
##                 ## Samples of marginal frequency distributions
##                 ## if(!continue){
##                 subsamples <- (if(totsamples=='all'){1:nitertot}else{round(seq(1, nitertot, length.out=totsamples))})
##                 ## showsubsample <- round(seq(1, length(subsamples), length.out=showsamples))
##                 ##
##                 cat('\nPlotting samples of frequency distributions')
##                 plotFsamples(file=paste0(dirname,'mcmcdistributions-',nameroot,'--',chainnumber,'_',achain,'-',acore),
##                              mcsamples=mcsubset(allmcsamples,subsamples),
##                              auxmetadata=auxmetadata,
##                              data=data,
##                              plotuncertainty='samples',
##                              uncertainty=showsamples,
##                              plotmeans=plotmeans,
##                              datahistogram=TRUE, datascatter=TRUE,
##                              useOquantiles=useOquantiles,
##                              parallel=FALSE, silent=TRUE
##                              )
##             }

            cat('\nMCMC+diagnostics time', printtime(Sys.time() - starttime), '\n')
            

#### Print estimated remaining time
            ertime <- (Sys.time()-starttime)/achain*(nchainspercore-achain+1)
            if(is.finite(ertime) && ertime > 0){
                printnull(paste0('\rSampling. Estimated remaining time ',
                                 printtime(ertime),
                                 '                  '),
                          outcon)
            }

            maxusedclusters <- max(maxusedclusters, usedclusters)
        } #### END LOOP OVER CHAINS (WITHIN ONE CORE)

        ##
        cat('\nTotal time', printtime(Sys.time() - starttime), '\n')

        cbind(maxusedclusters, allflagmc)
    } #### END FOREACH OVER CORES
    suppressWarnings(sink())
    suppressWarnings(sink(NULL,type='message'))

    maxusedclusters <- max(chaininfo[,1])
    nonfinitechains <- sum(chaininfo[,2])
    
############################################################
#### End of all MCMC
############################################################

    
############################################################
#### Join chains
############################################################
    ## mcsamples <- foreach(chainnumber=1:(ncores*nchainspercore), .combine=rbind)%do%{
    ##     readRDS(file=paste0(dirname,'_mcsamples-',nameroot,'--',chainnumber,'.rds'))
    ## }
    ## ## attr(mcsamples, 'rng') <- NULL
    ## ## attr(mcsamples, 'doRNG_version') <- NULL
    ## ## ## Remove extra chains
    ## ## mcsamples <- mcsamples[round(seq(1,nrow(mcsamples),length.out=nsamples)),-(1:3)]
    ## saveRDS(mcsamples,file=paste0(dirname,'Fdistribution-',nameroot,'.rds'))

    joinmc <- function(mc1, mc2){
        mapply(function(xx,yy){
            temp <- c(xx,yy)
            dx <- dim(xx)[-length(dim(xx))]
            dim(temp) <- c(dx, length(temp)/prod(dx))
            temp
        },
        mc1, mc2)
    }
    mcsamples <- foreach(chainnumber=1:(ncores*nchainspercore), .combine=joinmc, .multicombine=F)%do%{
        readRDS(file=paste0(dirname,'_mcsamples-',nameroot,'--',chainnumber,'.rds'))
    }
    
    saveRDS(mcsamples,file=paste0(dirname,'Fdistribution-',nameroot,'.rds'))

    traces <- foreach(chainnumber=1:(ncores*nchainspercore), .combine=rbind)%do%{
        readRDS(file=paste0(dirname,'_mctraces-',nameroot,'--',chainnumber,'.rds'))
    }
    ## traces <- mcsamples[round(seq(1,nrow(mcsamples),length.out=nsamples)),1:3]
    saveRDS(traces,file=paste0(dirname,'MCtraces-',nameroot,'.rds'))

    cat('\rFinished Monte Carlo sampling.\n')
    gc()

############################################################
#### Final joint diagnostics
############################################################
    ## cat('\nSome diagnostics:\n')
    traces <- traces[apply(traces,1,function(x){all(is.finite(x))}),]
    ## flagll <- nrow(traces) != nrow(traces)

    ## funMCSE <- function(x){LaplacesDemon::MCSE(x, method='batch.means')$se}
    ## diagnESS <- LaplacesDemon::ESS(traces)
    ## diagnIAT <- apply(traces, 2, function(x){LaplacesDemon::IAT(x)})
    ## diagnBMK <- LaplacesDemon::BMK.Diagnostic(traces[1:(4*trunc(nrow(traces)/4)),], batches=4)[,1]
    ## diagnMCSE <- 100*apply(traces, 2, function(x){funMCSE(x)/sd(x)})
    ## diagnStat <- apply(traces, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
    ## diagnBurn <- apply(traces, 2, function(x){LaplacesDemon::burnin(matrix(x[1:(10*trunc(length(x)/10))], ncol=1))})
    ## diagnBurn2 <- proposeburnin(traces, batches=10)
    ## diagnThin <- proposethinning(traces)
    ## ##
    ## cat('\nESSs:',paste0(round(diagnESS), collapse=', '))
    ## cat('\nIATs:',paste0(round(diagnIAT), collapse=', '))
    ## cat('\nBMKs:',paste0(round(diagnBMK,3), collapse=', '))
    ## cat('\nMCSEs:',paste0(round(diagnMCSE,2), collapse=', '))
    ## cat('\nStationary:',paste0(diagnStat, collapse=', '))
    ## cat('\nBurn-in I:',paste0(diagnBurn, collapse=', '))
    ## cat('\nBurn-in II:',diagnBurn2)
    ## cat('\n')
    ## ## cat(paste0('\nProposed thinning: ',paste0(diagnThin, collapse=', ')))
    ## ##       )
    ## ##         }
    ## ##

    ## Plot various info and traces
    cat('\nPlotting final Monte Carlo traces.\n')

    ##
    colpalette <- c(7,2,1)
    names(colpalette) <- colnames(traces)
    graphics.off()
    pdff(paste0(dirname,'MCtraces-',nameroot), apaper=4)
    ## Traces of likelihood and cond. probabilities
    for(avar in colnames(traces)){
        tplot(y=traces[,avar], type='l', lty=1, col=colpalette[avar],
              ## main=paste0(avar,
              ##             '\nESS = ', signif(diagnESS[avar], 3),
              ##             ' | IAT = ', signif(diagnIAT[avar], 3),
              ##             ' | BMK = ', signif(diagnBMK[avar], 3),
              ##             ' | MCSE = ', signif(diagnMCSE[avar], 3),
              ##             ' | stat: ', diagnStat[avar],
              ##             ' | burn I: ', diagnBurn[avar],
              ##             ' | burn II: ', diagnBurn2
              ##             ),
              ylab=paste0(avar,'/dHart'), xlab='sample', family=family
              )
    }

    cat('\nMax number of occupied clusters:',maxusedclusters,'\n')
    if(maxusedclusters > nclusters-5){
        cat('Too many clusters occupied\n')
    }

    if(nonfinitechains > 0){
        cat(nonfinitechains,'chains with some non-finite outputs\n')
    }

    cat('Plotting marginal samples.\n')
    plotFsamples(file=paste0(dirname,'plotsamples_Fdistribution-',nameroot),
                 mcsamples=mcsamples, auxmetadata=auxmetadata,
                 data=data,
                 plotuncertainty='samples',
                 uncertainty=showsamples, plotmeans=TRUE,
                 datahistogram=TRUE, datascatter=TRUE,
                 useOquantiles=useOquantiles,
                 parallel=TRUE, silent=TRUE
                 )

    cat('Plotting marginal samples with quantiles.\n')
    plotFsamples(file=paste0(dirname,'plotquantiles_Fdistribution-',nameroot),
                 mcsamples=mcsamples, auxmetadata=auxmetadata,
                 data=data,
                 plotuncertainty='quantiles',
                 uncertainty=showquantiles, plotmeans=TRUE,
                 datahistogram=TRUE, datascatter=TRUE,
                 useOquantiles=useOquantiles,
                 parallel=TRUE, silent=TRUE
                 )

    cat('\nTotal computation time', printtime(Sys.time() - timestart0), '\n')

    if(exists('cl')){
        cat('\nClosing connections to cores.\n')
        stopCluster(cl)
    }

    
#### remove partial files if required
    if(cleanup){
        cat('Removing temporary output files.\n')
        for(chainnumber in 1:(ncores*nchainspercore)){
            ## file.remove(paste0(dirname,'_mcsamples-',nameroot,'--',chainnumber,'.rds'))
            file.remove(paste0(dirname,'_mcsamples-',nameroot,'--',chainnumber,'.rds'))
            file.remove(paste0(dirname,'_mctraces-',nameroot,'--',chainnumber,'.rds'))
        }
    }
    
    cat('Finished.\n\n')

    if(output){mcsamples}
}
