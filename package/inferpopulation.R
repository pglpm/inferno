inferpopulation <- function(dataset, varinfoaux, outputdir, nsamples=4096, nsamplesperchain=4, nchains, ncores, saveallchains=F, plotallchains=F, seed=701){

    if(!missing(nsamples) && !missing(nchains) && missing(nsamplesperchain)){
        nsamplesperchain <- ceiling(nsamples/nchains)
    }else if(!missing(nsamples) && missing(nchains) && !missing(nsamplesperchain)){
        nchains <- ceiling(nsamples/nsamplesperchain)
    }else if(missing(nsamples) && !missing(nchains) && !missing(nsamplesperchain)){
        nsamples <- nchains * nsamplesperchain
    }else{
        stop('please specify exactly two among "nsamples", "nchains", "nsamplesperchain"')
    }

    ## set number of cores for parallel computation
    if(missing(ncores) || !is.numeric(ncores)){
        warning('The number of cores has not been given.\nIt is much preferable that it be set by the user.')
        ncores <- round(parallel::detectCores()/2)
        cat('\nTrying to use ',ncores,' cores\n')
    }

    if(missing(outputdir)){
        outputdir <- paste0('_inference_',format(Sys.time(), '%y%m%dT%H%M%S'))
    }

    nchainspercore <- ceiling(nchains/ncores)

    timestart0 <- Sys.time()

    cat('\n')
    require('data.table')
    require('LaplacesDemon', include.only=NULL)
    require('foreach')
    require('doParallel')
    require('doRNG')
    

    ## read dataset
    if(is.character(dataset) && file.exists(dataset)){
        dataset <- paste0(sub('.csv$', '', dataset), '.csv')
        dataset <- fread(dataset, na.strings='')
    }
    dataset <- as.data.table(dataset)

    ## varinfoaux
    if(is.character(varinfoaux) && file.exists(varinfoaux)){
        varinfoaux <- paste0(sub('.rds$', '', varinfoaux), '.rds')
        varinfoaux <- readRDS(varinfoaux)
    }

    ## ## list of predictand variates
    ## if(is.character(predictands) && file.exists(predictands)){
    ##     predictands <- as.vector(unlist(read.csv(predictands, header=F)))
    ## }else{
    ##     predictands <- unlist(predictands)
    ## }
    
    ## file to save output
    if(is.character(file) || (is.logical(file) && file)){ # must save to file
        if(is.character(file)){
            file <- paste0(sub('.rds$', '', file), '.rds')
        }else{
            if(file.exists('longrunsamples.rds')){
                file <- paste0('longrunsamples_',format(Sys.time(), '%y%m%dT%H%M%S'),'.rds')
            }else{
                file <- 'longrunsamples.rds'
            }
        }
    }

#### various internal parameters
    niter0 <- 1024L # initial iterations to try
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
    npoints <- nrow(dataset)

    ## other options
    showdata <- TRUE # 'histogram' 'scatter' FALSE TRUE
    plotmeans <- TRUE # plot frequency averages
    totsamples <- 'all' # 'all' number of samples if plotting frequency averages
    showsamples <- 100 # number of samples to show. Shown separately for posterior=F
    showhyperparametertraces <- FALSE ##
    showsamplertimes <- FALSE ##
    family <- 'Palatino'

    basename <- paste0(outputdir,'-V',nrow(varinfoaux),'-D',npoints,'-K',nclusters,'-I',nsamples)
    ##
    dirname <- paste0(basename,'/')
    dir.create(dirname)
    cat('\n',paste0(rep('*',max(nchar(dirname),26)),collapse=''),
        '\n Saving output in directory\n',dirname,'\n',
       paste0(rep('*',max(nchar(dirname),26)),collapse=''),'\n')


    ## Parameter and function to test MCMC convergence
    multcorr <- 2L

        source('pglpm_plotfunctions.R')
        source('vtransform.R')
        source('samplesFDistribution.R')
        source('proposeburnin.R')
        source('proposethinning.R')
        source('plotFsamples.R')

    vn <- list()
    vnames <- list()
    for(atype in c('R','C','D','O','N','B')){
        vn[[atype]] <- length(varinfoaux[mcmctype == atype, name])
        vnames[[atype]] <- varinfoaux[mcmctype == atype, name]
    }
    ## ## choose predictands if not chosen
    ## if(is.null(predictands) || !exists('predictands')){
    ##     predictands <- sapply(vnames,function(xx)xx[1])
    ##     cat('\nSelf-choosing predictand variates:\n', predictands, '\n')
    ## }

    if(vn$N > 0){
        Nmaxn <- max(varinfoaux[mcmctype == 'N', Nvalues])
        Nalpha0 <- matrix(1e-100, nrow=vn$N, ncol=Nmaxn)
        for(avar in 1:length(vnames$N)){
            Nalpha0[avar, 1:varinfoaux[name == vnames$N[avar], Nvalues]] <- 1
        }
    }

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
             Cleft = vtransform(dataset[,vnames$C, with=F], varinfoaux, Cout='left'),
             Cright = vtransform(dataset[,vnames$C, with=F], varinfoaux, Cout='right')
             ) },
    if(vn$D > 0){# discretized
        list(Dn = vn$D,
             Dmean1 = rep(0, 1),
             Dvarm1 = rep(1, 1),
             Dvar1 = rep(1, 1),
             Dshapelo = rep(Dshapelo, 1),
             Dshapehi = rep(Dshapehi, 1),
             Dleft = vtransform(dataset[,vnames$D, with=F], varinfoaux, Dout='left'),
             Dright = vtransform(dataset[,vnames$D, with=F], varinfoaux, Dout='right')
             ) },
    if(vn$O > 0){# ordinal
        list(On = vn$O,
             Omean1 = rep(0, 1),
             Ovarm1 = rep(1, 1),
             Ovar1 = rep(1, 1),
             Oshapelo = rep(Oshapelo, 1),
             Oshapehi = rep(Oshapehi, 1),
             Oleft = vtransform(dataset[,vnames$O, with=F], varinfoaux, Oout='left'),
             Oright = vtransform(dataset[,vnames$O, with=F], varinfoaux, Oout='right')
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


    datapoints <- c(
        if(vn$R > 0){# continuous
            list(
                Rdata = vtransform(dataset[,vnames$R, with=F], varinfoaux)
            ) },
        if(vn$C > 0){# censored
            list(
                Caux = vtransform(dataset[,vnames$C, with=F], varinfoaux, Cout='aux'),
                Clat = vtransform(dataset[,vnames$C, with=F], varinfoaux, Cout='lat')
            ) },
        if(vn$D > 0){# discretized
            list(
                Daux = vtransform(dataset[,vnames$D, with=F], varinfoaux, Dout='aux')
            ) },
        if(vn$O > 0){# ordinal
            list(
                Oaux = vtransform(dataset[,vnames$O, with=F], varinfoaux, Oout='aux')
            ) },
        if(vn$N > 0){# nominal
            list(
                Ndata = vtransform(dataset[,vnames$N,with=F], varinfoaux, Nout='numeric')                
            ) },
        if(vn$B > 0){# binary
            list(
                Bdata = vtransform(dataset[,vnames$B,with=F], varinfoaux, Bout='numeric')                
            ) }
    )

    cat('\nStarting Monte Carlo sampling with',nchains,'chains across', ncores, 'cores.\n')
    cat('Core logs are being saved in individual files.\n\nEst. Remaining Time...\r')
    ## stopCluster(cluster)
    stopImplicitCluster()
    registerDoSEQ()
    ## cl <- makePSOCKcluster(ncores)
    if(ncores > 1){
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
    }else{
        registerDoSEQ()
    }
    ## toexport <- c('constants', 'datapoints', 'vn', 'vnames', 'nalpha', 'nclusters')
    ## toexport <- c('vtransform','samplesFDistribution','proposeburnin','proposethinning','plotFsamples')

    mcsamples <- foreach(acore=1:ncores, .combine=rbind, .inorder=FALSE)%dorng%{

        outcon <- file(paste0(dirname,'_log-',basename,'-',acore,'.log'), open = "w")
        sink(outcon)
        sink(outcon, type = "message")
        suppressPackageStartupMessages(library('data.table'))
        suppressPackageStartupMessages(library('nimble'))

        source('pglpm_plotfunctions.R')
        source('vtransform.R')
        source('samplesFDistribution.R')
        source('proposeburnin.R')
        source('proposethinning.R')
        source('plotFsamples.R')

        thresholdfn <- function(diagnESS, diagnIAT, diagnBMK, diagnMCSE, diagnStat, diagnBurn, diagnBurn2, diagnThin){
            ceiling(2* max(diagnBurn2) + (nsamplesperchain-1L) * multcorr * ceiling(max(diagnIAT, diagnThin)))
    }


        
        ## predictors <- setdiff(unlist(vnames), predictands)

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



        initsfn <- function(){
            Alpha <- sample(1:nalpha, 1, prob=constants$probalpha0, replace=T)
            W <- 1/nclusters + 0*nimble::rdirch(n=1, alpha=constants$basealphas*2^Alpha)
            outlist <- list(
                Alpha = Alpha,
                W = W,
                K = rep(which(W>0)[1], npoints)
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
                                 Clat = vtransform(dataset[,vnames$C, with=F], varinfoaux, Cout='init') ## for data with boundary values
                             ))
            }
            if(vn$D > 0){# discretized
                Drate <- matrix(nimble::rinvgamma(n=vn$D*nclusters, shape=constants$Dshapehi, rate=constants$Dvar1), nrow=vn$D, ncol=nclusters)
                outlist <- c(outlist,
                             list(
                                 Dmean = matrix(rnorm(n=vn$D*nclusters, mean=constants$Dmean1, sd=sqrt(constants$Dvarm1)), nrow=vn$D, ncol=nclusters),
                                 Drate = Drate,
                                 Dvar = matrix(nimble::rinvgamma(n=vn$D*nclusters, shape=constants$Dshapelo, rate=Drate), nrow=vn$D, ncol=nclusters),
                                 Dlat = vtransform(dataset[,vnames$D, with=F], varinfoaux, Dout='init') ## for data with boundary values
                             ))
            }
            if(vn$O > 0){# ordinal
                Orate <- matrix(nimble::rinvgamma(n=vn$O*nclusters, shape=constants$Oshapehi, rate=constants$Ovar1), nrow=vn$O, ncol=nclusters)
                outlist <- c(outlist,
                             list(
                                 Omean = matrix(rnorm(n=vn$O*nclusters, mean=constants$Omean1, sd=sqrt(constants$Ovarm1)), nrow=vn$O, ncol=nclusters),
                                 Orate = Orate,
                                 Ovar = matrix(nimble::rinvgamma(n=vn$O*nclusters, shape=constants$Oshapelo, rate=Orate), nrow=vn$O, ncol=nclusters),
                                 Olat = vtransform(dataset[,vnames$O, with=F], varinfoaux, Oout='init') ## for data with boundary values
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
            monitors2=c( 'Alpha', 'K')
        )
        
        ## replace Alpha's cat-sampler and RW samplers with slice
        targetslist <- sapply(confnimble$getSamplers(), function(xx)xx$target)   
        nameslist <- sapply(confnimble$getSamplers(), function(xx)xx$name)
        for(asampler in c('Alpha', targetslist[nameslist == 'RW'])){
            confnimble$removeSamplers(asampler)
            confnimble$addSampler(target=asampler, type='slice')
        }
        ## call this to do a first reordering
        mcsampler <- buildMCMC(confnimble)

        ## change execution order for some variates
        samplerorder <- c('K', 
                          if(vn$R > 0){c('Rmean','Rrate','Rvar')},
                          if(vn$C > 0){c('Cmean','Crate','Cvar')},
                          if(vn$D > 0){c('Dmean','Drate','Dvar')},
                          if(vn$O > 0){c('Omean','Orate','Ovar')},
                          if(vn$N > 0){c('Nprob')},
                          if(vn$B > 0){c('Bprob')},
                          'W','Alpha')
        ##
        neworder <- foreach(var=samplerorder, .combine=c)%do%{grep(paste0('^',var,'(\\[.+\\])*$'), sapply(confnimble$getSamplers(), function(x)x$target))}
        ##
        confnimble$setSamplerExecutionOrder(c(setdiff(confnimble$getSamplerExecutionOrder(), neworder), neworder))
        print(confnimble)


        mcsampler <- buildMCMC(confnimble)
        Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

        
        cat('\nSetup', capture.output(print(Sys.time() - timecount)), '\n')

##################################################
        ## Monte Carlo sampler and plots of MC diagnostics
##################################################
        time0 <- Sys.time()
        nitertot <- 0L
        achain <- 0L
        continue <- TRUE
        newchain <- TRUE
        allmcsamples <- NULL
        maxusedclusters <- 0
        testdata <- t(varinfoaux[,paste0('mctest',1:3)])
        colnames(testdata) <- varinfoaux[,name]
        ## testdata <- rbind(varinfoaux[['Q2']], varinfo[['Q1']], varinfo[['Q3']]) #, varinfo[['plotmin']], varinfo[['plotmax']], varinfo[['datamin']], varinfo[['datamax']])

        calctime <- Sys.time()
        while(continue){
            if(newchain){
                niter <- nitertot <- niter0
                reset <- TRUE
                traces <- mcsamples <- prevmcsamples <- NULL
                achain <- achain + 1L
                ## if(TRUE){ # send message to user screen with est. remaining time
                    sink(NULL,type='message')
                    message(paste0('\rEst. Remaining ',
                                   capture.output(print((Sys.time()-calctime)/(achain-1)*(nchainspercore-achain+1)))), appendLF=FALSE)
                flush.console()
                    sink(outcon,type='message')
                ##}
                if(!(achain > nchainspercore)){
                    mcmcseed <- (acore-1L)*nchainspercore + achain
                    cat('Seed:', mcmcseed+seed, '\n')
                    set.seed(mcmcseed+seed)
                    Cfinitemixnimble$setInits(initsfn())
                }
            }else{
                prevmcsamples <- rbind(prevmcsamples,mcsamples)
                niter <- lengthmeasure - nitertot + 1L
                nitertot <- lengthmeasure
                reset <- FALSE
            }
            showsamplertimes0 <- showsamplertimes && (achain==1) && newchain
            showhyperparametertraces0 <- showhyperparametertraces && (achain==1) && newchain

            if(achain > nchainspercore){
                continue <- FALSE
                achain <- 'F'
                mcsamples <- allmcsamples
                usedclusters <- maxusedclusters
                ## saveRDS(mcsamples, file=paste0(dirname,'_mcsamples-R',basename,'--',mcmcseed,'-',achain,'.rds'))
            }else{
                ##
                cat('Iterations:', niter,'\n')
                cat('chain:', achain,'of',nchainspercore,'. Est. Remaining',
                capture.output(print((Sys.time()-calctime)/(achain-1)*(nchainspercore-achain+1))), '\n')
                Cmcsampler$run(niter=niter, thin=1, thin2=niter, nburnin=0, time=showsamplertimes0, reset=reset, resetMV=TRUE)
                mcsamples <- as.matrix(Cmcsampler$mvSamples)
                finalstate <- as.matrix(Cmcsampler$mvSamples2)
                cat('\nMCMC', capture.output(print(Sys.time() - calctime)), '\n')

                if(any(!is.finite(mcsamples))){cat('\nWARNING: SOME NON-FINITE OUTPUTS')}
                
                if(showsamplertimes){
                    samplertimes <- Cmcsampler$getTimes()
                    names(samplertimes) <- sapply(confnimble$getSamplers(),function(x)x$target)
                    sprefixes <- unique(sub('^([^[]+)(\\[.*\\])', '\\1', names(samplertimes)))
                    cat('\nSampler times:\n')
                    print(sort(sapply(sprefixes, function(x)sum(samplertimes[grepl(x,names(samplertimes))])),decreasing=T))
                }
                ##
                if(showhyperparametertraces){
                    occupations <- apply(finalstate[,grepl('^K\\[', colnames(finalstate)),drop=F], 1, function(xxx){length(unique(xxx))})
                    cat('\nSTATS OCCUPIED CLUSTERS:\n')
                    print(summary(occupations))
                    ##
                    pdff(paste0(dirname,'traces_hyperparameters-',mcmcseed,'-',achain), apaper=4)
                    tplot(y=occupations, ylab='occupied clusters',xlab=NA,ylim=c(0,nclusters))
                    histo <- thist(occupations,n='i')
                    tplot(x=histo$mids,y=histo$density,xlab='occupied clusters',ylab=NA)
                    tplot(y=log2(alpha0[finalstate[,'Alpha']]), ylab='Alpha index',xlab=NA)
                    histo <- thist(log2(alpha0[finalstate[,'Alpha']]))
                    tplot(x=histo$mids,y=histo$density,xlab='Alpha index',ylab='')
                    for(vtype in c('R','C','D','O','N','B')){
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
                cat('\nOCCUPIED CLUSTERS:', usedclusters, 'OF', nclusters,'\n')
            }
            ##
            ## Diagnostics
            ## Log-likelihood
            diagntime <- Sys.time()
            ## ll <- colMeans(log(samplesFDistribution(Y=data.matrix(data0), X=NULL, mcsamples=mcsamples, varinfo=varinfo, jacobian=FALSE)),na.rm=T) #- sum(log(invjacobian(data.matrix(data0), varinfo)), na.rm=T)
            ## lld <- colMeans(log(samplesFDistribution(Y=data.matrix(data0[,..predictands]), X=data.matrix(data0[,..predictors]), mcsamples=mcsamples, varinfo=varinfo, jacobian=FALSE)),na.rm=T) # - sum(log(invjacobian(data.matrix(data0[,..predictands]), varinfo)), na.rm=T)
            ## lli <- colMeans(log(samplesFDistribution(Y=data.matrix(data0[,..predictors]), X=data.matrix(data0[,..predictands]), mcsamples=mcsamples, varinfo=varinfo, jacobian=FALSE)),na.rm=T) #- sum(log(invjacobian(data.matrix(data0[,..predictors]), varinfo)), na.rm=T)
            ##
            ll <- t(
                log(samplesFDistribution(Y=testdata, X=NULL, mcsamples=mcsamples, varinfoaux=varinfoaux, subsamples=NULL, jacobian=FALSE, parallel=FALSE)) #- sum(log(invjacobian(data.matrix(data0), varinfo)), na.rm=T)
            )
            colnames(ll) <- paste0('log-',c('mid','lo','hi')) #,'pm','pM','dm','dM'))
            ## testdatalld <- log(samplesFDistribution(Y=testdata[,predictands,drop=F], X=testdata[,predictors,drop=F], mcsamples=mcsamples, varinfo=varinfo, jacobian=FALSE)) # - sum(log(invjacobian(data.matrix(data0[,..predictands]), varinfo)), na.rm=T)
            ## rownames(testdatalld) <- paste0(c('Q2','Q1','Q3'),'d')
            ## testdatalli <- log(samplesFDistribution(Y=testdata[,predictors,drop=F], X=testdata[,predictands,drop=F], mcsamples=mcsamples, varinfo=varinfo, jacobian=FALSE)) #- sum(log(invjacobian(data.matrix(data0[,..predictors]), varinfo)), na.rm=T)
            ## rownames(testdatalli) <- paste0(c('Q2','Q1','Q3'),'i')
            ## stestdatall <- colSums(testdatall, na.rm=T)
            ## stestdatalld <- colSums(testdatalld, na.rm=T)
            ## stestdatalli <- colSums(testdatalli, na.rm=T)
            ##
            traces <- rbind(traces, 10/log(10) * ll)
            ## pdff('testdifferenttraces2')
            ## tplot(y=traces, lwd=1, lty=1)
            ## for(i in 1:ncol(traces)){
            ##     tplot(y=traces[,i], main=colnames(traces)[i])
            ## }
            ## dev.off()
            traces2 <- traces[apply(traces,1,function(x){all(is.finite(x))}),]
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
            cat('\nESSs:',paste0(round(diagnESS), collapse=', '))
            cat('\nIATs:',paste0(round(diagnIAT), collapse=', '))
            cat('\nBMKs:',paste0(round(diagnBMK,3), collapse=', '))
            cat('\nMCSEs:',paste0(round(diagnMCSE,2), collapse=', '))
            cat('\nStationary:',paste0(diagnStat, collapse=', '))
            cat('\nBurn-in I:',paste0(diagnBurn, collapse=', '))
            cat('\nBurn-in II:',diagnBurn2)
            cat('\nProposed thinning:',paste0(diagnThin, collapse=', '),'\n')

            cat('\nDiagnostics', capture.output(print(Sys.time() - diagntime)), '\n')

#########################################
#### CHECK IF WE NEED TO SAMPLE MORE ####
#########################################
            if(continue){
                lengthmeasure <- thresholdfn(diagnESS=diagnESS, diagnIAT=diagnIAT, diagnBMK=diagnBMK, diagnMCSE=diagnMCSE, diagnStat=diagnStat, diagnBurn=diagnBurn, diagnBurn2=diagnBurn2, diagnThin=diagnThin)
                cat('\nNumber of iterations', nitertot, ', required', lengthmeasure,'\n')
                ##
                if(nitertot < lengthmeasure){
                    cat('Increasing by', lengthmeasure-nitertot+1L, '\n')
                    newchain <- FALSE
                }else{
                    mcsamples <- rbind(prevmcsamples, mcsamples)
                    allmcsamples <- rbind(allmcsamples, mcsamples[rev( nrow(mcsamples) - seq(from=0, length.out=nsamplesperchain, by=multcorr*ceiling(max(diagnIAT,diagnThin))) ),,drop=F])
                    maxusedclusters <- max(usedclusters, maxusedclusters)
                    newchain <- TRUE
                }
            }

#########################################
#### END CHECK                       ####
#########################################
            ##
            if(newchain && saveallchains){
                saveRDS(traces,file=paste0(dirname,'_mctraces-',basename,'--',mcmcseed,'-',achain,'.rds'))
                saveRDS(allmcsamples, file=paste0(dirname,'_mcsamples-',basename,'--',mcmcseed,'-','F','.rds'))
            }

            if(newchain && (plotallchains || !continue)){
                ##
###############
#### PLOTS ####
###############
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
                pdff(paste0(dirname,'_mcmcpartialtraces-',basename,'--',mcmcseed,'-',achain), apaper=4)
                ## Summary stats
                matplot(1:2, type='l', col='white', main=paste0('Stats chain ',achain), axes=FALSE, ann=FALSE)
                legendpositions <- c('topleft','topright','bottomleft','bottomright')
                for(alegend in 1:length(grouplegends)){
                    legend(x=legendpositions[alegend], bty='n', cex=1.5,
                           legend=grouplegends[[alegend]] )
                }
                legend(x='center', bty='n', cex=1,
                       legend=c(
                           paste0('Chain ', achain),
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
                                      ' | stat: ', diagnStat[avar]*1L,
                                      ' | burnI: ', diagnBurn[avar],
                                      ' | burnII: ', diagnBurn2,
                                      ' | thin: ', diagnThin[avar]
                                      ), 
                          ylab=paste0(avar,'/dHart'), xlab='sample', family=family
                          )
                }
                dev.off()
                ##
                ## Samples of marginal frequency distributions
                if(!continue){
                    subsamples <- (if(totsamples=='all'){1:nrow(mcsamples)}else{round(seq(1, nrow(mcsamples), length.out=totsamples))})
                    ## showsubsample <- round(seq(1, length(subsamples), length.out=showsamples))
                    ##
                    cat('\nPlotting samples of frequency distributions')

                    plotFsamples(file=paste0(dirname,'mcmcdistributions-',basename,'--',mcmcseed,'-',achain),
                                 mcsamples=mcsamples[subsamples,,drop=F],
                                 varinfoaux=varinfoaux,
                                 dataset=dataset,
                                 nsubsamples=showsamples,
                                 plotmeans=plotmeans, showdata='histogram',
                                 parallel=FALSE
                                 )
                }
                ##
                cat('\nMCMC+diagnostics', capture.output(print(Sys.time() - calctime)), '\n')
                ##

            }
        }
        ##
        cat('\nTotal', capture.output(print(Sys.time() - time0)), '\n')

        ## Final output of foreach
        cbind(traces,mcsamples)
    }
############################################################
        ## End MCMC
############################################################

    attr(mcsamples, 'rng') <- NULL
    attr(mcsamples, 'doRNG_version') <- NULL
    traces <- mcsamples[round(seq(1,nrow(mcsamples),length.out=nsamples)),1:3]
    mcsamples <- mcsamples[round(seq(1,nrow(mcsamples),length.out=nsamples)),-(1:3)]
    saveRDS(mcsamples,file=paste0(dirname,'Fdistribution-',basename,'-',nsamples,'.rds'))
    saveRDS(traces,file=paste0(dirname,'MCtraces-',basename,'-',nsamples,'.rds'))
    cat('\nFinished Monte Carlo sampling.\n')
    gc()

############################################################
        ## Final joint diagnostics
############################################################
cat('\nSome diagnostics:\n')
    traces2 <- traces[apply(traces,1,function(x){all(is.finite(x))}),]
    flagll <- nrow(traces) != nrow(traces2)


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
    cat('\nESSs:',paste0(round(diagnESS), collapse=', '))
    cat('\nIATs:',paste0(round(diagnIAT), collapse=', '))
    cat('\nBMKs:',paste0(round(diagnBMK,3), collapse=', '))
    cat('\nMCSEs:',paste0(round(diagnMCSE,2), collapse=', '))
    cat('\nStationary:',paste0(diagnStat, collapse=', '))
    cat('\nBurn-in I:',paste0(diagnBurn, collapse=', '))
    cat('\nBurn-in II:',diagnBurn2)
    cat('\n')
    ## cat(paste0('\nProposed thinning: ',paste0(diagnThin, collapse=', ')))
    ##       )
    ##         }
    colpalette <- c(7,2,1)
    names(colpalette) <- colnames(traces)
##

    ## Plot various info and traces
    cat('\nPlotting final Monte Carlo traces and marginal samples.\n')
    
    ##
    graphics.off()
    pdff(paste0(dirname,'MCtraces-',basename,'-',nsamples), apaper=4)
    ## Traces of likelihood and cond. probabilities
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
    
    plotFsamples(file=paste0(dirname,'Fdistribution-',basename,'-',nsamples),
                 mcsamples=mcsamples, varinfoaux=varinfoaux,
                 dataset=dataset,
                 nsubsamples=showsamples, plotmeans=TRUE,
                 showdata = 'histogram', parallel=TRUE)
    cat('\nClosing connections to cores.\n')
    registerDoSEQ()
    stopCluster(cl)

        cat('\nComputation', capture.output(print(Sys.time() - timestart0)), '\n')
    cat('Finished.\n\n')
    
    mcsamples
}
