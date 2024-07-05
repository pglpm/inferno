mcinputs <- readRDS(commandArgs(trailingOnly=TRUE))


function(dataset, varinfoaux, predictands, outputdir=TRUE, nsamples=4096, nsamplesperchain=4, nchains=4, ncores=NULL, seed=701){

    if(!missing(nsamples) && !missing(nchains) && !missing(nsamplesperchain)){
stop('specify only one of "nsamples", "nchains", "nsamplesperchain"')
    }
 
        require('data.table')
        require('LaplacesDemon', include.only=NULL)

    ## Read MCMC seed from command line
arguments <- as.integer(commandArgs(trailingOnly=TRUE))[1]
mcmcseed <- arguments
if(is.na(mcmcseed) || (!is.na(mcmcseed) && mcmcseed <= 0)){mcmcseed <- 1}
cat(paste0('\nMCMC seed = ',mcmcseed,'\n'))
##

    ## read dataset
    if(is.character(dataset) && file.exists(dataset)){dataset <- fread(dataset, na.strings='')}
    dataset <- as.data.table(dataset)

    ## varinfoaux
    if(is.character(varinfoaux) && file.exists(varinfoaux)){
        varinfoaux <- readRDS(varinfoaux)
    }

    ## list of predictand variates
    if(is.character(predictands) && file.exists(predictands)){
        predictands <- as.vector(unlist(read.csv(predictands, header=F)))
    }else{
        predictands <- unlist(predictands)
    }
    
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

    ## set number of cores for parallel computation
    if(is.null(cores)){
        warning('The number of cores has not been given.\nIt is much preferable that it be set by the user.')
        cores <- round(detectCores()/2)
        cat('\nTrying to use ',cores,' cores\n')
    }
    
#### various internal parameters
    niter0 <- 1024L # 3L # iterations to try
#### Hyperparameters
    nclusters <- 4L
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
    plottempdistributions <- FALSE # plot temporary sampled distributions
    showdata <- TRUE # 'histogram' 'scatter' FALSE TRUE
    plotmeans <- TRUE # plot frequency averages
    totsamples <- 'all' # 'all' number of samples if plotting frequency averages
    showsamples <- 100 # number of samples to show. Shown separately for posterior=F
    testLength <- TRUE
    nthreshold <- 2 # multiple of threshold for acceptable number of burn-in samples
    showhyperparametertraces <- FALSE ##
    showsamplertimes <- FALSE ##
    family <- 'Palatino'

    basename <- paste0(outputdir,'-V',nrow(varinfoaux),'-D',npoints,'-K',nclusters,'-I',nsamples)
##
    dirname <- paste0(basename,'/')
    dir.create(dirname)



    vn <- list()
    vnames <- list()
    for(atype in c('R','C','D','O','N','B')){
        vn[[atype]] <- length(varinfoaux[mcmctype == atype, name])
        vnames[[atype]] <- varinfoaux[mcmctype == atype, name]
    }
    ## choose predictands if unchosen
    if(is.null(predictands) || !exists('predictands')){
        predictands <- sapply(vnames,function(xx)xx[1])
        cat('\nSelf-choosing predictand variates:\n', predictands, '\n')
    }
    predictors <- setdiff(unlist(vnames), predictands)

    if(vn$N > 0){
        Nmaxn <- max(varinfoaux[mcmctype == 'N', Nvalues])
        Nalpha0 <- matrix(1e-100, nrow=vn$N, ncol=Nmaxn)
        for(avar in 1:length(vnames$N)){
            Nalpha0[avar, 1:varinfoaux[name == vnames$N[avar], Nvalues]] <- 1
        }
    }

stopCluster(cluster)
stopImplicitCluster()
registerDoSEQ()
## cl <- makePSOCKcluster(ncores)
cluster <- makeCluster(ncores, outfile='')
registerDoParallel(cluster)

    
    mcsamples <- foreach(chain=1:ncores, .combine=rbind, .packages=c('nimble','data.table'), .inorder=FALSE)%dorng%{

        ## print(predictands)
        ## print(predictors)
        ## print(dataset[,..predictands])
        ## print(dataset[,..predictors])

        source('vtransform.R')
        source('samplesFDistribution.R')
        source('proposeburnin.R')
        source('proposethinning.R')

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
            monitors2=c( 'Alpha', 'K','Clat') # ****remove Clat****
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


        mcsampler <- buildMCMC(confnimble)
        Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

        
        cat('\nSetup time: ')
        print(Sys.time() - timecount)

## #### ****remove tests below****
##         Cmcsampler$run(niter=1024, thin=1, thin2=64, nburnin=0, time=FALSE, reset=TRUE, resetMV=TRUE)

##         mcsamples <- as.matrix(Cmcsampler$mvSamples)
##         finalstate <- as.matrix(Cmcsampler$mvSamples2)

##         tplot(y=mcsamples2[,'Alpha'])

##         tplot(y=mcsamples[,'Rmean[1, 1]'])

##         tplot(y=mcsamples2[,'Clat[3, 1]'])

##         tplot(y=log(apply(mcsamples,1,function(xx){min(xx[paste0('W[',1:nclusters,']')])})))


##################################################
#### Monte Carlo sampler and plots of MC diagnostics
##################################################
        time0 <- Sys.time()
        traces <- NULL
        niterfinal <- niter0
        nitertot <- 0L
        stage <- 0L
        continue <- TRUE

        set.seed(chain+1000)
        Cfinitemixnimble$setInits(initsfn())
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
                set.seed(chain+1000)
            }

            cat(paste0('Iterations: ', niter),'\n')
            mcsamples <- finalstate <- NULL
            calctime <- Sys.time()
            for(achain in 1:nsamples0){
                cat(paste0('chain: ', achain,', est. remaining time: ')); print((Sys.time()-calctime)/(achain-1)*(nsamples0-achain+1))
                if(!continue){
                    ## set.seed(mcmcseed+achain+999)
                    Cfinitemixnimble$setInits(initsfn())
                }
                Cmcsampler$run(niter=niter, thin=1, thin2=1, nburnin=nburnin, time=showsamplertimes0, reset=reset, resetMV=TRUE)
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
                pdff(paste0(dirname,'traces_hyperparameters-',chain,'-',stage))
                tplot(y=occupations, ylab='occupied clusters',xlab=NA,ylim=c(0,nclusters))
                histo <- thist(occupations,n='i')
                tplot(x=histo$mids,y=histo$density,xlab='occupied clusters',ylab=NA)
                tplot(y=finalstate[,'Alpha'], ylab='Alpha-index',xlab=NA)
                histo <- thist(finalstate[,'Alpha'])
                tplot(x=histo$mids,y=histo$density,xlab='Alpha-index',ylab='')
                ## for(vtype in c('R','I',#'O',
                ##                'D')){
                ##     if(len[[vtype]] > 0){
                ##         for(togrep in c('varscaleindex')){
                ##             for(v in colnames(finalstate)[grepl(paste0('^',vtype,togrep,'\\['), colnames(finalstate))]){
                ##                 tplot(y=finalstate[,v],ylab=v,xlab=NA,ylim=c(1,(2*hwidth+1)))
                ##             }
                ##         }
                ##     }
                ## }
                dev.off()
            }

            finalstate <- c(mcsamples[nrow(mcsamples),], finalstate[nrow(finalstate),])
            ## Check how many "clusters" were occupied. Warns if too many
            occupations <- finalstate[grepl('^K\\[', names(finalstate))]
            usedclusters <- length(unique(occupations))
            if(usedclusters > nclusters-5){cat('\nWARNING: TOO MANY CLUSTERS OCCUPIED')}
            cat(paste0('\nOCCUPIED CLUSTERS: ', usedclusters, ' OF ', nclusters),'\n')


            ## Diagnostics
            ## Log-likelihood
            diagntime <- Sys.time()
            ll <- colSums(log(samplesFDistribution(Y=dataset, X=NULL, mcsamples=mcsamples, varinfoaux=varinfoaux, jacobian=FALSE)), na.rm=T) #- sum(log(invjacobian(data.matrix(dataset), varinfo)), na.rm=T)
            lld <- colSums(log(samplesFDistribution(Y=dataset[,..predictands], X=dataset[,..predictors], mcsamples=mcsamples, varinfoaux=varinfoaux, jacobian=FALSE)), na.rm=T) # - sum(log(invjacobian(data.matrix(dataset[,..predictands]), varinfo)), na.rm=T)
            lli <- colSums(log(samplesFDistribution(Y=dataset[,..predictors], X=dataset[,..predictands], mcsamples=mcsamples, varinfoaux=varinfoaux, jacobian=FALSE)), na.rm=T) #- sum(log(invjacobian(data.matrix(dataset[,..predictors]), varinfo)), na.rm=T)
            ##
            traces <- rbind(traces,
                            10/log(10)/npoints *
                            cbind(loglikelihood=ll,
                                  'mean of direct logprobabilities'=lld,
                                  'mean of inverse logprobabilities'=lli)
                            )
            traces2 <- traces[apply(traces,1,function(x){all(is.finite(x))}),]
            saveRDS(traces,file=paste0(dirname,'_mctraces-R',basename,'--',chain,'-',stage,'.rds'))
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
                saveRDS(mcsamples, file=paste0(dirname,'_mcsamples-R',basename,'--',chain,'-F.rds'))
                saveRDS(traces, file=paste0(dirname,'_mctraces-R',basename,'--',chain,'-F.rds'))
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
        pdff(paste0(dirname,'mcmcplottraces-R',basename,'--',chain,'-',stage),'a4')
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
        pdff(paste0(dirname,'mcmcdistributions-R',basename,'--',chain,'-',stage),'a4')
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
                    datum <- dataset[[v]]
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
                    datum <- dataset[[v]]
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
                    datum <- dataset[[v]]
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
mcsamples
        }
##
cat('\nTotal time: ')
print(Sys.time() - time0)

############################################################
## End MCMC
############################################################
registerDoSEQ()

stop('NONE. End of script')
