inferpopulation <- function(data, varinfoaux, predictands, nsamples=4096, file=TRUE, ncores=NULL){
    ## read data
    if(is.character(data) && file.exists(data)){data <- fread(data, na.strings='')}
    data <- as.data.table(data)
    ## varinfoaux
    if(is.character(varinfoaux) && file.exists(varinfoaux)){
        varinfoauxname <- varinfoaux
        varinfoaux <- readRDS(varinfoaux)
    }
    ## list of predictand variates
    if(is.character(predictands) && file.exists(predictands)){
        predictandsname <- predictands
        predictands <- as.vector(unlist(read.csv(predictands, header=F)))
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
    ##
#### various internal parameters
    niter0 <- 1024L # 3L # iterations to try
    testLength <- TRUE
    nthreshold <- 2 # multiple of threshold for acceptable number of burn-in samples
    casualinitvalues <- FALSE
    showhyperparametertraces <- FALSE ##
    showsamplertimes <- FALSE ##
    family <- 'Palatino'
#### Hyperparameters
    nclusters <- 64L
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

    vn <- list()
    vnames <- list()
    for(atype in c('R','C','D','O','N','B')){
        vn[[atype]] <- length(varinfoaux[mcmctype==atype, name])
        vnames[[atype]] <- varinfoaux[mcmctype==atype, name]
    }

    if(vn$N > 0){
        Nmaxn <- max(varinfoaux[mcmctype=='N', Nvalues])
        Nalpha0 <- matrix(1e-100, nrow=vn$N, ncol=Nmaxn)
        for(avar in 1:length(vnames$N)){
            Nalpha0[avar, 1:varinfoaux[name==vnames$N[avar], Nvalues]] <- 1
        }
    }
    
    mcsamples <- foreach(chain=1:ncores, .combine=rbind, .packages='nimble', .inorder=FALSE)%dorng%{
        ##
        ## hierarchical probability structure
        finitemix <- nimble::nimbleCode({
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
                        Rmean[v, k] ~ dnorm(mean=Rmean1[v], var=Rvarm1[v])
                        Rrate[v, k] ~ dinvgamma(shape=Rshapehi[v], rate=Rvar1[v])
                        Rvar[v, k] ~ dinvgamma(shape=Rshapelo[v], rate=Rrate[v, k])
                    }
                }
                if(vn$C > 0){# censored
                    for(v in 1:Cn){
                        Cmean[v, k] ~ dnorm(mean=Cmean1[v], var=Cvarm1[v])
                        Crate[v, k] ~ dinvgamma(shape=Cshapehi[v], rate=Cvar1[v])
                        Cvar[v, k] ~ dinvgamma(shape=Cshapelo[v], rate=Crate[v, k])
                    }
                }
                if(vn$D > 0){# discretized
                    for(v in 1:Dn){
                        Dmean[v, k] ~ dnorm(mean=Dmean1[v], var=Dvarm1[v])
                        Drate[v, k] ~ dinvgamma(shape=Dshapehi[v], rate=Dvar1[v])
                        Dvar[v, k] ~ dinvgamma(shape=Dshapelo[v], rate=Drate[v, k])
                    }
                }
                if(vn$O > 0){# ordinal
                    for(v in 1:On){
                        Omean[v, k] ~ dnorm(mean=Omean1[v], var=Ovarm1[v])
                        Orate[v, k] ~ dinvgamma(shape=Oshapehi[v], rate=Ovar1[v])
                        Ovar[v, k] ~ dinvgamma(shape=Oshapelo[v], rate=Orate[v, k])
                    }
                }
                if(vn$N > 0){# nominal
                    for(v in 1:Nn){
                        Nprob[v, k, 1:Nmaxn] ~ ddirch(alpha=Nalpha0[v, 1:Nmaxn])
                    }
                }
                if(vn$B > 0){# binary
                    for(v in 1:Bn){
                        Bprob[v, k] ~ dbeta(shape1=Bshapelo[v], shape2=Bshapehi[v])
                    }
                }
            }
            ## Probability of data
            for(d in 1:npoints){
                K[d] ~ dcat(prob=W[1:nclusters])
                ##
                if(vn$R > 0){# continuous
                    for(v in 1:Rn){
                        Rdata[d, k] ~ dnorm(mean=Rmean[v, K[d]], var=Rvar[v, K[d]])
                    }
                }
                if(vn$C > 0){# censored
                    for(v in 1:Cn){
                        Caux[v, d] ~ dconstraint(Clat[v, d] >= Cleft[v] & Clat[v, d] <= Cright[v])
                        Clat[v, d] ~ dnorm(mean=Cmean[v, K[d]], var=Cvar[v, K[d]])
                    }
                }
                if(vn$D > 0){# discretized
                    for(v in 1:Dn){
                        Daux[v, d] ~ dconstraint(Dlat[v, d] >= Dleft[v, d] & Dlat[v, d] < Dright[v, d])
                        Dlat[v, d] ~ dnorm(mean=Dmean[v, K[d]], var=Dvar[v, K[d]])
                    }
                }
                if(vn$O > 0){# ordinal
                    for(v in 1:On){
                        Oaux[v, d] ~ dconstraint(Olat[v, d] >= Oleft[v, d] & Olat[v, d] < Oright[v, d])
                        Olat[v, d] ~ dnorm(mean=Omean[v, K[d]], var=Ovar[v, K[d]])
                    }
                }
                if(vn$N > 0){# nominal
                    for(v in 1:Nn){
                        Ndata[v, d] ~ dcat(prob=Nprob[v, K[d], 1:Nmaxn])
                    }
                }
                if(vn$B > 0){# binary
                    for(v in 1:Bn){
                        Bdata[v, d] ~ dbern(prob=Bprob[v, K[d]])
                    }
                }
            }
        })
        ##
        ##
        constants <- c( list(
            nclusters = nclusters,
            npoints = npoints,
            nalpha = nalpha,
            probalpha0 = rep(1/nalpha, nalpha),
            basealphas = rep((2^(minalpha-1L))/nclusters, nclusters)
        ),
        if(vn$R > 0){# continuous
            list(Rn = vn$R,
                 Rmean1 = rep(0, vn$R),
                 Rvarm1 = rep(1, vn$R),
                 Rvar1 = rep(1, vn$R),
                 Rshapelo = rep(Rshapelo, vn$R),
                 Rshapehi = rep(Rshapehi, vn$R)
                 ) },
        if(vn$C > 0){# censored
            list(Cn = vn$C,
                 Cmean1 = rep(0, vn$C),
                 Cvarm1 = rep(1, vn$C),
                 Cvar1 = rep(1, vn$C),
                 Cshapelo = rep(Cshapelo, vn$C),
                 Cshapehi = rep(Cshapehi, vn$C),
                 Cleft = vtransform(data0[,vnames$C, with=F], varinfoaux, Cout='left'),
                 Cright = vtransform(data0[,vnames$C, with=F], varinfoaux, Cout='right')
                 ) },
        if(vn$D > 0){# discretized
            list(Dn = vn$D,
                 Dmean1 = rep(0, vn$D),
                 Dvarm1 = rep(1, vn$D),
                 Dvar1 = rep(1, vn$D),
                 Dshapelo = rep(Dshapelo, vn$D),
                 Dshapehi = rep(Dshapehi, vn$D),
                 Dleft = vtransform(data0[,vnames$D, with=F], varinfoaux, Dout='left'),
                 Dright = vtransform(data0[,vnames$D, with=F], varinfoaux, Dout='right')
                 ) },
        if(vn$O > 0){# ordinal
            list(On = vn$O,
                 Omean1 = rep(0, vn$O),
                 Ovarm1 = rep(1, vn$O),
                 Ovar1 = rep(1, vn$O),
                 Oshapelo = rep(Oshapelo, vn$O),
                 Oshapehi = rep(Oshapehi, vn$O),
                 Oleft = vtransform(data0[,vnames$O, with=F], varinfoaux, Oout='left'),
                 Oright = vtransform(data0[,vnames$O, with=F], varinfoaux, Oout='right')
                 ) },
        if(vn$N > 0){# nominal
            list(Nn = vn$N,
                 Nmaxn = Nmaxn,
                 Nalpha0 = Nalpha0
                 ) },
        if(vn$B > 0){# binary
            list(Bn = vn$B,
                 Bshapelo = rep(Bshapelo, vn$B),
                 Bshapehi = rep(Bshapehi, vn$B)
                 ) },
        )
        ##
        ##
        datapoints <- c(
            if(vn$R > 0){# continuous
                list(
                    Rdata = vtransform(data0[,vnames$R, with=F], varinfoaux)
                ) },
            if(vn$C > 0){# censored
                list(
                    Caux = vtransform(data0[,vnames$C, with=F], varinfoaux, Cout='aux'),
                    Clat = vtransform(data0[,vnames$C, with=F], varinfoaux, Cout='lat')
                ) },
            if(vn$D > 0){# discretized
                list(
                    Daux = vtransform(data0[,vnames$D, with=F], varinfoaux, Dout='aux')
                ) },
            if(vn$O > 0){# ordinal
                list(
                    Oaux = vtransform(data0[,vnames$O, with=F], varinfoaux, Oout='aux')
                ) },
            if(vn$N > 0){# nominal
                list(
                    Ndata = transf(data0[,variate$N,with=F], varinfo)                
                ) },
            if(vn$B > 0){# binary
                list(
                    Bdata = transf(data0[,variate$B,with=F], varinfo)                
                ) }
        )
        ##
        ##
        initsfn <- function(){
            Alpha <- sample(1:nalpha, 1, prob=constants$probalpha0, replace=T),
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
                                 Clat = vtransform(data0[,vnames$C, with=F], varinfoaux, Cout='init') ## for data with boundary values
                             ))
            }
            if(vn$D > 0){# discretized
                Drate <- matrix(nimble::rinvgamma(n=vn$D*nclusters, shape=constants$Dshapehi, rate=constants$Dvar1), nrow=vn$D, ncol=nclusters)
                outlist <- c(outlist,
                             list(
                                 Dmean = matrix(rnorm(n=vn$D*nclusters, mean=constants$Dmean1, sd=sqrt(constants$Dvarm1)), nrow=vn$D, ncol=nclusters),
                                 Drate = Drate,
                                 Dvar = matrix(nimble::rinvgamma(n=vn$D*nclusters, shape=constants$Dshapelo, rate=Drate), nrow=vn$D, ncol=nclusters),
                                 Dlat = vtransform(data0[,vnames$D, with=F], varinfoaux, Dout='init') ## for data with boundary values
                             ))
            }
            if(vn$O > 0){# ordinal
                Orate <- matrix(nimble::rinvgamma(n=vn$O*nclusters, shape=constants$Oshapehi, rate=constants$Ovar1), nrow=vn$O, ncol=nclusters)
                outlist <- c(outlist,
                             list(
                                 Omean = matrix(rnorm(n=vn$O*nclusters, mean=constants$Omean1, sd=sqrt(constants$Ovarm1)), nrow=vn$O, ncol=nclusters),
                                 Orate = Orate,
                                 Ovar = matrix(nimble::rinvgamma(n=vn$O*nclusters, shape=constants$Oshapelo, rate=Orate), nrow=vn$O, ncol=nclusters),
                                 Olat = vtransform(data0[,vnames$O, with=F], varinfoaux, Oout='init') ## for data with boundary values
                             ))
            }
            if(vn$N > 0){# nominal
                outlist <- c(outlist,
                             list(
                                 Nprob = aperm(array(sapply(1:vn$N, function(xx){sapply(1:nclusters, function(avar){rdirch(n=1, alpha=Nalpha0[avar,])})}), dim=c(Nmaxn,nclusters,vn$N)))
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
                       if(len$R > 0){c('Rmean', 'Rvar')},
                       if(len$C > 0){c('Cmean', 'Cvar')},
                       if(len$D > 0){c('Dmean', 'Dvar')},
                       if(len$O > 0){c('Omean', 'Ovar')},
                       if(len$N > 0){c('Nprob')},
                       if(len$B > 0){c('Bprob')},
                       ),
            monitors2=c( 'Alpha', 'K')
        )
        ## replace Alpha's cat-sampler and RW samplers with slice
        targetslist <- sapply(confnimble$getSamplers(), function(xx)xx$target)   
        rwlist <- targetslist[which(sapply(confnimble$getSamplers(), function(xx)xx$name) == 'RW')]
        for(asampler in c('Alpha', rwlist)){
            confnimble$removeSamplers(asampler)
            confnimble$addSampler(target=asampler, type='slice')
        }
        ## reorder samplers
        samplerorder <- c('K',
                          if(vn$R > 0){c('Rmean','Rrate','Rvar')},
                          if(vn$C > 0){c('Cmean','Crate','Cvar')},
                          if(vn$D > 0){c('Dmean','Drate','Dvar')},
                          if(vn$O > 0){c('Omean','Orate','Ovar')},
                          if(vn$N > 0){c('Nprob')},
                          if(vn$B > 0){c('Bprob')},
                          'W','Alpha')
        neworder <- foreach(var=samplerorder, .combine=c)%do%{grep(paste0('^',var,'(\\[.+\\])*$'),targetslist)}
        neworder <- c(neworder, setdiff(confnimble$getSamplerExecutionOrder(), neworder))
        confnimble$setSamplerExecutionOrder(neworder)

        mcsampler <- buildMCMC(confnimble)
        Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)
        
        cat('\nSetup time: ')
        print(Sys.time() - timecount)

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
                set.seed(mcmcseed+1000)
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
                pdff(paste0(dirname,'traces_hyperparameters-',mcmcseed,'-',stage))
                tplot(y=occupations, ylab='occupied clusters',xlab=NA,ylim=c(0,nclusters))
                histo <- thist(occupations,n='i')
                tplot(x=histo$mids,y=histo$density,xlab='occupied clusters',ylab=NA)
                tplot(y=finalstate[,'Alpha'], ylab='Alpha-index',xlab=NA)
                histo <- thist(finalstate[,'Alphaindex'])
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


            

        }
