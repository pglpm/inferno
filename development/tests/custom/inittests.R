source('/home/pglpm/repos/inferno/R/util_buildauxmetadata.R')
source('/home/pglpm/repos/inferno/R/util_vtransform.R')
logprob <- function(data, metadata, initmethod, hyperparams = NULL){
    #### hyperparams
    if(is.null(hyperparams)){
        hyperparams <- list(
            ncomponents = 64,
            minalpha = -4,
            maxalpha = 4,
            byalpha = 1,
            Rshapelo = 0.5,
            Rshapehi = 0.5,
            Rvarm1 = 3^2,
            Cshapelo = 0.5,
            Cshapehi = 0.5,
            Cvarm1 = 3^2,
            Dshapelo = 0.5,
            Dshapehi = 0.5,
            Dvarm1 = 3^2,
            Bshapelo = 1,
            Bshapehi = 1,
            Dthreshold = 1,
            tscalefactor = 4.266,
            avoidzeroW = FALSE,
            initmethod = initmethod
            ## precluster, prior, allcentre
        )
    }

    for(aname in names(hyperparams)){
        assign(aname, hyperparams[[aname]])
    }


#### Metadata
    if (is.character(metadata) && file.exists(metadata)) {
        metadata <- read.csv(metadata,
            na.strings = '', stringsAsFactors = FALSE,
            colClasses=c(
                name = 'character',
                type = 'character',
                domainmin = 'numeric',
                domainmax = 'numeric',
                datastep = 'numeric',
                minincluded = 'character',
                maxincluded = 'character'
            ))
    }
    metadata <- as.data.frame(metadata)

    ## eliminate possible empty V-columns
    for(i in intersect(paste0('V', 11:3), colnames(metadata))){
        if(all(is.na(metadata[, i]))){
            metadata <- metadata[, -which(colnames(metadata) == i), drop=FALSE]
        }
    }

#### Dataset
    ## Check if 'data' is given
    if(!(missing(data) || is.null(data))) {
        ## Check if 'data' argument is an existing file
        ## otherwise we assume it is an object
        datafile <- NULL
        if (is.character(data)) {
            ## add '.csv' if missing
            datafile <- paste0(sub('.csv$', '', data), '.csv')
            if (file.exists(datafile)) {
                data <- read.csv(datafile,
                    na.strings = '', stringsAsFactors = FALSE, tryLogical = FALSE)
            } else {
                stop('Cannot find data file')
            }
        }
        data <- as.data.frame(data)
        rownames(data) <- NULL

        ## convert factors to strings if necessary
        if(any(sapply(data, is.factor))){
            cat('Converting factors to characters\n')
            . <- sapply(data, is.factor)
            data[, .] <- lapply(data[, ., drop = FALSE], as.character)
        }

        ## Consistency checks for data
        ## They should be moved to an external function

        ## Are data missing variates?
        if (!all(metadata[['name']] %in% colnames(data))) {
            stop('Missing variates in data. Check data- and metadata-files.')
        }

        ## Drop variates in data that are not in the metadata file
        if (!all(colnames(data) %in% metadata[['name']])) {
            cat('Warning: data have additional variates. Dropping them.\n')
            subvar <- intersect(colnames(data), metadata[['name']])
            data <- data[, subvar, drop = FALSE]
            rm(subvar)
        }

        ## Remove empty datapoints
        tokeep <- which(apply(data, 1, function(xx) { !all(is.na(xx)) }))
        if(length(tokeep) == 0 && !prior) {
            stop('Data are given but empty')
        } else if(length(tokeep) < nrow(data)) {
            cat('Warning: data contain empty datapoints. Dropping them.\n')
            data <- data[tokeep, , drop = FALSE]
        }
        rm(tokeep)

        npoints <- nrow(data)

        
    } else {
        ## data not given: we assume user wants prior calculation
        message('Missing data')
        prior <- TRUE
        npoints <- 0
    }


    auxmetadata <- buildauxmetadata(
        data = data,
        metadata = metadata,
        Dthreshold = hyperparams$Dthreshold,
        tscalefactor = hyperparams$tscalefactor
    )




    nalpha <- length(seq(minalpha, maxalpha, by = byalpha))

    vn <- vnames <- list(R=NULL, C=NULL, D=NULL, O=NULL, N=NULL, B=NULL)

    for (atype in names(vn)) {
        vnames[[atype]] <- auxmetadata[auxmetadata$mcmctype == atype, 'name']
        vn[[atype]] <- length(vnames[[atype]])
    }


#### CONSTANTS OF NIMBLE MODEL
    probalpha0 <- (1:nalpha)^2.25
    probalpha0 <- probalpha0/sum(probalpha0)
    constants <- c(
        list(
            ncomponents = ncomponents,
            npoints = npoints,
            nalpha = nalpha,
            alphabase = sqrt(2),
            probalpha0 = probalpha0,
            dirchalphas = rep((2^(minalpha - 0.5)) / ncomponents, ncomponents)
        ),
        if (vn$R > 0) { # continuous open domain
            list(
                Rn = vn$R, # This indexing variable is needed internally
                Rmean1 = rep(0, 1),
                Rvarm1 = rep(Rvarm1, 1),
                Rvar1 = rep(1, 1),
                Rshapelo = rep(Rshapelo, 1),
                Rshapehi = rep(Rshapehi, 1)
            )
        },
        if (vn$C > 0) { # continuous closed domain
            list(
                Cn = vn$C, # This indexing variable is needed internally
                Cmean1 = rep(0, 1),
                Cvarm1 = rep(Cvarm1, 1),
                Cvar1 = rep(1, 1),
                Cshapelo = rep(Cshapelo, 1),
                Cshapehi = rep(Cshapehi, 1),
                Cleft = as.matrix(vtransform(data[, vnames$C, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Cout = 'left', logjacobianOr = NULL)),
                Cright = as.matrix(vtransform(data[, vnames$C, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Cout = 'right', logjacobianOr = NULL)),
                Clatinit = as.matrix(vtransform(data[, vnames$C, drop = FALSE],
                    auxmetadata,
                    Cout = 'init', logjacobianOr = NULL))
            )
            ## Cleft & Cright are as many as the datapoints
            ## so we do not create copies outside of Nimble
            ## to save RAM
        },
        if (vn$D > 0) { # continuous rounded
            list(
                Dn = vn$D, # This indexing variable is needed internally
                Dmean1 = rep(0, 1),
                Dvarm1 = rep(Dvarm1, 1),
                Dvar1 = rep(1, 1),
                Dshapelo = rep(Dshapelo, 1),
                Dshapehi = rep(Dshapehi, 1),
                Dleft = as.matrix(vtransform(data[, vnames$D, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Dout = 'left', logjacobianOr = NULL)),
                Dright = as.matrix(vtransform(data[, vnames$D, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Dout = 'right', logjacobianOr = NULL)),
                Dlatinit = as.matrix(vtransform(data[, vnames$D, drop = FALSE],
                    auxmetadata,
                    Dout = 'init', logjacobianOr = NULL))
            )
        },
        ## if (vn$L > 0) { # latent
        ##     list(
        ##         Ln = vn$L, # This indexing variable is needed internally
        ##         Lmean1 = rep(0, 1),
        ##         Lvarm1 = rep(Lvarm1, 1),
        ##         Lvar1 = rep(1, 1),
        ##         Lshapelo = rep(Lshapelo, 1),
        ##         Lshapehi = rep(Lshapehi, 1),
        ##         Lleft = as.matrix(vtransform(data[, vnames$L, drop = FALSE],
        ##             auxmetadata = auxmetadata,
        ##             Lout = 'left')),
        ##         Lright = as.matrix(vtransform(data[, vnames$L, drop = FALSE],
        ##             auxmetadata = auxmetadata,
        ##             Lout = 'right')),
        ##         Llatinit = as.matrix(vtransform(data[, vnames$L, drop = FALSE],
        ##             auxmetadata,
        ##             Lout = 'init'))
        ##     )
        ## },
        if (vn$O > 0) { # ordinal
            Omaxn <- max(auxmetadata[auxmetadata$mcmctype == 'O', 'Nvalues'])
            Oalpha0 <- matrix(1e-100, nrow = vn$O, ncol = Omaxn)
            for (avar in seq_along(vnames$O)) {
                nvalues <- auxmetadata[auxmetadata$name == vnames$O[avar], 'Nvalues']
                ## ## we choose a flatter hyperprior for ordinal variates
                ## we choose a Hadamard-like hyperprior for nominal variates
                Oalpha0[avar, 1:nvalues] <- 1/nvalues
            }
            ##
            list(
                On = vn$O, # This indexing variable is needed internally
                Omaxn = Omaxn,
                Oalpha0 = Oalpha0
            )
        },
        if (vn$N > 0) { # nominal
            Nmaxn <- max(auxmetadata[auxmetadata$mcmctype == 'N', 'Nvalues'])
            Nalpha0 <- matrix(1e-100, nrow = vn$N, ncol = Nmaxn)
            for (avar in seq_along(vnames$N)) {
                nvalues <- auxmetadata[auxmetadata$name == vnames$N[avar], 'Nvalues']
                ## we choose a Hadamard-like hyperprior for nominal variates
                Nalpha0[avar, 1:nvalues] <- 1 / nvalues
            }
            ##
            list(
                Nn = vn$N, # This indexing variable is needed internally
                Nmaxn = Nmaxn,
                Nalpha0 = Nalpha0
            )
        },
        if (vn$B > 0) { # binary
            list(
                Bn = vn$B, # This indexing variable is needed internally
                Bshapelo = rep(Bshapelo, 1),
                Bshapehi = rep(Bshapehi, 1)
            )
        }
    ) # End constants

#### DATAPOINTS
    datapoints <- c(
        if (vn$R > 0) { # continuous open domain
            list(
                Rdata = as.matrix(vtransform(data[, vnames$R, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Rout = 'normalized', logjacobianOr = NULL))
            )
        },
        if (vn$C > 0) { # continuous closed domain
            list(
                Caux = as.matrix(vtransform(data[, vnames$C, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Cout = 'aux', logjacobianOr = NULL)),
                Clat = as.matrix(vtransform(data[, vnames$C, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Cout = 'lat', logjacobianOr = NULL))
            )
        },
        if (vn$D > 0) { # continuous rounded
            list(
                Daux = as.matrix(vtransform(data[, vnames$D, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Dout = 'aux', logjacobianOr = NULL))
            )
        },
        ## if (vn$L > 0) { # latent
        ##     list(
        ##         Laux = as.matrix(vtransform(data[, vnames$L, drop = FALSE],
        ##             auxmetadata = auxmetadata,
        ##             Lout = 'aux'))
        ##     )
        ## },
        if (vn$O > 0) { # nominal
            list(
                Odata = as.matrix(vtransform(data[, vnames$O, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Oout = 'numeric', logjacobianOr = NULL))
            )
        },
        if (vn$N > 0) { # nominal
            list(
                Ndata = as.matrix(vtransform(data[, vnames$N, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Nout = 'numeric', logjacobianOr = NULL))
            )
        },
        if (vn$B > 0) { # binary
            list(
                Bdata = as.matrix(vtransform(data[, vnames$B, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Bout = 'numeric', logjacobianOr = NULL))
            )
        }
    ) # End datapoints

####
#### INITIAL-VALUE FUNCTION
        ## init functions defined in 'util_mcmcinit.R'
    if(initmethod == 'precluster'){
                    initsfn <- function() {
                ## Create components centres
                ## distance function
                ## NB: all variances will be initialized to 1
                lpnorm <- function(xx){abs(xx)}
                distances <- matrix(0, nrow = constants$npoints,
                    ncol = constants$ncomponents)
                if (vn$R > 0) { # continuous open domain
                    Rmeans <- matrix(rnorm(
                        n = vn$R * constants$ncomponents,
                        mean = constants$Rmean1,
                        sd = sqrt(constants$Rvarm1)
                    ), nrow = vn$R, ncol = constants$ncomponents)
                    Rvars <- matrix(1, nrow = vn$R, ncol = constants$ncomponents)
                    ## distances from datapoints
                    distances <- distances + apply(Rmeans, 2, function(ameans){
                        colSums(lpnorm(t(datapoints$Rdata) - ameans), na.rm = TRUE)
                    })
                }
                if (vn$C > 0) { # continuous closed domain
                    Cmeans <- matrix(rnorm(
                        n = vn$C * constants$ncomponents,
                        mean = constants$Cmean1,
                        sd = sqrt(constants$Cvarm1)
                    ), nrow = vn$C, ncol = constants$ncomponents)
                    Cvars <- matrix(1, nrow = vn$C, ncol = constants$ncomponents)
                    ## distances from datapoints
                    distances <- distances + apply(Cmeans, 2, function(ameans){
                        colSums(lpnorm(t(datapoints$Clat) - ameans), na.rm = TRUE)
                    })
                }
                if (vn$D > 0) { # discrete
                    Dmeans <- matrix(rnorm(
                        n = vn$D * constants$ncomponents,
                        mean = constants$Dmean1,
                        sd = sqrt(constants$Dvarm1)
                    ), nrow = vn$D, ncol = constants$ncomponents)
                    Dvars <- matrix(1, nrow = vn$D, ncol = constants$ncomponents)
                    ## distances from datapoints
                    distances <- distances + apply(Dmeans, 2, function(ameans){
                        colSums(lpnorm(t(constants$Dlatinit) - ameans), na.rm = TRUE)
                    })
                }
                ## if (vn$L > 0) { # 
                ##     Lmeans <- matrix(rnorm(
                ##         n = vn$L * constants$ncomponents,
                ##         mean = constants$Lmean1,
                ##         sd = sqrt(constants$Lvarm1)
                ##     ), nrow = vn$L, ncol = constants$ncomponents)
                ##     ## distances from datapoints
                ##     distances <- distances + apply(Lmeans, 2, function(ameans){
                ##         colSums(lpnorm(t(constants$Llatinit) - ameans), na.rm = TRUE)
                ##     })
                ## }
                ## if (vn$B > 0) {
                ##     Bprobs <- matrix(rbeta(
                ##         n = vn$B * constants$ncomponents,
                ##         shape1 = Bshapelo,
                ##         shape2 = Bshapehi,
                ##         ), nrow = vn$B, ncol = constants$ncomponents)
                ##     ## distances from datapoints
                ##     distances <- distances + apply(Bprobs, 2, function(ameans){
                ##         colSums(lpnorm(t(datapoints$Bdata) - ameans), na.rm = TRUE)
                ##     })
                ## }

                ## assign datapoints to component with closest centre
                K <- apply(distances, 1, which.min)
                occupied <- unique(K)

                ## recalculate components centres according to their points
                if (vn$R > 0) {
                    Rmeans[, occupied] <- sapply(occupied, function(acomponent){
                        colMeans(
                            datapoints$Rdata[which(K == acomponent), , drop = FALSE],
                            na.rm = TRUE)
                    })
                    Rmeans[, -occupied] <- 0
                    Rvars[, occupied] <- sapply(occupied, function(acomponent){
                        apply(
                            datapoints$Rdata[which(K == acomponent), , drop = FALSE],
                            2, sd, na.rm = TRUE)
                    })^2
                    Rvars[is.na(Rvars)] <- 1
                }
                if (vn$C > 0) {
                    Cmeans[, occupied] <- sapply(occupied, function(acomponent){
                        colMeans(
                            datapoints$Clat[which(K == acomponent), , drop = FALSE],
                            na.rm = TRUE)
                    })
                    Cmeans[, -occupied] <- 0
                    Cvars[, occupied] <- sapply(occupied, function(acomponent){
                        apply(
                            datapoints$Clat[which(K == acomponent), , drop = FALSE],
                            2, sd, na.rm = TRUE)
                    })^2
                    Cvars[is.na(Cvars)] <- 1
                }
                if (vn$D > 0) {
                    Dmeans[, occupied] <- sapply(occupied, function(acomponent){
                        colMeans(
                            constants$Dlatinit[which(K == acomponent), , drop = FALSE],
                            na.rm = TRUE)
                    })
                    Dmeans[, -occupied] <- 0
                    Dvars[, occupied] <- sapply(occupied, function(acomponent){
                        apply(constants$Dlatinit[which(K == acomponent), , drop = FALSE],
                            2, sd, na.rm = TRUE)
                    })^2
                    Dvars[is.na(Dvars)] <- 1
                }
                ## if (vn$L > 0) { # continuous open domain
                ##     Lmeans[, occupied] <- sapply(occupied, function(acomponent){
                ##         colMeans(constants$Llatinit[which(K == acomponent), , drop = FALSE],
                ##             na.rm = TRUE)
                ##     })
                ##     Lmeans[, -occupied] <- 0
                ## }
                ## if (vn$B > 0) {
                ##     Bprobs[, occupied] <- sapply(occupied, function(acomponent){
                ##         colMeans(datapoints$Bdata[which(K == acomponent), , drop = FALSE],
                ##             na.rm = TRUE)
                ##     })
                ##     Bprobs[, -occupied] <- 0.5
                ## }
                ## Alpha <- sample(1:constants$nalpha, 1, prob = constants$probalpha0, replace = TRUE)
                ## W <- c(rep(rempoints, miconstants$npoints), rep(1, constants$ncomponents - miconstants$npoints))
                ## W <- W/sum(W)

                outlist <- list(
                    Alpha = round(constants$nalpha/2),
                    W = rep(1/constants$ncomponents, constants$ncomponents),
                    ## ## Assign every point to the closest component centre
                    K = K
                    ## ## Other assignment methods:
                    ## ## A. assign all points to an unsystematically chosen component
                    ## K = rep(sample(rep(which(W > 0), 2), 1), nepoints)
                    ## ## B. distribute points unsystematically among components
                    ## K = sample(rep(which(W > 0), 2), constants$npoints, replace = TRUE)
                    ## ## or:
                    ## ## C. assign all points to the most probable component
                    ## K = rep(which.max(W), constants$npoints)
                    ## ## or:
                    ## ## D. assign all points to the least probable component
                    ## K = rep(which(W == min(W[W > 0]))[1], constants$npoints)
                    ## ## or:
                    ## ## E. distribute points unsystematically among M=2 components
                    ## K = sample(sample(rep(which(W > 0), 2), 2, replace = TRUE),
                    ##           constants$npoints, replace = TRUE)
                    ## ## F. mix methods A. and B.
                    ## K = (if(achain %% 2 == Ksample) {
                    ##        ## ## assign all points to an unsystematically chosen component
                    ##        rep(sample(rep(which(W > 0), 2), 1), constants$npoints)
                    ##      } else {
                    ##        ## distribute points unsystematically among components
                    ##        sample(rep(which(W > 0), 2), constants$npoints, replace = TRUE)
                    ##      })
                )
                ##
                if (vn$R > 0) { # continuous open domain
                    outlist <- c(
                        outlist,
                        list(
                            Rmean = Rmeans,
                            Rrate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Rshapehi,
                                    rate = constants$Rvar1),
                                nrow = vn$R, ncol = constants$ncomponents
                            ),
                            Rvar = Rvars
                        )
                    )
                }
                if (vn$C > 0) { # continuous closed domain
                    outlist <- c(
                        outlist,
                        list(
                            Cmean = Cmeans,
                            Crate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Cshapehi,
                                    rate = constants$Cvar1),
                                nrow = vn$C, ncol = constants$ncomponents
                            ),
                            Cvar = Cvars,
                            ## for data with boundary values
                            Clat = constants$Clatinit
                            ## Clat = vtransform(data[, vnames$C, with = FALSE],
                            ##   auxmetadata, Cout = 'init')
                        )
                    )
                }
                if (vn$D > 0) { # continuous rounded
                    outlist <- c(
                        outlist,
                        list(
                            Dmean = Dmeans,
                            Drate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Dshapehi,
                                    rate = constants$Dvar1),
                                nrow = vn$D, ncol = constants$ncomponents
                            ),
                            Dvar = Dvars,
                            ## for data with boundary values
                            Dlat = constants$Dlatinit
                            ## Dlat = vtransform(data[, vnames$D, with = FALSE],
                            ##   auxmetadata, Dout = 'init')
                        )
                    )
                }
                ## if (vn$L > 0) { # latent
                ##     outlist <- c(
                ##         outlist,
                ##         list(
                ##             Lmean = Lmeans,
                ##             Lrate = matrix(
                ##                 nimble::qinvgamma(p = 0.5,
                ##                     shape = constants$Lshapehi,
                ##                     rate = constants$Lvar1),
                ##                 nrow = vn$L, ncol = constants$ncomponents
                ##             ),
                ##             Lvar = matrix(1,
                ##                 nrow = vn$L, ncol = constants$ncomponents),
                ##             ## for data with boundary values
                ##             Llat = constants$Llatinit
                ##             ## Llat = vtransform(data[, vnames$L, with = FALSE],
                ##             ##   auxmetadata, Lout = 'init')
                ##         )
                ##     )
                ## }
                if (vn$O > 0) { # ordinal
                    outlist <- c(
                        outlist,
                        list(
                            Oprob = aperm(array(sapply(1:vn$O, function(avar) {
                                sapply(1:constants$ncomponents, function(aclus) {
                                    Oalpha0[avar, ]/sum(Oalpha0[avar, ])
                                    ## nimble::rdirch(n = 1, alpha = Oalpha0[avar, ])
                                })
                            }), dim = c(Omaxn, constants$ncomponents, vn$O)))
                        )
                    )
                }
                if (vn$N > 0) { # nominal
                    outlist <- c(
                        outlist,
                        list(
                            Nprob = aperm(array(sapply(1:vn$N, function(avar) {
                                sapply(1:constants$ncomponents, function(aclus) {
                                    constants$Nalpha0[avar, ]/sum(constants$Nalpha0[avar, ])
                                    ## nimble::rdirch(n = 1, alpha = constants$Nalpha0[avar, ])
                                })
                            }), dim = c(Nmaxn, constants$ncomponents, vn$N)))
                        )
                    )
                }
                if (vn$B > 0) { # binary
                    outlist <- c(
                        outlist,
                        list(
                            ## Bprob = Bprobs
                            Bprob = matrix(0.5, nrow = vn$B, ncol = constants$ncomponents)
                        )
                    )
                }
                ##
                outlist
            } # end precluster

        } else if(initmethod == 'precluster0'){
            ## pre-clustering, k-mean style
            initsfn <- function() {
                ## Create components centres
                ## distance function
                ## NB: all variances will be initialized to 1
                lpnorm <- function(xx){abs(xx)}
                distances <- matrix(0, nrow = constants$npoints,
                    ncol = constants$ncomponents)
                if (vn$R > 0) { # continuous open domain
                    Rmeans <- matrix(rnorm(
                        n = vn$R * constants$ncomponents,
                        mean = constants$Rmean1,
                        sd = sqrt(constants$Rvarm1)
                    ), nrow = vn$R, ncol = constants$ncomponents)
                    ## distances from datapoints
                    distances <- distances + apply(Rmeans, 2, function(ameans){
                        colSums(lpnorm(t(datapoints$Rdata) - ameans), na.rm = TRUE)
                    })
                }
                if (vn$C > 0) { # continuous closed domain
                    Cmeans <- matrix(rnorm(
                        n = vn$C * constants$ncomponents,
                        mean = constants$Cmean1,
                        sd = sqrt(constants$Cvarm1)
                    ), nrow = vn$C, ncol = constants$ncomponents)
                    ## distances from datapoints
                    distances <- distances + apply(Cmeans, 2, function(ameans){
                        colSums(lpnorm(t(datapoints$Clat) - ameans), na.rm = TRUE)
                    })
                }
                if (vn$D > 0) { # discrete
                    Dmeans <- matrix(rnorm(
                        n = vn$D * constants$ncomponents,
                        mean = constants$Dmean1,
                        sd = sqrt(constants$Dvarm1)
                    ), nrow = vn$D, ncol = constants$ncomponents)
                    ## distances from datapoints
                    distances <- distances + apply(Dmeans, 2, function(ameans){
                        colSums(lpnorm(t(constants$Dlatinit) - ameans), na.rm = TRUE)
                    })
                }
                ## if (vn$L > 0) { # 
                ##     Lmeans <- matrix(rnorm(
                ##         n = vn$L * constants$ncomponents,
                ##         mean = constants$Lmean1,
                ##         sd = sqrt(constants$Lvarm1)
                ##     ), nrow = vn$L, ncol = constants$ncomponents)
                ##     ## distances from datapoints
                ##     distances <- distances + apply(Lmeans, 2, function(ameans){
                ##         colSums(lpnorm(t(constants$Llatinit) - ameans), na.rm = TRUE)
                ##     })
                ## }
                ## if (vn$B > 0) {
                ##     Bprobs <- matrix(rbeta(
                ##         n = vn$B * constants$ncomponents,
                ##         shape1 = Bshapelo,
                ##         shape2 = Bshapehi,
                ##         ), nrow = vn$B, ncol = constants$ncomponents)
                ##     ## distances from datapoints
                ##     distances <- distances + apply(Bprobs, 2, function(ameans){
                ##         colSums(lpnorm(t(datapoints$Bdata) - ameans), na.rm = TRUE)
                ##     })
                ## }

                ## assign datapoints to component with closest centre
                K <- apply(distances, 1, which.min)
                occupied <- unique(K)

                ## recalculate components centres according to their points
                if (vn$R > 0) {
                    Rmeans[, occupied] <- sapply(occupied, function(acomponent){
                        colMeans(datapoints$Rdata[which(K == acomponent), , drop = FALSE],
                            na.rm = TRUE)
                    })
                    Rmeans[, -occupied] <- 0
                }
                if (vn$C > 0) {
                    Cmeans[, occupied] <- sapply(occupied, function(acomponent){
                        colMeans(datapoints$Clat[which(K == acomponent), , drop = FALSE],
                            na.rm = TRUE)
                    })
                    Cmeans[, -occupied] <- 0
                }
                if (vn$D > 0) {
                    Dmeans[, occupied] <- sapply(occupied, function(acomponent){
                        colMeans(constants$Dlatinit[which(K == acomponent), , drop = FALSE],
                            na.rm = TRUE)
                    })
                    Dmeans[, -occupied] <- 0
                }
                ## if (vn$L > 0) { # continuous open domain
                ##     Lmeans[, occupied] <- sapply(occupied, function(acomponent){
                ##         colMeans(constants$Llatinit[which(K == acomponent), , drop = FALSE],
                ##             na.rm = TRUE)
                ##     })
                ##     Lmeans[, -occupied] <- 0
                ## }
                ## if (vn$B > 0) {
                ##     Bprobs[, occupied] <- sapply(occupied, function(acomponent){
                ##         colMeans(datapoints$Bdata[which(K == acomponent), , drop = FALSE],
                ##             na.rm = TRUE)
                ##     })
                ##     Bprobs[, -occupied] <- 0.5
                ## }
                ## Alpha <- sample(1:constants$nalpha, 1, prob = constants$probalpha0, replace = TRUE)
                ## W <- c(rep(rempoints, miconstants$npoints), rep(1, constants$ncomponents - miconstants$npoints))
                ## W <- W/sum(W)

                outlist <- list(
                    Alpha = round(constants$nalpha/2),
                    W = rep(1/constants$ncomponents, constants$ncomponents),
                    ## ## Assign every point to the closest component centre
                    K = K
                    ## ## Other assignment methods:
                    ## ## A. assign all points to an unsystematically chosen component
                    ## K = rep(sample(rep(which(W > 0), 2), 1), nepoints)
                    ## ## B. distribute points unsystematically among components
                    ## K = sample(rep(which(W > 0), 2), constants$npoints, replace = TRUE)
                    ## ## or:
                    ## ## C. assign all points to the most probable component
                    ## K = rep(which.max(W), constants$npoints)
                    ## ## or:
                    ## ## D. assign all points to the least probable component
                    ## K = rep(which(W == min(W[W > 0]))[1], constants$npoints)
                    ## ## or:
                    ## ## E. distribute points unsystematically among M=2 components
                    ## K = sample(sample(rep(which(W > 0), 2), 2, replace = TRUE),
                    ##           constants$npoints, replace = TRUE)
                    ## ## F. mix methods A. and B.
                    ## K = (if(achain %% 2 == Ksample) {
                    ##        ## ## assign all points to an unsystematically chosen component
                    ##        rep(sample(rep(which(W > 0), 2), 1), constants$npoints)
                    ##      } else {
                    ##        ## distribute points unsystematically among components
                    ##        sample(rep(which(W > 0), 2), constants$npoints, replace = TRUE)
                    ##      })
                )
                ##
                if (vn$R > 0) { # continuous open domain
                    outlist <- c(
                        outlist,
                        list(
                            Rmean = Rmeans,
                            Rrate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Rshapehi,
                                    rate = constants$Rvar1),
                                nrow = vn$R, ncol = constants$ncomponents
                            ),
                            Rvar = matrix(1,
                                nrow = vn$R, ncol = constants$ncomponents)
                        )
                    )
                }
                if (vn$C > 0) { # continuous closed domain
                    outlist <- c(
                        outlist,
                        list(
                            Cmean = Cmeans,
                            Crate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Cshapehi,
                                    rate = constants$Cvar1),
                                nrow = vn$C, ncol = constants$ncomponents
                            ),
                            Cvar = matrix(1,
                                nrow = vn$C, ncol = constants$ncomponents),
                            ## for data with boundary values
                            Clat = constants$Clatinit
                            ## Clat = vtransform(data[, vnames$C, with = FALSE],
                            ##   auxmetadata, Cout = 'init')
                        )
                    )
                }
                if (vn$D > 0) { # continuous rounded
                    outlist <- c(
                        outlist,
                        list(
                            Dmean = Dmeans,
                            Drate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Dshapehi,
                                    rate = constants$Dvar1),
                                nrow = vn$D, ncol = constants$ncomponents
                            ),
                            Dvar = matrix(1,
                                nrow = vn$D, ncol = constants$ncomponents),
                            ## for data with boundary values
                            Dlat = constants$Dlatinit
                            ## Dlat = vtransform(data[, vnames$D, with = FALSE],
                            ##   auxmetadata, Dout = 'init')
                        )
                    )
                }
                ## if (vn$L > 0) { # latent
                ##     outlist <- c(
                ##         outlist,
                ##         list(
                ##             Lmean = Lmeans,
                ##             Lrate = matrix(
                ##                 nimble::qinvgamma(p = 0.5,
                ##                     shape = constants$Lshapehi,
                ##                     rate = constants$Lvar1),
                ##                 nrow = vn$L, ncol = constants$ncomponents
                ##             ),
                ##             Lvar = matrix(1,
                ##                 nrow = vn$L, ncol = constants$ncomponents),
                ##             ## for data with boundary values
                ##             Llat = constants$Llatinit
                ##             ## Llat = vtransform(data[, vnames$L, with = FALSE],
                ##             ##   auxmetadata, Lout = 'init')
                ##         )
                ##     )
                ## }
                if (vn$O > 0) { # ordinal
                    outlist <- c(
                        outlist,
                        list(
                            Oprob = aperm(array(sapply(1:vn$O, function(avar) {
                                sapply(1:constants$ncomponents, function(aclus) {
                                    Oalpha0[avar, ]/sum(Oalpha0[avar, ])
                                    ## nimble::rdirch(n = 1, alpha = Oalpha0[avar, ])
                                })
                            }), dim = c(Omaxn, constants$ncomponents, vn$O)))
                        )
                    )
                }
                if (vn$N > 0) { # nominal
                    outlist <- c(
                        outlist,
                        list(
                            Nprob = aperm(array(sapply(1:vn$N, function(avar) {
                                sapply(1:constants$ncomponents, function(aclus) {
                                    constants$Nalpha0[avar, ]/sum(constants$Nalpha0[avar, ])
                                    ## nimble::rdirch(n = 1, alpha = constants$Nalpha0[avar, ])
                                })
                            }), dim = c(Nmaxn, constants$ncomponents, vn$N)))
                        )
                    )
                }
                if (vn$B > 0) { # binary
                    outlist <- c(
                        outlist,
                        list(
                            ## Bprob = Bprobs
                            Bprob = matrix(0.5, nrow = vn$B, ncol = constants$ncomponents)
                        )
                    )
                }
                ##
                outlist
            } # end precluster


        } else if(initmethod == 'prior'){
            ## values chosen from prior
            initsfn <- function() {
                Alpha <- sample(1:constants$nalpha, 1, prob = probalpha0[1:constants$nalpha])
                W <- nimble::rdirch(n = 1,
                    alpha = constants$dirchalphas[1:constants$ncomponents] *
                        constants$alphabase^Alpha)
                outlist <- list(
                    Alpha = Alpha,
                    W = W,
                    K = sample(rep(which(W > 0), 2), constants$npoints, replace = TRUE)
                )
                ##
                if (vn$R > 0) { # continuous open domain
                    Rrate <- matrix(
                        nimble::rinvgamma(
                            n = vn$R * constants$ncomponents,
                            shape = constants$Rshapehi,
                            rate = constants$Rvar1),
                        nrow = vn$R, ncol = constants$ncomponents
                    )
                    outlist <- c(
                        outlist,
                        list(
                            Rmean = matrix(rnorm(
                                n = vn$R * constants$ncomponents,
                                mean = constants$Rmean1,
                                sd = sqrt(constants$Rvarm1)
                            ), nrow = vn$R, ncol = constants$ncomponents),
                            Rrate = Rrate,
                            Rvar = matrix(
                                nimble::rinvgamma(
                                    n = vn$R * constants$ncomponents,
                                    shape = constants$Rshapelo,
                                    rate = Rrate),
                                nrow = vn$R, ncol = constants$ncomponents
                            )
                        )
                    )
                }
                if (vn$C > 0) { # continuous closed domain
                    Crate <- matrix(
                        nimble::rinvgamma(
                            n = vn$C * constants$ncomponents,
                            shape = constants$Cshapehi,
                            rate = constants$Cvar1),
                        nrow = vn$C, ncol = constants$ncomponents
                    )
                    outlist <- c(
                        outlist,
                        list(
                            Cmean = matrix(rnorm(
                                n = vn$C * constants$ncomponents,
                                mean = constants$Cmean1,
                                sd = sqrt(constants$Cvarm1)
                            ), nrow = vn$C, ncol = constants$ncomponents),
                            Crate = Crate,
                            Cvar = matrix(
                                nimble::rinvgamma(
                                    n = vn$C * constants$ncomponents,
                                    shape = constants$Cshapelo,
                                    rate = Crate),
                                nrow = vn$C, ncol = constants$ncomponents
                            ),
                            Clat = constants$Clatinit
                        )
                    )
                }
                if (vn$D > 0) { # continuous rounded
                    Drate <- matrix(
                        nimble::rinvgamma(
                            n = vn$D * constants$ncomponents,
                            shape = constants$Dshapehi,
                            rate = constants$Dvar1),
                        nrow = vn$D, ncol = constants$ncomponents
                    )
                    outlist <- c(
                        outlist,
                        list(
                            Dmean = matrix(rnorm(
                                n = vn$D * constants$ncomponents,
                                mean = constants$Dmean1,
                                sd = sqrt(constants$Dvarm1)
                            ), nrow = vn$D, ncol = constants$ncomponents),
                            Drate = Drate,
                            Dvar = matrix(
                                nimble::rinvgamma(
                                    n = vn$D * constants$ncomponents,
                                    shape = constants$Dshapelo,
                                    rate = Drate),
                                nrow = vn$D, ncol = constants$ncomponents
                            ),
                            Dlat = constants$Dlatinit
                        )
                    )
                }
                if (vn$O > 0) { # ordinal
                    outlist <- c(
                        outlist,
                        list(
                            Oprob = aperm(array(sapply(1:vn$O, function(avar) {
                                sapply(1:constants$ncomponents, function(aclus) {
                                    nimble::rdirch(n = 1, alpha = Oalpha0[avar, ])
                                })
                            }), dim = c(Omaxn, constants$ncomponents, vn$O)))
                        )
                    )
                }
                if (vn$N > 0) { # nominal
                    outlist <- c(
                        outlist,
                        list(
                            Nprob = aperm(array(sapply(1:vn$N, function(avar) {
                                sapply(1:constants$ncomponents, function(aclus) {
                                    nimble::rdirch(n = 1, alpha = constants$Nalpha0[avar, ])
                                })
                            }), dim = c(Nmaxn, constants$ncomponents, vn$N)))
                        )
                    )
                }
                if (vn$B > 0) { # nominal
                    outlist <- c(
                        outlist,
                        list(
                            Bprob = matrix(
                                rbeta(n = vn$B * constants$ncomponents,
                                    shape1 = Bshapelo, shape2 = Bshapehi),
                                nrow = vn$B, ncol = constants$ncomponents)
                        )
                    )
                }
                ##
                outlist
            } # end prior

        } else if(initmethod == 'priorbutfirst'){
            ## values chosen from prior
            ## except first, fit to data
            initsfn <- function() {
                Alpha <- sample(1:constants$nalpha, 1, prob = probalpha0[1:constants$nalpha])
                W <- nimble::rdirch(n = 1,
                    alpha = constants$dirchalphas[1:constants$ncomponents] *
                        constants$alphabase^Alpha)
                outlist <- list(
                    Alpha = Alpha,
                    W = W,
                    K = rep(1, constants$npoints)
                )
                ##
                if (vn$R > 0) { # continuous open domain
                    Rrate <- matrix(
                        nimble::rinvgamma(
                            n = vn$R * constants$ncomponents,
                            shape = constants$Rshapehi,
                            rate = constants$Rvar1),
                        nrow = vn$R, ncol = constants$ncomponents
                    )
                    outlist <- c(
                        outlist,
                        list(
                            Rmean = matrix(
                                c(rep(0, vn$R), rnorm(
                                n = vn$R * (constants$ncomponents -1),
                                mean = constants$Rmean1,
                                sd = sqrt(constants$Rvarm1)
                            )), nrow = vn$R, ncol = constants$ncomponents),
                            Rrate = Rrate,
                            Rvar = matrix(
                                c(rep(1, vn$R), nimble::rinvgamma(
                                    n = vn$R * (constants$ncomponents - 1),
                                    shape = constants$Rshapelo,
                                    rate = Rrate[-1, ])),
                                nrow = vn$R, ncol = constants$ncomponents
                            )
                        )
                    )
                }
                if (vn$C > 0) { # continuous closed domain
                    Crate <- matrix(
                        nimble::rinvgamma(
                            n = vn$C * constants$ncomponents,
                            shape = constants$Cshapehi,
                            rate = constants$Cvar1),
                        nrow = vn$C, ncol = constants$ncomponents
                    )
                    outlist <- c(
                        outlist,
                        list(
                            Cmean = matrix(
                                c(rep(0, vn$C), rnorm(
                                n = vn$C * (constants$ncomponents - 1),
                                mean = constants$Cmean1,
                                sd = sqrt(constants$Cvarm1)
                            )), nrow = vn$C, ncol = constants$ncomponents),
                            Crate = Crate,
                            Cvar = matrix(
                                c(rep(1, vn$C), nimble::rinvgamma(
                                    n = vn$C * (constants$ncomponents - 1),
                                    shape = constants$Cshapelo,
                                    rate = Crate[-1, ])),
                                nrow = vn$C, ncol = constants$ncomponents
                            ),
                            Clat = constants$Clatinit
                        )
                    )
                }
                if (vn$D > 0) { # continuous rounded
                    Drate <- matrix(
                        nimble::rinvgamma(
                            n = vn$D * constants$ncomponents,
                            shape = constants$Dshapehi,
                            rate = constants$Dvar1),
                        nrow = vn$D, ncol = constants$ncomponents
                    )
                    outlist <- c(
                        outlist,
                        list(
                            Dmean = matrix(
                                c(rep(0, vn$D), rnorm(
                                n = vn$D * (constants$ncomponents - 1),
                                mean = constants$Dmean1,
                                sd = sqrt(constants$Dvarm1)
                            )), nrow = vn$D, ncol = constants$ncomponents),
                            Drate = Drate,
                            Dvar = matrix(
                                c(rep(1, vn$D), nimble::rinvgamma(
                                    n = vn$D * (constants$ncomponents - 1),
                                    shape = constants$Dshapelo,
                                    rate = Drate[-1, ])),
                                nrow = vn$D, ncol = constants$ncomponents
                            ),
                            Dlat = constants$Dlatinit
                        )
                    )
                }
                if (vn$O > 0) { # ordinal
                    outlist <- c(
                        outlist,
                        list(
                            Oprob = aperm(array(sapply(1:vn$O, function(avar) {
                                sapply(1:constants$ncomponents, function(aclus) {
                                    nimble::rdirch(n = 1, alpha = Oalpha0[avar, ])
                                })
                            }), dim = c(Omaxn, constants$ncomponents, vn$O)))
                        )
                    )
                }
                if (vn$N > 0) { # nominal
                    outlist <- c(
                        outlist,
                        list(
                            Nprob = aperm(array(sapply(1:vn$N, function(avar) {
                                sapply(1:constants$ncomponents, function(aclus) {
                                    nimble::rdirch(n = 1, alpha = constants$Nalpha0[avar, ])
                                })
                            }), dim = c(Nmaxn, constants$ncomponents, vn$N)))
                        )
                    )
                }
                if (vn$B > 0) { # nominal
                    outlist <- c(
                        outlist,
                        list(
                            Bprob = matrix(
                                rbeta(n = vn$B * constants$ncomponents,
                                    shape1 = Bshapelo, shape2 = Bshapehi),
                                nrow = vn$B, ncol = constants$ncomponents)
                        )
                    )
                }
                ##
                outlist
            } # end priorbutfirst


        } else if(initmethod == 'allcentre'){
            ## all components equal, all points to first component
            initsfn <- function() {
                initvar <- 6
                outlist <- list(
                    Alpha = round(constants$nalpha/2),
                    W = rep(1/constants$ncomponents, constants$ncomponents),
                    ## ## Assign every point to first component
                    K = rep(1, constants$npoints)
                )
                ##
                if (vn$R > 0) { # continuous open domain
                    outlist <- c(
                        outlist,
                        list(
                            Rmean =  matrix(0, nrow = vn$R, ncol = constants$ncomponents),
                            Rrate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Rshapehi,
                                    rate = constants$Rvar1),
                                nrow = vn$R, ncol = constants$ncomponents
                            ),
                            Rvar = matrix(initvar, nrow = vn$R, ncol = constants$ncomponents)
                        )
                    )
                }
                if (vn$C > 0) { # continuous closed domain
                    outlist <- c(
                        outlist,
                        list(
                            Cmean = matrix(0, nrow = vn$C, ncol = constants$ncomponents),
                            Crate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Cshapehi,
                                    rate = constants$Cvar1),
                                nrow = vn$C, ncol = constants$ncomponents
                            ),
                            Cvar = matrix(initvar, nrow = vn$C, ncol = constants$ncomponents),
                            ## for data with boundary values
                            Clat = constants$Clatinit
                        )
                    )
                }
                if (vn$D > 0) { # continuous rounded
                    outlist <- c(
                        outlist,
                        list(
                            Dmean = matrix(0, nrow = vn$D, ncol = constants$ncomponents),
                            Drate = matrix(
                                nimble::qinvgamma(p = 0.5,
                                    shape = constants$Dshapehi,
                                    rate = constants$Dvar1),
                                nrow = vn$D, ncol = constants$ncomponents
                            ),
                            Dvar = matrix(initvar, nrow = vn$D, ncol = constants$ncomponents),
                            ## for data with boundary values
                            Dlat = constants$Dlatinit
                        )
                    )
                }
                if (vn$O > 0) { # ordinal
                    outlist <- c(
                        outlist,
                        list(
                            Oprob = aperm(array(sapply(1:vn$O, function(avar) {
                                sapply(1:constants$ncomponents, function(aclus) {
                                    Oalpha0[avar, ]/sum(Oalpha0[avar, ])
                                })
                            }), dim = c(Omaxn, constants$ncomponents, vn$O)))
                        )
                    )
                }
                if (vn$N > 0) { # nominal
                    outlist <- c(
                        outlist,
                        list(
                            Nprob = aperm(array(sapply(1:vn$N, function(avar) {
                                sapply(1:constants$ncomponents, function(aclus) {
                                    constants$Nalpha0[avar, ]/sum(constants$Nalpha0[avar, ])
                                })
                            }), dim = c(Nmaxn, constants$ncomponents, vn$N)))
                        )
                    )
                }
                if (vn$B > 0) { # binary
                    outlist <- c(
                        outlist,
                        list(
                            Bprob = matrix(0.5, nrow = vn$B, ncol = constants$ncomponents)
                        )
                    )
                }
                ##
                outlist
            } # end allcentre

        } else if(initmethod == 'allinone'){
            ## Components from prior, data all in one
            initsfn <- function() {
                outlist <- list(
                    Alpha = round(constants$nalpha/2),
                    W = rep(1/constants$ncomponents, constants$ncomponents),
                    ## ## Assign every point to first component
                    K = rep(1, constants$npoints)
                )
                ##
                if (vn$R > 0) { # continuous open domain
                    Rrate <- matrix(
                        nimble::rinvgamma(
                            n = vn$R,
                            shape = constants$Rshapehi,
                            rate = constants$Rvar1),
                        nrow = vn$R, ncol = constants$ncomponents
                    )
                    outlist <- c(
                        outlist,
                        list(
                            Rmean = matrix(rnorm(
                                n = vn$R,
                                mean = constants$Rmean1,
                                sd = sqrt(constants$Rvarm1)
                            ), nrow = vn$R, ncol = constants$ncomponents),
                            Rrate = Rrate,
                            Rvar = matrix(
                                nimble::rinvgamma(
                                    n = vn$R,
                                    shape = constants$Rshapelo,
                                    rate = Rrate),
                                nrow = vn$R, ncol = constants$ncomponents
                            )
                        )
                    )
                }
                if (vn$C > 0) { # continuous closed domain
                    Crate <- matrix(
                        nimble::rinvgamma(
                            n = vn$C,
                            shape = constants$Cshapehi,
                            rate = constants$Cvar1),
                        nrow = vn$C, ncol = constants$ncomponents
                    )
                    outlist <- c(
                        outlist,
                        list(
                            Cmean = matrix(rnorm(
                                n = vn$C,
                                mean = constants$Cmean1,
                                sd = sqrt(constants$Cvarm1)
                            ), nrow = vn$C, ncol = constants$ncomponents),
                            Crate = Crate,
                            Cvar = matrix(
                                nimble::rinvgamma(
                                    n = vn$C,
                                    shape = constants$Cshapelo,
                                    rate = Crate),
                                nrow = vn$C, ncol = constants$ncomponents
                            ),
                            Clat = constants$Clatinit
                        )
                    )
                }
                if (vn$D > 0) { # continuous rounded
                    Drate <- matrix(
                        nimble::rinvgamma(
                            n = vn$D,
                            shape = constants$Dshapehi,
                            rate = constants$Dvar1),
                        nrow = vn$D, ncol = constants$ncomponents
                    )
                    outlist <- c(
                        outlist,
                        list(
                            Dmean = matrix(rnorm(
                                n = vn$D,
                                mean = constants$Dmean1,
                                sd = sqrt(constants$Dvarm1)
                            ), nrow = vn$D, ncol = constants$ncomponents),
                            Drate = Drate,
                            Dvar = matrix(
                                nimble::rinvgamma(
                                    n = vn$D,
                                    shape = constants$Dshapelo,
                                    rate = Drate),
                                nrow = vn$D, ncol = constants$ncomponents
                            ),
                            Dlat = constants$Dlatinit
                        )
                    )
                }
                if (vn$O > 0) { # ordinal
                    outlist <- c(
                        outlist,
                        list(
                            Oprob = aperm(array(sapply(1:vn$O, function(avar) {
                                ## sapply(1:constants$ncomponents, function(aclus) {
                                nimble::rdirch(n = 1, alpha = Oalpha0[avar, ])
                                ## })
                            }), dim = c(Omaxn, vn$O, constants$ncomponents)), c(2, 3, 1))
                        )
                    )
                }
                if (vn$N > 0) { # nominal
                    outlist <- c(
                        outlist,
                        list(
                            Nprob = aperm(array(sapply(1:vn$N, function(avar) {
                                ## sapply(1:constants$ncomponents, function(aclus) {
                                nimble::rdirch(n = 1, alpha = constants$Nalpha0[avar, ])
                                ## })
                            }), dim = c(Nmaxn, vn$N, constants$ncomponents)), c(2, 3, 1))
                        )
                    )
                }
                if (vn$B > 0) { # nominal
                    outlist <- c(
                        outlist,
                        list(
                            Bprob = matrix(
                                rbeta(n = vn$B,
                                    shape1 = Bshapelo, shape2 = Bshapehi),
                                nrow = vn$B, ncol = constants$ncomponents)
                        )
                    )
                }
                ##
                outlist
            } # end allinone

        } else {stop('Unknown "initmethod"')}

####
    
    vals <- with(c(hyperparams, constants), {
        initsfn()
    })
    ## print(str(vals))

    out <- with(c(constants,vals,datapoints), {
        ## print(
        ##     sum(na.rm=TRUE, sapply(1:npoints, function(d){
        ##     ##     nimble::dcat(K[d] , prob = W[1:ncomponents], log=TRUE)
        ##     ## }))
        ##     ## +
        ##             ##
        ##             ## (if (vn$R > 0) { # continuous open domain
        ##             ##     sum(na.rm=TRUE, sapply(1:Rn, function(v){
        ##             ##         dnorm(Rdata[d, v], mean = Rmean[v, K[d]],  sd = sqrt(Rvar[v, K[d]]), log=TRUE)
        ##             ##     }))}else{0}) +
        ##             ## (if (vn$C > 0) { # continuous closed domain
        ##             ##     sum(na.rm=TRUE, sapply(1:Cn, function(v){
        ##             ##         nimble::dconstraint(Caux[d, v], Clat[d, v] >= Cleft[d, v] &
        ##             ##                                             Clat[d, v] <= Cright[d, v], log=TRUE) +
        ##             ##             dnorm(Clat[d, v], mean = Cmean[v, K[d]],  sd = sqrt(Cvar[v, K[d]]), log=TRUE)
        ##             ##     }))}else{0}) +
        ##             (if (vn$D > 0) { # continuous rounded
        ##                 sum(na.rm=TRUE, sapply(1:Dn, function(v){
        ##                     nimble::dconstraint(Daux[d, v], Dlat[d, v] >= Dleft[d, v] &
        ##                                                         Dlat[d, v] < Dright[d, v], log=TRUE) +
        ##                         dnorm(Dlat[d, v], mean = Dmean[v, K[d]],  sd = sqrt(Dvar[v, K[d]]), log=TRUE)
        ##                 }))}else{0})
        ##     ##         +
        ##     ##         (if (vn$O > 0) { # nominal
        ##     ##             sum(na.rm=TRUE, sapply(1:On, function(v){
        ##     ##                 nimble::dcat(Odata[d, v], prob = Oprob[v, K[d], 1:Omaxn], log=TRUE)
        ##     ##             }))}else{0}) +
        ##     ##         (if (vn$N > 0) { # nominal
        ##     ##             sum(na.rm=TRUE, sapply(1:Nn, function(v){
        ##     ##                 nimble::dcat(Ndata[d, v], prob = Nprob[v, K[d], 1:Nmaxn], log=TRUE)
        ##     ##             }))}else{0}) +
        ##     ##         (if (vn$B > 0) { # binary
        ##     ##             sum(na.rm=TRUE, sapply(1:Bn, function(v){
        ##     ##                 dbinom(Bdata[d, v], size = 1, prob = Bprob[v, K[d]], log=TRUE) 
        ##     ##             }))}else{0})
        ##     }))
        ##     )
        nimble::dcat(Alpha, prob = probalpha0[1:nalpha], log=TRUE) +
            nimble::ddirch(W,
                alpha = dirchalphas[1:ncomponents] * alphabase^Alpha,
                log=TRUE) +
            sum(na.rm=TRUE, sapply(1:ncomponents, function(k) {
                (if (vn$R > 0) { # continuous open domain
                    sum(na.rm=TRUE, sapply(1:Rn, function(v){
                        dnorm(Rmean[v, k] , mean = Rmean1,  sd = sqrt(Rvarm1), log=TRUE) +
                            dinvgamma(Rrate[v, k] , shape = Rshapehi, rate = Rvar1, log=TRUE) +
                            dinvgamma(Rvar[v, k] , shape = Rshapelo, rate = Rrate[v, k], log=TRUE)
                    }))}else{0}) +
                    (if (vn$C > 0) { # continuous closed domain
                        sum(na.rm=TRUE, sapply(1:Cn, function(v){
                            dnorm(Cmean[v, k] , mean = Cmean1,  sd = sqrt(Cvarm1), log=TRUE) +
                                nimble::dinvgamma(Crate[v, k] , shape = Cshapehi, rate = Cvar1, log=TRUE) +
                                nimble::dinvgamma(Cvar[v, k] , shape = Cshapelo, rate = Crate[v, k], log=TRUE)
                        }))}else{0}) +
                    (if (vn$D > 0) { # continuous rounded
                        sum(na.rm=TRUE, sapply(1:Dn, function(v){
                            dnorm(Dmean[v, k] , mean = Dmean1,  sd = sqrt(Dvarm1), log=TRUE) +
                                nimble::dinvgamma(Drate[v, k] , shape = Dshapehi, rate = Dvar1, log=TRUE) +
                                nimble::dinvgamma(Dvar[v, k] , shape = Dshapelo, rate = Drate[v, k], log=TRUE)
                        }))}else{0}) +
                    (if (vn$O > 0) { # ordinal
                        sum(na.rm=TRUE, sapply(1:On, function(v){
                            nimble::ddirch(Oprob[v, k, 1:Omaxn] , alpha = Oalpha0[v, 1:Omaxn], log=TRUE)
                        }))}else{0}) +
                    (if (vn$N > 0) { # ordinal
                        sum(na.rm=TRUE, sapply(1:Nn, function(v){
                            nimble::ddirch(Nprob[v, k, 1:Nmaxn] , alpha = Nalpha0[v, 1:Nmaxn], log=TRUE)
                        }))}else{0}) +
                    (if (vn$B > 0) { # ordinal
                        sum(na.rm=TRUE, sapply(1:Bn, function(v){
                            dbeta(Bprob[v, k] , shape1 = Bshapelo, shape2 = Bshapehi, log=TRUE)
                        }))}else{0})
            })) +
            ##
            sum(na.rm=TRUE, sapply(1:npoints, function(d){
                nimble::dcat(K[d] , prob = W[1:ncomponents], log=TRUE) +
                    ##
                    (if (vn$R > 0) { # continuous open domain
                        sum(na.rm=TRUE, sapply(1:Rn, function(v){
                            dnorm(Rdata[d, v], mean = Rmean[v, K[d]],  sd = sqrt(Rvar[v, K[d]]), log=TRUE)
                        }))}else{0}) +
                    (if (vn$C > 0) { # continuous closed domain
                        sum(na.rm=TRUE, sapply(1:Cn, function(v){
                            nimble::dconstraint(Caux[d, v], Clat[d, v] >= Cleft[d, v] &
                                                                Clat[d, v] <= Cright[d, v], log=TRUE) +
                                dnorm(Clat[d, v], mean = Cmean[v, K[d]],  sd = sqrt(Cvar[v, K[d]]), log=TRUE)
                        }))}else{0}) +
                    (if (vn$D > 0) { # continuous rounded
                        sum(na.rm=TRUE, sapply(1:Dn, function(v){
                            nimble::dconstraint(Daux[d, v], Dlat[d, v] >= Dleft[d, v] &
                                                                Dlat[d, v] < Dright[d, v], log=TRUE) +
                                dnorm(Dlat[d, v], mean = Dmean[v, K[d]],  sd = sqrt(Dvar[v, K[d]]), log=TRUE)
                        }))}else{0}) +
                    (if (vn$O > 0) { # nominal
                        sum(na.rm=TRUE, sapply(1:On, function(v){
                            nimble::dcat(Odata[d, v], prob = Oprob[v, K[d], 1:Omaxn], log=TRUE)
                        }))}else{0}) +
                    (if (vn$N > 0) { # nominal
                        sum(na.rm=TRUE, sapply(1:Nn, function(v){
                            nimble::dcat(Ndata[d, v], prob = Nprob[v, K[d], 1:Nmaxn], log=TRUE)
                        }))}else{0}) +
                    (if (vn$B > 0) { # binary
                        sum(na.rm=TRUE, sapply(1:Bn, function(v){
                            dbinom(Bdata[d, v], size = 1, prob = Bprob[v, K[d]], log=TRUE) 
                        }))}else{0})
            }))
    })
    return(out)
}


## > inim <- 'allcentre' ; summary(sapply(1:100, function(ii){max(sapply(1:10,function(i)logprob(data='/home/pglpm/repos/inferno/development/downloads/penguin_data.csv', metadata='/home/pglpm/repos/inferno/development/downloads/penguin_metadata.csv', initmethod=inim)))}))
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -18543  -11807  -10639  -10975   -9792   -7866
## 
## > inim <- 'prior' ; summary(sapply(1:100, function(ii){max(sapply(1:10,function(i)logprob(data='/home/pglpm/repos/inferno/development/downloads/penguin_data.csv', metadata='/home/pglpm/repos/inferno/development/downloads/penguin_metadata.csv', initmethod=inim)))}))
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -31174  -12945  -10978  -12044  -10119   -8330
## 
## > inim <- 'allinone' ; summary(sapply(1:100, function(ii){max(sapply(1:10,function(i)logprob(data='/home/pglpm/repos/inferno/development/downloads/penguin_data.csv', metadata='/home/pglpm/repos/inferno/development/downloads/penguin_metadata.csv', initmethod=inim)))}))
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -24717  -12004  -10615  -11113   -9803   -8251
## 
## > inim <- 'precluster' ; summary(sapply(1:100, function(ii){max(sapply(1:10,function(i)logprob(data='/home/pglpm/repos/inferno/development/downloads/penguin_data.csv', metadata='/home/pglpm/repos/inferno/development/downloads/penguin_metadata.csv', initmethod=inim)))}))
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  -28001  -11789  -10827  -11147  -10016   -8049 





## > inim <- 'precluster' ; summary(foreach(ii=1:360, .combine=c, .inorder = F)%dorng%{max(sapply(1:10,function(i)logprob(data='/home/pglpm/repos/inferno/development/downloads/penguin_data.csv', metadata='/home/pglpm/repos/inferno/development/downloads/penguin_metadata.csv', initmethod=inim)))})
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  -44986  -12270  -10963  -11649   -9894   -7570 
## >
## > inim <- 'prior' ; summary(foreach(ii=1:360, .combine=c, .inorder = F)%dorng%{max(sapply(1:10,function(i)logprob(data='/home/pglpm/repos/inferno/development/downloads/penguin_data.csv', metadata='/home/pglpm/repos/inferno/development/downloads/penguin_metadata.csv', initmethod=inim)))})
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  -25088  -11902  -10628  -11210   -9694   -7895 
## >
## > inim <- 'allinone' ; summary(foreach(ii=1:360, .combine=c, .inorder = F)%dorng%{max(sapply(1:10,function(i)logprob(data='/home/pglpm/repos/inferno/development/downloads/penguin_data.csv', metadata='/home/pglpm/repos/inferno/development/downloads/penguin_metadata.csv', initmethod=inim)))})
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  -29009  -12051  -10725  -11330   -9834   -7988 
## 
## > inim <- 'allcentre' ; summary(foreach(ii=1:360, .combine=c, .inorder = F)%dorng%{max(sapply(1:10,function(i)logprob(data='/home/pglpm/repos/inferno/development/downloads/penguin_data.csv', metadata='/home/pglpm/repos/inferno/development/downloads/penguin_metadata.csv', initmethod=inim)))})
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  -27662  -12021  -10733  -11432   -9826   -8139 
## > 




#### varm 8, tsf 1

## > inim <- 'prior' ; . <- (foreach(ii=1:360, .combine=c, .inorder = F)%dorng%{max(sapply(1:10,function(i)logprob(data='/home/pglpm/repos/inferno/development/downloads/penguin_data.csv', metadata='/home/pglpm/repos/inferno/development/downloads/penguin_metadata.csv', initmethod=inim)))}) ; summary(.) ; summary(.[is.finite(.)])
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -5406854      Inf      Inf      Inf      Inf      Inf 
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -5406854 -1768044 -1566540 -1882457 -1017649  -580003 
## > sum(!is.finite(.))
## [1] 349

## > inim <- 'precluster' ; summary(foreach(ii=1:360, .combine=c, .inorder = F)%dorng%{max(sapply(1:10,function(i)logprob(data='/home/pglpm/repos/inferno/development/downloads/penguin_data.csv', metadata='/home/pglpm/repos/inferno/development/downloads/penguin_metadata.csv', initmethod=inim)))})
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   -6460   -6356   -6346   -6349   -6337   -6306 

## > inim <- 'allinone' ; summary(foreach(ii=1:360, .combine=c, .inorder = F)%dorng%{max(sapply(1:10,function(i)logprob(data='/home/pglpm/repos/inferno/development/downloads/penguin_data.csv', metadata='/home/pglpm/repos/inferno/development/downloads/penguin_metadata.csv', initmethod=inim)))})
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -108175  -22705  -16265  -19936  -12584   -8280 

## > inim <- 'allcentre' ; summary(foreach(ii=1:360, .combine=c, .inorder = F)%dorng%{max(sapply(1:10,function(i)logprob(data='/home/pglpm/repos/inferno/development/downloads/penguin_data.csv', metadata='/home/pglpm/repos/inferno/development/downloads/penguin_metadata.csv', initmethod=inim)))})
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   -8157   -8157   -8157   -8157   -8157   -8157 

## > inim <- 'priorbutfirst' ; . <- (foreach(ii=1:360, .combine=c, .inorder = F)%dorng%{max(sapply(1:10,function(i)logprob(data='/home/pglpm/repos/inferno/development/downloads/penguin_data.csv', metadata='/home/pglpm/repos/inferno/development/downloads/penguin_metadata.csv', initmethod=inim)))}) ; summary(.) ; summary(.[is.finite(.)])
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##  -28022     Inf     Inf     Inf     Inf     Inf      89 
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  -28022  -21096  -15041  -16497  -12619   -8079 


#### varm 8, tsf 1.35
## > inim <- 'precluster' ; summary(foreach(ii=1:360, .combine=c, .inorder = F)%dorng%{max(sapply(1:10,function(i)logprob(data='/home/pglpm/repos/inferno/development/downloads/penguin_data.csv', metadata='/home/pglpm/repos/inferno/development/downloads/penguin_metadata.csv', initmethod=inim)))})
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   -6590   -6435   -6414   -6417   -6394   -6340 
