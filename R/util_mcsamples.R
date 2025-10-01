#' Worker function called by learn()
#'
#' This worker function is defined outside of learn.R in order to avoid import of spurious objects into the parallel workers, with waste of memory
#'
#' @keywords internal
workerfun <- function(
    acore,
    dirname,
    dashnameroot,
    avoidzeroW,
    initmethod,
    constants,
    datapoints,
    vn,
    showAlphatraces,
    Alphatoslice,
    Ktoslice,
    RWtoslice,
    changeSamplerOrder,
    minchainspercore,
    coreswithextrachain,
    nchains,
    maxhours,
    timestart0,
    showsamplertimes,
    startupMCiterations,
    maxMCiterations,
    showKtraces,
    ncomponents,
    plottraces,
    Qlo,
    Qhi,
    Qerror,
    minESS,
    initES,
    nsamplesperchain,
    minMCiterations,
    printtimediff
) {
    ## functions to format printing of time
    printtimeend <- function(tim) {
        format(Sys.time() + tim, format='%Y-%m-%d %H:%M')
    }
    ## We need to send some messages to the log files, others to the user.
    ## This is done by changing output sink:
    print2user <- function(msg, outcon) {
        flush(outcon)
        sink(file = NULL, type = 'message')
        message(msg, appendLF = FALSE)
        sink(file = outcon, type = 'message')
    }

    ## Create log file
    ## Redirect diagnostics and service messages there
    outcon <- file(file.path(dirname,
        paste0('log', dashnameroot,
            '-', acore, '.log')
    ), open = 'w')
    sink(file = outcon, type = 'output')
    sink(file = outcon, type = 'message')
    ## Leave this FALSE to bypass message bug in doParallel
    if(TRUE){
        closecons <- function(){
            ## Close output to log files
            flush(outcon)
            sink(file = NULL, type = 'output')
            sink(file = NULL, type = 'message')
            close(outcon)
        }
        on.exit(closecons())
    }

    usedmem <- sum(gc()[,6])

    ## Timer
    headertimestart <- Sys.time()

    cat('Log core', acore)
    cat(' - Current time:',
        strftime(as.POSIXlt(headertimestart), '%Y-%m-%d %H:%M:%S'))
    cat('\n')

    suppressPackageStartupMessages(require('nimble'))
    ## requireNamespace("nimble", quietly = TRUE)
    ##library('nimble')


#### Alternative Dirichlet sampler, to deal with zero-W preserving conjugacy
#### from Daniel Turek 2025-06-29, nimble-users group
    if(is.null(avoidzeroW)) {
        conjugate_dirch_plus_eps <- nimbleFunction(
            contains = sampler_BASE,
            setup = function(model, mvSaved, target, control) {
                ## control list extraction
                eps <- extractControlElement(control,
                    'eps', .Machine$double.xmin)
                ## node list generation
                target <- model$expandNodeNames(target)
                calcNodes <- model$getDependencies(target)
                depNodes <- model$getDependencies(target, self = FALSE)
                ## numeric value generation
                d <- length(model[[target]])
                Ndeps <- length(depNodes)
                ## checks
                if(length(target) > 1) stop('conjugate_dirch_plus_eps only operates on a single node')
                if(model$getDistribution(target) != 'ddirch') stop('conjugate_dirch_plus_eps target node must have ddirch distribution')
                depDists <- sapply(depNodes, function(x) model$getDistribution(x))
                if(!all(depDists == 'dcat')) stop('conjugate_dirch_plus_eps all dependencies should have dcat distributions')
            },
            run = function() {
                posterior_alpha <- model$getParam(target, 'alpha')
                for(i in 1:Ndeps) {
                    depValue <- model$getParam(depNodes[i], 'value')
                    posterior_alpha[depValue] <- posterior_alpha[depValue] + 1
                }
                newTargetValue <- rdirch(1, posterior_alpha)
                newTargetValue <- newTargetValue + eps
                model[[target]] <<- newTargetValue
                model$calculate(calcNodes)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
            },
            methods = list(
                reset = function() {}
            )
        )
    }

#### COMPONENT REPRESENTATION OF FREQUENCY SPACE
#### Dirichlet-process mixture of product-kernels

    ## hierarchical probability structure
    finitemix <- nimbleCode({
        ## Component weights
        Alpha ~ dcat(prob = probalpha0[1:nalpha])
        alphas[1:ncomponents] <- dirchalphas[1:ncomponents] * alphabase^Alpha
        W[1:ncomponents] ~ ddirch(alpha = alphas[1:ncomponents])
        if(isTRUE(avoidzeroW)){
            W0[1:ncomponents] <- W[1:ncomponents] + epsd
        }
        ## Probability density for the parameters of the components
        for (k in 1:ncomponents) {
            ## Probability distributions of parameters
            ## of the different variate types
            if (vn$R > 0) { # continuous open domain
                for (v in 1:Rn) {
                    Rmean[v, k] ~ dnorm(mean = Rmean1, var = Rvarm1)
                    Rrate[v, k] ~ dinvgamma(shape = Rshapehi, rate = Rvar1)
                    Rvar[v, k] ~ dinvgamma(shape = Rshapelo, rate = Rrate[v, k])
                }
            }
            if (vn$C > 0) { # continuous closed domain
                for (v in 1:Cn) {
                    Cmean[v, k] ~ dnorm(mean = Cmean1, var = Cvarm1)
                    Crate[v, k] ~ dinvgamma(shape = Cshapehi, rate = Cvar1)
                    Cvar[v, k] ~ dinvgamma(shape = Cshapelo, rate = Crate[v, k])
                }
            }
            if (vn$D > 0) { # continuous rounded
                for (v in 1:Dn) {
                    Dmean[v, k] ~ dnorm(mean = Dmean1, var = Dvarm1)
                    Drate[v, k] ~ dinvgamma(shape = Dshapehi, rate = Dvar1)
                    Dvar[v, k] ~ dinvgamma(shape = Dshapelo, rate = Drate[v, k])
                }
            }
            ## if (vn$L > 0) { # latent
            ##     for (v in 1:Ln) {
            ##         Lmean[v, k] ~ dnorm(mean = Lmean1, var = Lvarm1)
            ##         Lrate[v, k] ~ dinvgamma(shape = Lshapehi, rate = Lvar1)
            ##         Lvar[v, k] ~ dinvgamma(shape = Lshapelo, rate = Lrate[v, k])
            ##     }
            ## }
            if (vn$O > 0) { # ordinal
                for (v in 1:On) {
                    Oprob[Oi[v]:Of[v], k] ~ ddirch(alpha = Oalpha0[Oi[v]:Of[v]])
                }
            }
            if (vn$N > 0) { # nominal
                for (v in 1:Nn) {
                    Nprob[Ni[v]:Nf[v], k] ~ ddirch(alpha = Nalpha0[Ni[v]:Nf[v]])
                }
            }
            if (vn$B > 0) { # binary
                for (v in 1:Bn) {
                    Bprob[v, k] ~ dbeta(shape1 = Bshapelo, shape2 = Bshapehi)
                }
            }
        }
        ## Probability of data
        for (d in 1:npoints) {
            if(isTRUE(avoidzeroW)){
                K[d] ~ dcat(prob = W0[1:ncomponents])
            } else {
                K[d] ~ dcat(prob = W[1:ncomponents])
            }
            ##
            if (vn$R > 0) { # continuous open domain
                for (v in 1:Rn) {
                    Rdata[d, v] ~ dnorm(mean = Rmean[v, K[d]], var = Rvar[v, K[d]])
                }
            }
            if (vn$C > 0) { # continuous closed domain
                for (v in 1:Cn) {
                    Caux[d, v] ~ dconstraint(Clat[d, v] >= Cleft[d, v] &
                                                 Clat[d, v] <= Cright[d, v])
                    Clat[d, v] ~ dnorm(mean = Cmean[v, K[d]], var = Cvar[v, K[d]])
                }
            }
            if (vn$D > 0) { # continuous rounded
                for (v in 1:Dn) {
                    Daux[d, v] ~ dconstraint(Dlat[d, v] >= Dleft[d, v] &
                                                 Dlat[d, v] < Dright[d, v])
                    Dlat[d, v] ~ dnorm(mean = Dmean[v, K[d]], var = Dvar[v, K[d]])
                }
            }
            ## if (vn$L > 0) { # latent
            ##     for (v in 1:Ln) {
            ##         Laux[d, v] ~ dconstraint(Llat[d, v] >= Lleft[d, v] &
            ##                                  Llat[d, v] < Lright[d, v])
            ##         Llat[d, v] ~ dnorm(mean = Lmean[v, K[d]], var = Lvar[v, K[d]])
            ##     }
            ## }
            if (vn$O > 0) { # nominal
                for (v in 1:On) {
                    Odata[d, v] ~ dcat(prob = Oprob[Oi[v]:Of[v], K[d]])
                }
            }
            if (vn$N > 0) { # nominal
                for (v in 1:Nn) {
                    Ndata[d, v] ~ dcat(prob = Nprob[Ni[v]:Nf[v], K[d]])
                }
            }
            if (vn$B > 0) { # binary
                for (v in 1:Bn) { # Bprob is the probability that Bdata=1
                    Bdata[d, v] ~ dbern(prob = Bprob[v, K[d]])
                }
            }
        }
    }) # end finitemix NimbleCode


#### INITIAL-VALUE FUNCTION
    ## init functions defined in 'util_mcmcinit.R'
    if(initmethod == 'precluster'){
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
                ## Rvars[, occupied] <- sapply(occupied, function(acomponent){
                ##     apply(
                ##         datapoints$Rdata[which(K == acomponent), , drop = FALSE],
                ##         2, sd, na.rm = TRUE)
                ## })^2
                ## Rvars[is.na(Rvars)] <- 1
            }
            if (vn$C > 0) {
                Cmeans[, occupied] <- sapply(occupied, function(acomponent){
                    colMeans(
                        datapoints$Clat[which(K == acomponent), , drop = FALSE],
                        na.rm = TRUE)
                })
                Cmeans[, -occupied] <- 0
                ## Cvars[, occupied] <- sapply(occupied, function(acomponent){
                ##     apply(
                ##         datapoints$Clat[which(K == acomponent), , drop = FALSE],
                ##         2, sd, na.rm = TRUE)
                ## })^2
                ## Cvars[is.na(Cvars)] <- 1
            }
            if (vn$D > 0) {
                Dmeans[, occupied] <- sapply(occupied, function(acomponent){
                    colMeans(
                        constants$Dlatinit[which(K == acomponent), , drop = FALSE],
                        na.rm = TRUE)
                })
                Dmeans[, -occupied] <- 0
                ## Dvars[, occupied] <- sapply(occupied, function(acomponent){
                ##     apply(constants$Dlatinit[which(K == acomponent), , drop = FALSE],
                ##         2, sd, na.rm = TRUE)
                ## })^2
                ## Dvars[is.na(Dvars)] <- 1
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
                        Oprob = matrix(unlist(sapply(Ocards, function(acard){
                            rep(1 / acard, acard)
                        })), nrow = sum(Ocards), ncol = ncomponents)
                    )
                )
            }
            if (vn$N > 0) { # nominal
                outlist <- c(
                    outlist,
                    list(
                        Nprob = matrix(unlist(sapply(Ncards, function(acard){
                            rep(1 / acard, acard)
                        })), nrow = sum(Ncards), ncol = ncomponents)
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

    } else if(initmethod == 'datacentre'){
        ## pre-clustering, centering on data
        initsfn <- function() {
            ## assign each cluster to a datapoint
            iK <- sample(seq_len(constants$npoints), constants$ncomponents,
                replace = (constants$npoints < constants$ncomponents))

            ## Create components centres
            ## distance function
            ## NB: all variances will be initialized to 1
            lpnorm <- function(xx){abs(xx)}
            distances <- matrix(0, nrow = constants$npoints,
                ncol = constants$ncomponents)
            if (vn$R > 0) { # continuous open domain
                Rmeans <- t(datapoints$Rdata[iK, , drop = FALSE])
                ## distances from datapoints
                distances <- distances + apply(Rmeans, 2, function(ameans){
                    colSums(lpnorm(t(datapoints$Rdata) - ameans), na.rm = TRUE)
                })
            }
            if (vn$C > 0) { # continuous closed domain
                Cmeans <- t(datapoints$Clat[iK, , drop = FALSE])
                ## distances from datapoints
                distances <- distances + apply(Cmeans, 2, function(ameans){
                    colSums(lpnorm(t(datapoints$Clat) - ameans), na.rm = TRUE)
                })
            }
            if (vn$D > 0) { # discrete
                Dmeans <- t(constants$Dlatinit[iK, , drop = FALSE])
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
                        ## Clat = inferno:::vtransform(data[, vnames$C, with = FALSE],
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
                        ## Dlat = inferno:::vtransform(data[, vnames$D, with = FALSE],
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
            ##             ## Llat = inferno:::vtransform(data[, vnames$L, with = FALSE],
            ##             ##   auxmetadata, Lout = 'init')
            ##         )
            ##     )
            ## }
            if (vn$O > 0) { # ordinal
                outlist <- c(
                    outlist,
                    list(
                        Oprob = matrix(unlist(sapply(Ocards, function(acard){
                            rep(1 / acard, acard)
                        })), nrow = sum(Ocards), ncol = ncomponents)
                    )
                )
            }
            if (vn$N > 0) { # nominal
                outlist <- c(
                    outlist,
                    list(
                        Nprob = matrix(unlist(sapply(Ncards, function(acard){
                            rep(1 / acard, acard)
                        })), nrow = sum(Ncards), ncol = ncomponents)
                    )
                )
            }
            if (vn$B > 0) { # binary
                outlist <- c(
                    outlist,
                    list(
                        ## Bprob = Bprobs
                        Bprob = matrix(0.5,
                            nrow = vn$B, ncol = constants$ncomponents)
                    )
                )
            }
            ##
            outlist
        } # end datacentre

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
                        Oprob = matrix(unlist(sapply(Ocards, function(acard){
                            rep(1 / acard, acard)
                        })), nrow = sum(Ocards), ncol = ncomponents)
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
                        Oprob = matrix(unlist(sapply(Ocards, function(acard){
                            rep(1 / acard, acard)
                        })), nrow = sum(Ocards), ncol = ncomponents)
                    )
                )
            }
            if (vn$N > 0) { # nominal
                outlist <- c(
                    outlist,
                    list(
                        Nprob = matrix(unlist(sapply(Ncards, function(acard){
                            rep(1 / acard, acard)
                        })), nrow = sum(Ncards), ncol = ncomponents)
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
                        Oprob = matrix(unlist(sapply(seq_along(Ocards),
                            function(v){nimble::rdirch(n = 1,
                                alpha = Oalpha0[Oi[v]:Of[v]])
                            })), nrow = sum(Ocards), ncol = ncomponents)
                    )
                )
            }
            if (vn$N > 0) { # nominal
                outlist <- c(
                    outlist,
                    list(
                        Nprob = matrix(unlist(sapply(seq_along(Ncards),
                            function(v){nimble::rdirch(n = 1,
                                alpha = Nalpha0[Ni[v]:Nf[v]])
                            })), nrow = sum(Ncards), ncol = ncomponents)
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

#################################################
#### NIMBLE SETUP
##################################################

    finitemixnimble <- nimbleModel(
        code = finitemix,
        name = 'finitemixnimble1',
        constants = constants,
        data = datapoints,
        ## dimensions = c(
        ##     ##     if (vn$R > 0) {
        ##     ##         list(
        ##     ##             Rdata = c(npoints, vn$R)
        ##     ##             )
        ##     ##     },
        ##     ##     if (vn$C > 0) {
        ##     ##         list(
        ##     ##             Caux = c(npoints, vn$C),
        ##     ##             Clat = c(npoints, vn$C),
        ##     ##             Cleft = c(npoints, vn$C),
        ##     ##             Cright = c(npoints, vn$C),
        ##     ##             Clatinit = c(npoints, vn$C)
        ##     ##             )
        ##     ##     },
        ##     ##     if (vn$D > 0) {
        ##     ##         list(
        ##     ##             Daux = c(npoints, vn$D),
        ##     ##             Dleft = c(npoints, vn$D),
        ##     ##             Dright = c(npoints, vn$D),
        ##     ##             Dlatinit = c(npoints, vn$D)
        ##     ##             )
        ##     ##     },
        ##     ## if (vn$O > 0) {
        ##     ##     list(
        ##     ##         Oi = length(Oi),
        ##     ##         Of = length(Of)
        ##     ##         ##             Odata = c(npoints, vn$O)
        ##     ##     )
        ##     ## },
        ##     ## if (vn$N > 0) {
        ##     ##     list(
        ##     ##         Ni = length(Ni),
        ##     ##         Nf = length(Nf)
        ##     ##         ##             Ndata = c(npoints, vn$N)
        ##     ##     )
        ##     ## }
        ##     ##     if (vn$B > 0) {
        ##     ##         list(
        ##     ##             Bdata = c(npoints, vn$B)
        ##     ##             )
        ##     ##     }
        ## ),
        inits = initsfn()
    )

    Cfinitemixnimble <- compileNimble(finitemixnimble,
        showCompilerOutput = FALSE)
    usedmem <- max(usedmem, sum(gc()[,6])) #garbage collection

    confnimble <- configureMCMC(
        Cfinitemixnimble, # nodes = NULL
        monitors = c(
            'W',
            if (vn$R > 0) {
                c('Rmean', 'Rvar')
            },
            if (vn$C > 0) {
                c('Cmean', 'Cvar')
            },
            if (vn$D > 0) {
                c('Dmean', 'Dvar')
            },
            ## if (vn$L > 0) {
            ##     c('Lmean', 'Lvar')
            ## },
            if (vn$O > 0) {
                c('Oprob')
            },
            if (vn$N > 0) {
                c('Nprob')
            },
            if (vn$B > 0) {
                c('Bprob')
            }
        ),
        ## It is necessary to monitor K to see if all components were used
        ## if 'showAlphatraces' is true,
        ## then the Alpha-parameter trace is also recorded and shown
        monitors2 = c(if (showAlphatraces) { 'Alpha' },
            'K')
    )

    ## ## Uncomment to debug Nimble (in case of Nimble updates)
    ## print(confnimble$getUnsampledNodes())
    ## cat('\nEX1\n')
    ## confnimble$printSamplers(executionOrder=TRUE)

    targetslist <- sapply(confnimble$getSamplers(), function(xx) xx$target)
    nameslist <- sapply(confnimble$getSamplers(), function(xx) xx$name)
    ## cat('\nNAMESLIST', nameslist, '\n')

    ## replace W's sampler with the one from Daniel Turek
    if(is.null(avoidzeroW)) {
        confnimble$replaceSampler('W', conjugate_dirch_plus_eps)
    }

    ## replace Alpha's cat-sampler with slice
    if (Alphatoslice &&
            !('Alpha' %in% targetslist[nameslist == 'posterior_predictive'])) {
        confnimble$replaceSampler(target='Alpha', type='slice')
        ## ## Old replacement method, didn't work in previous Nimble
        ## confnimble$removeSamplers('Alpha')
        ## confnimble$addSampler(target = 'Alpha', type = 'slice')
    }

    ## replace K's cat-sampler with slice
    if (Ktoslice) {
        for (asampler in grep('^K\\[', targetslist, value = TRUE)) {
            if (!(asampler %in% targetslist[nameslist == 'posterior_predictive'])) {
                confnimble$replaceSampler(target=asampler, type='slice')
                ## ## Old replacement method, didn't work in previous Nimble
                ## confnimble$removeSamplers(asampler)
                ## confnimble$addSampler(target = asampler, type = 'slice')
            }
        }
    }

    ## replace all RW samplers with slice
    ## Should find a way to do this faster
    ## testreptime <- Sys.time()
    if (RWtoslice) {
        for (asampler in targetslist[nameslist == 'RW']) {
            ## ## New replacement method, didn't work in previous Nimble
            confnimble$replaceSampler(target=asampler, type='slice')
            ## ## Old replacement method:
            ## confnimble$removeSamplers(asampler)
            ## confnimble$addSampler(target = asampler, type = 'slice')
        }
    }

    ## ## Uncomment when debugging Nimble
    ## print(confnimble$getUnsampledNodes())
    ## cat('\nEX2\n')
    ## confnimble$printSamplers(executionOrder=TRUE)

#### change execution order for some variates
    if (changeSamplerOrder) {
        ## call this to do a first reordering of the samplers
        ## it places posterior-predictive nodes last
        mcsampler <- buildMCMC(confnimble)

        samplerorder <- c(
            'K',
            'W',
            'Alpha',
            if (vn$R > 0) {
                c('Rmean', 'Rrate', 'Rvar')
            },
            if (vn$C > 0) {
                c('Cmean', 'Crate', 'Cvar')
            },
            if (vn$D > 0) {
                c('Dmean', 'Drate', 'Dvar')
            },
            ## if (vn$L > 0) {
            ##     c('Lmean', 'Lrate', 'Lvar')
            ## },
            if (vn$O > 0) {
                c('Oprob')
            },
            if (vn$N > 0) {
                c('Nprob')
            },
            if (vn$B > 0) {
                c('Bprob')
            }
        )
        ##
        neworder <- unlist(lapply(samplerorder, function(var){
            grep(
                paste0('^', var, '(\\[.+\\])*$'),
                sapply(confnimble$getSamplers(), function(x) {
                    if (!(x$name == 'posterior_predictive')) {
                        x$target
                    } else {
                        NULL
                    }
                })
            )
        }))
        ## ## Uncomment for debugging
        ## cat('\nNEW ORDER',neworder,'\n')
        confnimble$setSamplerExecutionOrder(c(setdiff(
            confnimble$getSamplerExecutionOrder(),
            neworder
        ), neworder))
        ## cat('\nEX3\n')
        ## confnimble$printSamplers(executionOrder=TRUE)

    }

#### Compile Monte Carlo sampler
    print(confnimble)
    mcsampler <- buildMCMC(confnimble)
    ## print(confnimble$getUnsampledNodes())
    Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

    cat('\nSetup time',
        printtimediff(difftime(Sys.time(), headertimestart, units = 'auto')),
        '\n')
    headertime <- difftime(Sys.time(), headertimestart, units = 'secs')


    ## Inform user that compilation is done, if core 1:
    if (acore == 1) {
        print2user(paste0('\rCompiled core ', acore, '. ',
            'Number of samplers: ',
            length(confnimble$samplerExecutionOrder), '.             \n',
            'Estimating remaining time, please be patient...'),
            outcon)
    }

    ## cat('Loop over chains')
##################################################
#### LOOP OVER CHAINS (WITHIN ONE CORE)
##################################################
    ## Start timer
    MCtimestart <- Sys.time()

    maxusedcomponents <- 0
    maxiterations <- 0
    ## keep count of chains having non-finite outputs
    nonfinitechains <- 0L
    ## keep count of chains stopped prematurely
    stoppedchains <- 0L
    usedmem <- max(usedmem, sum(gc()[,6])) #garbage collection
#### LOOP OVER CHAINS IN CORE
    nchainsperthiscore <- minchainspercore + (acore <= coreswithextrachain)
    ## print2user(paste0('\ncore ',acore,': ',nchainsperthiscore,'\n'), outcon)

    for (achain in 1:nchainsperthiscore) {

        chainnumber <- achain + minchainspercore * (acore - 1) +
            min(coreswithextrachain, acore - 1)
        padchainnumber <- sprintf(paste0('%0', nchar(nchains), 'i'), chainnumber)

        ## calculate how the remaining time is allotted to remaining chains
        timeleft <- (maxhours -
                         as.double(Sys.time() - timestart0, units = 'hours')) /
            (nchainsperthiscore - achain + 1)

        showsamplertimes0 <- showsamplertimes && (achain == 1)
        ## showAlphatraces0 <- showAlphatraces && (achain==1)
        niter <- min(startupMCiterations, maxMCiterations)
        ## ## Experimental: decrease number of iterations based on previous chain
        ## if(achain == 1){
        ##   niter <- startupMCiterations
        ## }else{
        ##  niter <- max(min(startupMCiterations,remainiter*2), 128)
        ##  }

        ## ## nitertot: total number of MC samples from chain start
        nitertot <- 0L
        remainiter <- +Inf
        reset <- TRUE
        allmcsamples <- NULL
        ## allmcsamplesKA <- list(Alpha = NULL, K = NULL)
        allmcsamplesKA <- NULL
        flagll <- FALSE
        flagnonfinite <- FALSE
        cat('\n###########################################################\nChain #',
            chainnumber,
            '(chain', achain, 'of', nchainsperthiscore, 'for this core)\n'
        )
        ## Read data to be used in log-likelihood
        testdata <- readRDS(file = file.path(dirname,
            paste0('___testdata_', chainnumber, '.rds')))
        cat('\nDatapoints for testing chain behaviour:\n',
            paste0('#', testdata$pointsid), '\n')

        ## will contain the MC traces of the test points
        traces <- matrix(NA_real_, nrow = 0, ncol = 1,
            dimnames = list(NULL, 'gmean'))


        ## Initial values for this chain
        ## random seed is taken care of by %doRNG%
        Cfinitemixnimble$setInits(initsfn())

        subiter <- 1L
        savedchunks <- 0L
#### WHILE-LOOP CONTINUING UNTIL STOP CRITERION MET
        while (remainiter > 0) {
            cat('\nIterations:', niter, '\n')

            ## MONTE-CARLO CALL
            ## If reporting Alpha or K traces,
            ## then save them more frequently
            ## Otherwise just save the last value
            Cmcsampler$run(
                niter = niter,
                thin = 1,
                thin2 = (if (showAlphatraces || showKtraces) {
                    max(1, round(niter / ncomponentsamples))
                } else {
                    max(2, floor(niter / 2))
                }),
                nburnin = 0,
                time = showsamplertimes0,
                reset = reset,
                resetMV = TRUE
            )

            ## iterationAsLastIndex: See sect 7.7 of Nimble manual
            mcsamples <- as.list(Cmcsampler$mvSamples,
                iterationAsLastIndex = TRUE)
            mcsamplesKA <- as.list(Cmcsampler$mvSamples2,
                iterationAsLastIndex = FALSE)

            ## 'mcsamplesKA$K' contains the component identity
            ## of each training datapoint, but we only want
            ## the number of distinct components used:
            mcsamplesKA$K <- apply(mcsamplesKA$K, 1,
                function(xx){length(unique(xx))})

            if (showAlphatraces) {
                dim(mcsamplesKA$Alpha) <- NULL # from matrix to vector
            }

            cat('\nCurrent time:',
                format(Sys.time(), '%Y-%m-%d %H:%M:%S'))
            cat('\nMCMC time',
                printtimediff(difftime(Sys.time(), MCtimestart, units = 'auto')),
                '\n')

#### Remove iterations with non-finite values
            if(any(!is.finite(unlist(mcsamples)))) {
                toRemove <- sort(unique(unlist(lapply(mcsamples, function(xx) {
                    temp <- which(!is.finite(xx), arr.ind = TRUE)
                    temp[, ncol(temp)]
                }))))

                cat('\nWARNING:', length(toRemove), 'NON-FINITE SAMPLES\n')
                ##
                flagnonfinite <- TRUE
                nonfinitechains <- TRUE
                saveRDS(mcsamples, file = file.path(dirname,
                    paste0('NONFINITEmcsamples',
                        dashnameroot, '--', padchainnumber,
                        '_', achain, '-',
                        acore, '-i', nitertot, '.rds')
                ))
                if (length(toRemove) == ncol(mcsamples$W)) {
                    cat('\n...TOO MANY NON-FINITE OUTPUTS!\n')
                    ## ## registerDoSEQ()
                    ## if(exists('cl')){ parallel::stopCluster(cl) }
                    ## stop('...TOO MANY NON-FINITE OUTPUTS. ABORTING')
                    mcsamples <- NULL
                    niter <- 0
                } else {
                    mcsamples <- mcsubset(mcsamples, -toRemove)
                    niter <- niter - length(toRemove)
                }
            }

            nitertot <- nitertot + niter

            ##
            if (showsamplertimes0) {
                samplertimes <- Cmcsampler$getTimes()
                names(samplertimes) <- sapply(
                    confnimble$getSamplers(),
                    function(x) x$target
                )
                sprefixes <- unique(sub(
                    '^([^[]+)(\\[.*\\])', '\\1',
                    names(samplertimes)
                ))
                cat('\nSampler times:\n')
                print(sort(sapply(sprefixes,
                    function(x) sum(samplertimes[grepl(x, names(samplertimes))])),
                    decreasing = TRUE
                ))
            }

            ## Check how many components were used at the last step
            ## usedcomponents <- mcsamplesKA$K[length(mcsamplesKA$K)]
            usedcomponents <- max(mcsamplesKA$K)
            cat('\nUSED COMPONENTS:', usedcomponents, 'OF', ncomponents, '\n')

#### Update traces
            ## Log-likelihood
            diagntime <- Sys.time()
            ##
            ll <- util_Pcheckpoints(
                testdata = testdata,
                learnt = mcsamples
            )

            ll <- cbind(
                exp(rowMeans(log(ll), na.rm = TRUE)), # geometric mean
                ll
            )
            colnames(ll) <- c('gmean', testdata$pointsid)

            if(plottraces){
                saveRDS(ll,
                    file = file.path(dirname,
                        paste0('____tempPtraces-',
                            padchainnumber, '-',
                            subiter, '.rds'))
                )
            }

            traces <- rbind(traces, ll[, 1, drop = FALSE])

            toRemove <- which(!is.finite(traces))
            if (length(toRemove) > 0) {
                flagll <- TRUE
                oktraces <- traces[-unique(toRemove), , drop = FALSE]
            } else {
                oktraces <- traces
            }

            ## ## version before geom.mean
            ## toRemove <- which(!is.finite(traces), arr.ind = TRUE)
            ## if (length(toRemove) > 0) {
            ##     flagll <- TRUE
            ##     oktraces <- traces[-unique(toRemove[, 1]), , drop = FALSE]
            ## } else {
            ##     oktraces <- traces
            ## }

            ## ## for debugging
            ## saveRDS(traces,
            ##     file = file.path(dirname,
            ##         paste0('___traces',
            ##             dashnameroot, '--',
            ##             padchainnumber,'-',subiter, '.rds')
            ##     ))
            ## saveRDS(oktraces,
            ##     file = file.path(dirname,
            ##         paste0('___oktraces',
            ##             dashnameroot, '--',
            ##             padchainnumber,'-',subiter, '.rds')
            ##     ))


            ## ##############################
            ## ## MCSE, ESS, THINNING, ETC ##
            ## ##############################

            N <- nrow(oktraces)

            ## quantiles to monitor
            Xlo <- quantile(x = oktraces, probs = Qlo,
                na.rm = FALSE, names = FALSE, type = 6)
            Xhi <- quantile(x = oktraces, probs = Qhi,
                na.rm = FALSE, names = FALSE, type = 6)
            ## quantile width
            width <- Xhi - Xlo

            ## CIs for lower and upper quantiles
            temp <- funMCEQ(x = oktraces, prob = c(Qlo, Qhi), Qpair = Qerror)
            wXlo <- temp[2, 1] - temp[1, 1]
            wXhi <- temp[2, 2] - temp[1, 2]

            ## Transform samples to normalized ranks, as in Vehtari et al. 2021
            essnrmean <- funESS3(qnorm(
            (rank(oktraces, na.last = NA, ties.method = 'average') - 0.5) / N
            ))

            ## We check: relative error of quantiles and ess of norm-rank-mean
            relmcse <- c(1 / sqrt(essnrmean), wXlo / width, wXhi / width)

            autothinning <- N * max(relmcse)^2

            ## Output available diagnostics
            toprint <- list(
                'rel. CI error' = max(relmcse[-1]),
                'ESS' = essnrmean,
                'needed thinning' = autothinning,
                'average' = mean(oktraces),
                'width' = width
            )

####
            for(i in names(toprint)) {
                thisdiagn <- toprint[[i]]
                cat(paste0('\n', i, ':'),
                    if(length(thisdiagn) > 1){
                        paste(signif(range(thisdiagn), 3), collapse = ' to ')
                    } else {
                        signif(thisdiagn, 3)
                    }
                )
            }

            rm(oktraces)

            cat('\nDiagnostics time',
                printtimediff(difftime(Sys.time(), diagntime, units = 'auto')),
                '\n')

#### concatenate samples with those of previous chunk, if existing
            ## add overall iteration index `MCindex` to mcsamples
            mcsamples$MCindex <- seq(to = nitertot,
                length.out = ncol(mcsamples$W))
            ## need to add 1D to behave well with mcsubset()
            dim(mcsamples$MCindex) <- ncol(mcsamples$W)

            ## concatenate samples
            if (is.null(allmcsamples)) {
                allmcsamples <- mcsamples
            } else {
                allmcsamples <- mcjoin(allmcsamples, mcsamples)
            }

            if (is.null(allmcsamplesKA)) {
                allmcsamplesKA <- mcsamplesKA
            } else {
                ## Concatenate samples of K and Alpha
                ## if (showAlphatraces || showKtraces) {
                allmcsamplesKA <- mapply(
                    function(xx, yy) {
                        c(xx, yy)
                    },
                    allmcsamplesKA, mcsamplesKA,
                    SIMPLIFY = FALSE
                )
                ## }
            }

            ## ######################################
            ## ## CHECK IF CHAIN MUST BE CONTINUED ##
            ## ######################################

            lowess <- (1 / max(relmcse))^2

            neededsamples <- max(ceiling(autothinning * (minESS + initES)),
                nsamplesperchain)

            ## if(relmcse > minrelMCSE)
            if(lowess < minESS + initES) {
                ## too small ESS
                reqiter <- neededsamples - nitertot
            } else {
                reqiter <- 0
            }

            ## respect min and max iterations chosen by user
            remainiter <- max(minMCiterations - nitertot, reqiter)
            remainiter <- min(maxMCiterations - nitertot, remainiter)

            cat('\nTotal iterations', nitertot,
                '- require further', reqiter,
                '- continue for', remainiter, '\n')

            if(as.double(Sys.time() - timestart0, units = 'hours') >= timeleft &&
                   remainiter > 0 && nitertot >= nsamplesperchain
            ) {
                cat('but stopping chain owing to lack of time\n')
                stoppedchains <- stoppedchains + 1L
                remainiter <- 0
            }

            if (remainiter > 0) { # This chain is going to continue

                ## Save cumulated mcsamples to save memory
                if(niter >= startupMCiterations){
                    savedchunks <- savedchunks + 1L
                    saveRDS(
                        mcsubset(allmcsamples, seq_len(startupMCiterations)),
                        file = file.path(dirname,
                            paste0('____tempmcsamples-',
                                padchainnumber, '-',
                                savedchunks, '.rds'))
                    )

                    if(ncol(allmcsamples$W) < startupMCiterations){
                        ## discard saved samples
                        allmcsamples <- mcsubset(allmcsamples,
                            -seq_len(startupMCiterations))
                    } else {
                        ## all samples saved
                        allmcsamples <- NULL
                    }
                }
                gc()

                ## limit number of iterations per loop, to save memory
                niter <- min(remainiter + 1L, startupMCiterations)
                subiter <- subiter + 1L
                cat('\nChain #', chainnumber, '- chunk', subiter,
                    '(chain', achain, 'of', nchainsperthiscore,
                    'for this core): increasing by', niter, '\n'
                )
            }

            ## ###########
            ## ## PLOTS ##
            ## ###########

#### Plot diagnostic traces of current chain
            if (plottraces) {
                cat('\nPlotting traces and samples.\n')

                ## Plot various info and traces
                ## colpalette <- 1:6 #seq_len(ncol(traces))
                ## names(colpalette) <- colnames(traces)[1:min(6, ncol(traces))]
                graphics.off()
                pdf(file.path(dirname,
                    paste0('___mcpartialtraces', dashnameroot, '--',
                        padchainnumber, '_', achain, '-', acore, '.pdf')
                ), height = 8.27, width = 11.69)

                ## Summary stats
                matplot(1:2, type = 'l', col = 'white',
                    main = paste0('Stats chain ', achain),
                    axes = FALSE, ann = FALSE)
                ## Legends
                ## legendpositions <- c('topleft', 'topright', 'bottomleft', 'bottomright')
                ## for (alegend in seq_along(grouplegends)) {
                ##     legend(x = legendpositions[alegend], bty = 'n', cex = 1.5,
                ##         legend = grouplegends[[alegend]]
                ##     )
                ## }
                legend(x = 'left', bty = 'n', cex = 0.75,
                    legend = c(
                        paste0('Chain ', chainnumber, ' - ',
                            achain, ' of core ', acore),
                        ##
                        paste0('Test points ',
                            paste0('#', testdata$pointsid, collapse=' ')
                        ),
                        ##
                        paste0('Iterations: ', nitertot),
                        ##
                        paste0('Used components: ', usedcomponents,
                            ' of ', ncomponents),
                        ##
                        'NOTES:',
                        if (flagnonfinite) {
                            'some non-finite MC outputs'
                        },
                        if (usedcomponents > ncomponents - 5) {
                            'too many components used'
                        },
                        if (flagll) {
                            'non-finite values in diagnostics'
                        }
                    )
                )
                ## Traces of likelihood and cond. probabilities
                par(mfrow = c(1, 1))
                for (avar in 1:ncol(traces)) {
                    flexiplot(
                        y = 10*log10(traces[is.finite(traces[, avar]), avar]),
                        type = 'l', lty = 1, col = 1,
                        main = paste0('#', colnames(traces)[avar], ': ',
                            paste(
                                names(toprint),
                                sapply(toprint, function(xx){
                                    signif(xx[avar], 3)
                                }),
                                collapse = ' | ', sep = ': '
                            )),
                        cex.main = 1.25,
                        ylab = paste0('log_F(#',
                            colnames(traces)[avar],
                            ')/dHart'),
                        xlab = 'Monte Carlo sample',
                        family = family
                        ## mar = c(NA, 6, NA, NA)
                    )
                }
                dev.off()
            }

#### Plot Alpha and component occupation, if required
            if (showAlphatraces || showKtraces) {
                cat('Plotting component and Alpha information.\n')
                pdf(file = file.path(dirname,
                    paste0('hyperparams_traces', dashnameroot, '--',
                        padchainnumber, '_', achain, '-', acore, '.pdf')),
                    height = 8.27, width = 11.69)

                if (showKtraces) {
                    cat('\nSTATS USED COMPONENTS:\n')
                    print(summary(allmcsamplesKA$K))
                    ##
                    flexiplot(y = allmcsamplesKA$K, ylab = 'used components',
                        xlab = 'iteration', ylim = c(0, ncomponents))
                    flexiplot(x = 0:ncomponents,
                        y = tabulate(allmcsamplesKA$K + 1, nbins = ncomponents + 1),
                        type = 'l', xlab = 'used components', ylab = NA,
                        ylim = c(0, NA))
                }
                if (showAlphatraces) {
                    cat('\nSTATS alpha:\n')
                    print(summary(allmcsamplesKA$Alpha, na.rm = TRUE))
                    flexiplot(y = allmcsamplesKA$Alpha,
                        ylab = bquote(alpha), xlab = 'iteration',
                        ylim = c(1, nalpha))
                    flexiplot(x = seq(minalpha, maxalpha, by = byalpha),
                        y = tabulate(allmcsamplesKA$Alpha, nbins = nalpha),
                        type = 'l', xlab = bquote(alpha), ylab = '',
                        ylim = c(0, NA))
                }
                dev.off()
            }

            ## ###############
            ## ## END PLOTS ##
            ## ###############

#### Print estimated end time
            fracchain <- achain - remainiter / (remainiter + nitertot)
            endTime <- Sys.time() + 180 +
                ( (nchainsperthiscore + (acore > coreswithextrachain) - fracchain) *
                      difftime(Sys.time(), MCtimestart) / fracchain )
            print2user(
                paste0(
                    '\rSampling. Core ', acore, ' estimated end time: ',
                    format(endTime, format='%Y-%m-%d %H:%M'),
                    '   '
                ),
                outcon
            )

            reset <- FALSE
        }
#### END WHILE-LOOP OVER CHUNKS OF ONE CHAIN


        ## #####################################
        ## ## BUILD AND SAVE CHAIN QUANTITIES ##
        ## #####################################

#### Determine which MC samples should be saved
        cat('\nKeeping ', nsamplesperchain, '\n')

        tokeep <- round(seq(from = autothinning * 2L,
            to = nitertot,
            length.out = nsamplesperchain
        ))
        if(length(tokeep) > length(unique(tokeep))){
            cat('\nWARNING: have to reduce thinning owing to time constraints\n')
            tokeep <- round(seq(from = max(1, nitertot - nsamplesperchain + 1L),
                to = nitertot,
                length.out = nsamplesperchain))
        }

        allmcsamples <- mcsubset(allmcsamples,
            which(allmcsamples$MCindex %in% tokeep)
        )

        for(chunk in rev(seq_len(savedchunks))){
            tempmcsamples <- readRDS(file = file.path(
                dirname,
                paste0('____tempmcsamples-',
                    padchainnumber, '-',
                    chunk, '.rds')
            ))
            tempmcsamples <- mcsubset(tempmcsamples,
                which(tempmcsamples$MCindex %in% tokeep)
            )

            allmcsamples <- mcjoin(tempmcsamples, allmcsamples)
        }

        ## Save thinned total chain
        saveRDS(allmcsamples,
            file = file.path(dirname,
                paste0('___mcsamples',
                    dashnameroot, '--',
                    padchainnumber, '.rds'))
        )

        saveRDS(allmcsamplesKA,
            file = file.path(dirname,
                paste0('___hyperparams_traces',
                    dashnameroot, '--',
                    padchainnumber, '.rds'))
        )

        ## put 'tokeep' in first slot if saving only nsamplesperchain
        saveRDS(traces[ , , drop = FALSE],
            file = file.path(dirname,
                paste0('___mcpartialtraces',
                    dashnameroot, '--',
                    padchainnumber, '.rds')
            ))

        ## possibly increase the count of chains with non-finite outputs
        nonfinitechains <- nonfinitechains + flagnonfinite

        cat('\nCurrent time:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'))
        cat('\nMCMC + diagnostics time',
            printtimediff(difftime(Sys.time(), MCtimestart, units = 'auto')),
            '\n')

        maxusedcomponents <- max(maxusedcomponents, usedcomponents)
        maxiterations <- max(maxiterations, nitertot)
    }
#### END LOOP OVER CHAINS (WITHIN ONE CORE)

    ##
    cat('\nCurrent time:',
        format(Sys.time(), '%Y-%m-%d %H:%M:%S'))

    cat('\nTotal time',
        printtimediff(difftime(Sys.time(), headertimestart, units = 'auto')),
        '\n')

#### Inform user that this core has finished
    ## ## timing text, for white spaces:
    ## ## "Sampling. Core 2 estimated end time: xxxx-xx-xx xx:xx   "
    print2user(paste0('\rCore ', acore,
        ' finished.                                        \n'),
        outcon)

    ## output information from a core,
    ## passed to the originally calling process
    cbind(
        maxusedcomponents = maxusedcomponents,
        maxiterations = maxiterations,
        nonfinitechains = nonfinitechains,
        stoppedchains = stoppedchains,
        headertime = headertime,
        MCtime = difftime(Sys.time(), MCtimestart, units = 'secs'),
        usedmem = max(usedmem, sum(gc()[,6]))
    )
}



#' Eliminate samples from a 'learnt' object
#'
#' @keywords internal
mcsubset <- function(learnt, subsamples) {
    lapply(learnt, function(xx) {
        do.call('[', c(
            list(xx),
            rep(TRUE, length(dim(xx)) - 1),
            list(subsamples),
            list(drop = FALSE))
        )
    })
}

#' Concatenate mcsample objects
#'
#' @keywords internal
mcjoin <- function(x, y){
    mapply(
        function(xx, yy) {
            temp <- c(xx, yy)
            dx <- dim(yy)[-length(dim(yy))]
            dim(temp) <- c(dx, length(temp) / prod(dx))
            temp
        },
        x, y,
        SIMPLIFY = FALSE
    )
}

#' Bind 3D arrays by first dimension
#'
#' @keywords internal
learnbind <- function(x, y) {
    if(is.null(x)) {
        y
    } else {
        nrx <- dim(x)[1]
        nry <- dim(y)[1]
        out <- array(data = NA, dim = c(nrx + nry, dim(x)[-1]),
            dimnames = NULL)
        out[seq_len(nrx), ,] <- x
        out[nrx + seq_len(nry), ,] <- y
        out
    }
}
## ## old, slower variant
## learnbind2 <- function(x, y) {
##     out <- c(aperm(x), aperm(y))
##     dim(out) <- c(rev(dim(x)[-1]), dim(x)[1] + dim(y)[1])
##     aperm(out)
## }


#' Cumulative sum along first dimension
#'
#' @keywords internal
rowcumsum <- function(x){
    for(i in 2:(dim(x)[1])){
        x[i,,] <- x[i,,] + x[i-1,,]
    }
    x
}

#' Inverse cumulative sum along first dimension
#'
#' @keywords internal
rowinvcumsum <- function(x){
    for(i in (dim(x)[1] - 1):1){
        x[i,,] <- x[i,,] + x[i+1,,]
    }
    x
}
