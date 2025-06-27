#### MCMC initialization function
#' @keywords internal
initsfnPrecluster <- function() {
    ## Create components centres
    ## distance function
    ## NB: all variances will be initialized to 1
    lpnorm <- function(xx){abs(xx)}
    distances <- matrix(0, nrow = constants$npoints, ncol = constants$ncomponents)
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
}

#### MCMC initialization function
#' @keywords internal
initsfnPrior<- function() {
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
}


#### MCMC initialization function
#' @keywords internal
initsfnPriormax<- function() {
}


#### MCMC initialization function
#' @keywords internal
initsfnAllcentre<- function() {
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
}



#### MCMC initialization function
#' @keywords internal
initsfnAllinone <- function() {
    initsfn <<- function(){
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
}}


#### Calculate the log-probability of a MCMC sample (produced by initsfn)
#' @keywords internal
mcmclogprob <- function(initvalues){
    with(initvalues, {
        ##
        nimble::dcat(Alpha, prob = probalpha0[1:constants$nalpha], log = TRUE) +
            ##
            nimble::ddirch(W,
                alpha = dirchalphas[1:constants$ncomponents] * alphabase^Alpha,
                log = TRUE) +
            ##
            ##
            sum(sapply(1:constants$ncomponents, function(k) {
                ##
                (if (vn$R > 0) { # continuous open domain
                    sum(sapply(1:Rn, function(v){
                        dnorm(Rmean[v, k],
                            mean = Rmean1,
                            sd = sqrt(Rvarm1),
                            log = TRUE) +
                            dinvgamma(Rrate[v, k],
                                shape = Rshapehi,
                                rate = Rvar1,
                                log = TRUE) +
                            dinvgamma(Rvar[v, k],
                                shape = Rshapelo,
                                rate = Rrate[v, k],
                                log = TRUE)
                    }), na.rm = TRUE)} else {0}) +
                    ##
                    (if (vn$C > 0) { # continuous closed domain
                        sum(sapply(1:Cn, function(v){
                            dnorm(Cmean[v, k],
                                mean = Cmean1,
                                sd = sqrt(Cvarm1),
                                log = TRUE) +
                                nimble::dinvgamma(Crate[v, k],
                                    shape = Cshapehi,
                                    rate = Cvar1,
                                    log = TRUE) +
                                nimble::dinvgamma(Cvar[v, k],
                                    shape = Cshapelo,
                                    rate = Crate[v, k],
                                    log = TRUE)
                        }), na.rm = TRUE)} else {0}) +
                    ##
                    (if (vn$D > 0) { # continuous rounded
                        sum(sapply(1:Dn, function(v){
                            dnorm(Dmean[v, k],
                                mean = Dmean1,
                                sd = sqrt(Dvarm1),
                                log = TRUE) +
                                nimble::dinvgamma(Drate[v, k],
                                    shape = Dshapehi,
                                    rate = Dvar1,
                                    log = TRUE) +
                                nimble::dinvgamma(Dvar[v, k],
                                    shape = Dshapelo,
                                    rate = Drate[v, k],
                                    log = TRUE)
                        }), na.rm = TRUE)} else {0}) +
                    ##
                    (if (vn$O > 0) { # ordinal
                        sum(sapply(1:On, function(v){
                            nimble::ddirch(Oprob[v, k, 1:Omaxn],
                                alpha = Oalpha0[v, 1:Omaxn],
                                log = TRUE)
                        }), na.rm = TRUE)} else {0}) +
                    ##
                    (if (vn$N > 0) { # ordinal
                        sum(sapply(1:Nn, function(v){
                            nimble::ddirch(Nprob[v, k, 1:Nmaxn],
                                alpha = constants$Nalpha0[v, 1:Nmaxn],
                                log = TRUE)
                        }), na.rm = TRUE)} else {0}) +
                    ##
                    (if (vn$B > 0) { # binary
                        sum(sapply(1:Bn, function(v){
                            dbeta(Bprob[v, k],
                                shape1 = Bshapelo,
                                shape2 = Bshapehi,
                                log = TRUE)
                        }), na.rm = TRUE)} else {0})
                ##
            }), na.rm = TRUE) +
            ##
            ##
            sum(sapply(1:constants$npoints, function(d){
                ##
                nimble::dcat(K[d], prob = W[1:constants$ncomponents], log = TRUE) +
                    ##
                    (if (vn$R > 0) { # continuous open domain
                        sum(sapply(1:Rn, function(v){
                            dnorm(Rdata[d, v],
                                mean = Rmean[v, K[d]],
                                sd = sqrt(Rvar[v, K[d]]),
                                log = TRUE)
                        }), na.rm = TRUE)} else {0}) +
                    ##
                    (if (vn$C > 0) { # continuous closed domain
                        sum(sapply(1:Cn, function(v){
                            nimble::dconstraint(
                                Caux[d, v], Clat[d, v] >= Cleft[d, v] &
                                                Clat[d, v] <= Cright[d, v],
                                log = TRUE) +
                                dnorm(Clat[d, v],
                                    mean = Cmean[v, K[d]],
                                    sd = sqrt(Cvar[v, K[d]]),
                                    log = TRUE)
                        }), na.rm = TRUE)} else {0}) +
                    ##
                    (if (vn$D > 0) { # continuous rounded
                        sum(sapply(1:Dn, function(v){
                            nimble::dconstraint(
                                Daux[d, v], Dlat[d, v] >= Dleft[d, v] &
                                                Dlat[d, v] < Dright[d, v],
                                log = TRUE) +
                                dnorm(Dlat[d, v],
                                    mean = Dmean[v, K[d]],
                                    sd = sqrt(Dvar[v, K[d]]),
                                    log = TRUE)
                        }), na.rm = TRUE)} else {0}) +
                    ##
                    (if (vn$O > 0) { # nominal
                        sum(sapply(1:On, function(v){
                            nimble::dcat(Odata[d, v],
                                prob = Oprob[v, K[d], 1:Omaxn],
                                log = TRUE)
                        }), na.rm = TRUE)} else {0}) +
                    ##
                    (if (vn$N > 0) { # nominal
                        sum(sapply(1:Nn, function(v){
                            nimble::dcat(Ndata[d, v],
                                prob = Nprob[v, K[d], 1:Nmaxn],
                                log = TRUE)
                        }), na.rm = TRUE)} else {0}) +
                    ##
                    (if (vn$B > 0) { # binary
                        sum(sapply(1:Bn, function(v){
                            dbinom(Bdata[d, v],
                                size = 1,
                                prob = Bprob[v, K[d]],
                                log = TRUE)
                        }), na.rm = TRUE)} else {0})
                ##
            }), na.rm = TRUE)
    })
}
