#' Prepare arguments for util_lprobsyx from data
#'
#' @keywords internal
util_lprobsargsyx <- function(
    x, auxmetadata, learnt, tails = NULL
) {
    Xv <- colnames(x)
    nX <- nrow(x)
    tailsv <- names(tails)

    xV0 <- xV1 <- xV2 <- xVN <- xVB <- matrix(NA_real_, 0, nX)

###
### point probability density
###
    nV0 <- FALSE
    V0mean <- V0sd <- NULL

### R-variates not in 'tails'
    toselect <- which((auxmetadata$name %in% Xv) &
                          !(auxmetadata$name %in% tailsv) &
                          (auxmetadata$mcmctype == 'R'))
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        nV0 <- TRUE
        V0mean <- learnbind(V0mean,
            learnt$Rmean[aux$id, , , drop = FALSE])
        V0sd <- learnbind(V0sd,
            sqrt(learnt$Rvar[aux$id, , , drop = FALSE]))
        xV0 <- rbind(xV0,
            t(as.matrix(vtransform(
                x[, aux$name, drop = FALSE],
                auxmetadata = auxmetadata,
                Rout = 'normalized',
                logjacobianOr = NULL
            )))
        )
    }

### C-variates not in 'tails' and with some non-boundary value
    toselect <- which((auxmetadata$name %in% Xv) &
                          !(auxmetadata$name %in% tailsv) &
                          (auxmetadata$mcmctype == 'C'))
    if(length(toselect) > 0) {
        toselect <- toselect[sapply(toselect, function(i){
            xx <- x[,auxmetadata$name[i]]
            any(xx > auxmetadata$domainmin[i] & xx < auxmetadata$domainmax[i],
                na.rm = TRUE)
        })]
    }
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        nV0 <- TRUE
        V0mean <- learnbind(V0mean,
            learnt$Cmean[aux$id, , , drop = FALSE])
        V0sd <- learnbind(V0sd,
            sqrt(learnt$Cvar[aux$id, , , drop = FALSE]))
        xV0 <- rbind(xV0,
            t(as.matrix(vtransform(
                x[, aux$name, drop = FALSE],
                auxmetadata = auxmetadata,
                Cout = 'boundisna',
                logjacobianOr = NULL
            )))
        )
    }

###
### tail probability
###
    nV1 <- FALSE
    V1mean <- V1sd <- NULL

### R-variates in 'tails'
    toselect <- which((auxmetadata$name %in% Xv) &
                          (auxmetadata$name %in% tailsv) &
                          (auxmetadata$mcmctype == 'R'))
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        nV1 <- TRUE
        V1mean <- learnbind(V1mean,
                tails[aux$name] * learnt$Rmean[aux$id, , , drop = FALSE])
        V1sd <- learnbind(V1sd,
            sqrt(learnt$Rvar[aux$id, , , drop = FALSE]))
        xV1 <- rbind(xV1,
            tails[aux$name] *
                t(as.matrix(vtransform(
                    x[, aux$name, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Rout = 'normalized',
                    logjacobianOr = NULL
                )))
        )
    }

### C-variates in 'tails'
    toselect <- which((auxmetadata$name %in% Xv) &
                          (auxmetadata$name %in% tailsv) &
                          (auxmetadata$mcmctype == 'C'))
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        nV1 <- TRUE
        V1mean <- learnbind(V1mean,
                tails[aux$name] * learnt$Cmean[aux$id, , , drop = FALSE])
        V1sd <- learnbind(V1sd,
            sqrt(learnt$Cvar[aux$id, , , drop = FALSE]))
        xV1 <- rbind(xV1,
            tails[aux$name] *
                t(as.matrix(vtransform(
                    x[, aux$name, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Cout = 'boundnormalized',
                    logjacobianOr = NULL
                )))
        )
    }

### D-variates in 'tails'
    toselect <- which((auxmetadata$name %in% Xv) &
                          (auxmetadata$name %in% tailsv) &
                          (auxmetadata$mcmctype == 'D'))
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        nV1 <- TRUE
        V1mean <- learnbind(V1mean,
                tails[aux$name] * learnt$Dmean[aux$id, , , drop = FALSE])
        V1sd <- learnbind(V1sd,
            sqrt(learnt$Dvar[aux$id, , , drop = FALSE]))
        xV1 <- rbind(xV1,
            tails[aux$name] *
                t(as.matrix(vtransform(
                    x[, aux$name, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Dout = 'boundnormalized',
                    logjacobianOr = NULL
                ))) +
                aux$halfstep / aux$tscale
        )
    }

### C-variates not in 'tails' and with left boundary values
    toselect <- which((auxmetadata$name %in% Xv) &
                          !(auxmetadata$name %in% tailsv) &
                          (auxmetadata$mcmctype == 'C'))
    if(length(toselect) > 0) {
        toselect <- toselect[sapply(toselect, function(i){
            any(x[,auxmetadata$name[i]] <= auxmetadata$domainmin[i],
                na.rm = TRUE)
        })]
    }
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        nV1 <- TRUE
        V1mean <- learnbind(V1mean,
                learnt$Cmean[aux$id, , , drop = FALSE])
        V1sd <- learnbind(V1sd,
            sqrt(learnt$Cvar[aux$id, , , drop = FALSE]))
        xV1 <- rbind(xV1,
                t(as.matrix(vtransform(
                    x[, aux$name, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Cout = 'leftbound',
                    logjacobianOr = NULL
                )))
        )
    }

### C-variates not in 'tails' and with right boundary values
    toselect <- which((auxmetadata$name %in% Xv) &
                          !(auxmetadata$name %in% tailsv) &
                          (auxmetadata$mcmctype == 'C'))
    if(length(toselect) > 0) {
        toselect <- toselect[sapply(toselect, function(i){
            any(x[,auxmetadata$name[i]] >= auxmetadata$domainmax[i],
                na.rm = TRUE)
        })]
    }
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        nV1 <- TRUE
        V1mean <- learnbind(V1mean,
                - learnt$Cmean[aux$id, , , drop = FALSE]) # minus sign
        V1sd <- learnbind(V1sd,
            sqrt(learnt$Cvar[aux$id, , , drop = FALSE]))
        xV1 <- rbind(xV1,
               - t(as.matrix(vtransform( # minus sign
                    x[, aux$name, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Cout = 'rightbound',
                    logjacobianOr = NULL
                )))
        )
    }

### D-variates not in 'tails' and with left boundary values
    toselect <- which((auxmetadata$name %in% Xv) &
                          !(auxmetadata$name %in% tailsv) &
                          (auxmetadata$mcmctype == 'D'))
    if(length(toselect) > 0) {
        toselect <- toselect[sapply(toselect, function(i){
            any(x[,auxmetadata$name[i]] <= auxmetadata$domainminplushs[i],
                na.rm = TRUE)
        })]
    }
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        nV1 <- TRUE
        V1mean <- learnbind(V1mean,
                learnt$Dmean[aux$id, , , drop = FALSE])
        V1sd <- learnbind(V1sd,
            sqrt(learnt$Dvar[aux$id, , , drop = FALSE]))
        xV1 <- rbind(xV1,
                t(as.matrix(vtransform(
                    x[, aux$name, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Dout = 'leftbound',
                    logjacobianOr = NULL
                ))) +
                aux$halfstep / aux$tscale
        )
    }

### D-variates not in 'tails' and with right boundary values
    toselect <- which((auxmetadata$name %in% Xv) &
                          !(auxmetadata$name %in% tailsv) &
                          (auxmetadata$mcmctype == 'D'))
    if(length(toselect) > 0) {
        toselect <- toselect[sapply(toselect, function(i){
            any(x[,auxmetadata$name[i]] >= auxmetadata$domainmaxminushs[i],
                na.rm = TRUE)
        })]
    }
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        nV1 <- TRUE
        V1mean <- learnbind(V1mean,
                - learnt$Dmean[aux$id, , , drop = FALSE]) # minus sign
        V1sd <- learnbind(V1sd,
            sqrt(learnt$Dvar[aux$id, , , drop = FALSE]))
        xV1 <- rbind(xV1,
               - t(as.matrix(vtransform( # minus sign
                    x[, aux$name, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Dout = 'rightbound',
                    logjacobianOr = NULL
                ))) +
                aux$halfstep / aux$tscale
        )
    }

###
### interval probability
###
    nV2 <- FALSE
    V2mean <- V2sd <- V2steps <- NULL

### D-variates not in 'tails'
    toselect <- which((auxmetadata$name %in% Xv) &
                          !(auxmetadata$name %in% tailsv) &
                          (auxmetadata$mcmctype == 'D'))
    if(length(toselect) > 0) {
        toselect <- toselect[sapply(toselect, function(i){
            xx <- x[,auxmetadata$name[i]]
            any(xx > auxmetadata$domainminplushs[i] &
                    xx < auxmetadata$domainmaxminushs[i],
                na.rm = TRUE)
        })]
    }
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        nV2 <- TRUE
        V2mean <- learnbind(V2mean,
                learnt$Dmean[aux$id, , , drop = FALSE])
        V2sd <- learnbind(V2sd,
            sqrt(learnt$Dvar[aux$id, , , drop = FALSE]))
        V2steps <- aux$halfstep / aux$tscale
        xV2 <- rbind(xV2,
                t(as.matrix(vtransform(
                    x[, aux$name, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Dout = 'boundisna',
                    logjacobianOr = NULL
                )))
        )
    }

###
### discrete case
###
    nVN <- FALSE
    VNprobs <- NULL
    Nshift <- 0L
### O-variates not in 'tails'
    toselect <- which((auxmetadata$name %in% Xv) &
                          !(auxmetadata$name %in% tailsv) &
                          (auxmetadata$mcmctype == 'O'))
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        nVN <- TRUE
        ## Nindices <- unlist(lapply(seq_len(nrow(aux)), function(i) {
        ##     aux$indexpos[i] + seq_len(aux$Nvalues[i])
        ## }))
        Nindices <- unlist(mapply(FUN = function(i, n) {i + seq_len(n)},
            aux$indexpos, aux$Nvalues,
            SIMPLIFY = FALSE))
        VNprobs <- learnbind(VNprobs,
            learnt$Oprob[Nindices, , , drop = FALSE])
        xVN <- rbind(xVN,
                t(as.matrix(vtransform(
                    x[, aux$name, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Oout = 'numeric',
                    logjacobianOr = NULL
                ))) +
                    Nshift + c(0, cumsum(aux$Nvalues[-1]))
        )
        Nshift <- Nshift + length(Nindices)
    }
### N-variates
    toselect <- which((auxmetadata$name %in% Xv) &
                          (auxmetadata$mcmctype == 'N'))
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        nVN <- TRUE
        ## Nindices <- unlist(lapply(seq_len(nrow(aux)), function(i) {
        ##     aux$indexpos[i] + seq_len(aux$Nvalues[i])
        ## }))
        Nindices <- unlist(mapply(FUN = function(i, n) {i + seq_len(n)},
            aux$indexpos, aux$Nvalues,
            SIMPLIFY = FALSE))
        VNprobs <- learnbind(VNprobs,
            learnt$Nprob[Nindices, , , drop = FALSE])
        xVN <- rbind(xVN,
                t(as.matrix(vtransform(
                    x[, aux$name, drop = FALSE],
                    auxmetadata = auxmetadata,
                    Nout = 'numeric',
                    logjacobianOr = NULL
                ))) +
                    Nshift + c(0, cumsum(aux$Nvalues[-1]))
        )
        Nshift <- Nshift + length(Nindices)
    }

### O-variates in 'tails'
    toselect <- which((auxmetadata$name %in% Xv) &
                          (auxmetadata$name %in% tailsv) &
                          (auxmetadata$mcmctype == 'O'))
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        nVN <- TRUE
        for(i in seq_len(nrow(aux))) {
            VNprobs <- learnbind(VNprobs,
                if(tails[aux$name[i]] < 0) {
                    rowcumsum(learnt$Oprob[
                        aux$indexpos[i] + seq_len(aux$Nvalues[i]), , , drop = FALSE
                    ])
                } else {
                    rowinvcumsum(learnt$Oprob[
                        aux$indexpos[i] + seq_len(aux$Nvalues[i]), , , drop = FALSE
                    ])
                }
            )
        }
        xVN <- rbind(xVN,
            t(as.matrix(vtransform(
                x[, aux$name, drop = FALSE],
                auxmetadata = auxmetadata,
                Oout = 'numeric',
                logjacobianOr = NULL
            ))) +
                Nshift + c(0, cumsum(aux$Nvalues[-1]))
        )
    }

###
### binary case
###
    nVB <- FALSE
    VBprobs <- NULL
### B-variates
    toselect <- which((auxmetadata$name %in% Xv) &
                          (auxmetadata$mcmctype == 'B'))
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        nVB <- TRUE
        VBprobs <- learnbind(VBprobs,
            learnt$Bprob[aux$id, , , drop = FALSE])
        xVB <- rbind(xVB,
            t(as.matrix(vtransform(
                x[, aux$name, drop = FALSE],
                auxmetadata = auxmetadata,
                Bout = 'numeric',
                logjacobianOr = NULL
            )))
        )
    }

    list(
        params = list(
            nV0 = nV0, V0mean = V0mean, V0sd = V0sd,
            nV1 = nV1, V1mean = V1mean, V1sd = V1sd,
            nV2 = nV2, V2mean = V2mean, V2sd = V2sd,
            V2steps = V2steps,
            nVN = nVN, VNprobs = VNprobs,
            nVB = nVB, VBprobs = VBprobs
        ),
        xVs = lapply(seq_len(nX), function(i){
            list(
                ii = i,
                xV0 = xV0[,i],
                xV1 = xV1[,i],
                xV2 = xV2[,i],
                xVN = xVN[,i],
                xVB = xVB[,i]
            )})
    )
}



#' Calculate collection of log-probabilities for different components and samples
#' @return Matrix with as many rows as components and as many cols as samples
#' @keywords internal
util_lprobsbase <- function(
    xVs, params, logW,
    temporarydir = NULL, lab = ''
) {
    with(c(xVs, params), {
    out <- logW
    ## point probability density
    if(nV0) {
        out <- out + colSums(
            x = dnorm(x = xV0, mean = V0mean, sd = V0sd, log = TRUE),
            na.rm = TRUE, dims = 1)
    }
    ## tail probability
    if(nV1) {
        out <- out + colSums(
            x = pnorm(q = xV1, mean = V1mean, sd = V1sd,
                log.p = TRUE, lower.tail = TRUE),
            na.rm = TRUE, dims = 1)
    }
    ## interval probability
    if(nV2) {
        pright <- pnorm(q = xV2 + V2steps, mean = V2mean, sd = V2sd,
            log.p = TRUE, lower.tail = TRUE)
        ##
        out <- out + colSums(
            x = pright + log(-expm1(
                pnorm(q = xV2 - V2steps, mean = V2mean, sd = V2sd,
                    log.p = TRUE, lower.tail = TRUE) - pright
            )),
            na.rm = TRUE, dims = 1)
        ## ## this alternate form leads to infinities in some cases
        ## pleft <- pnorm(q = xV2 - V2steps, mean = V2mean, sd = V2sd,
        ##     log.p = TRUE, lower.tail = TRUE)
        ## ##
        ## out <- out + colSums(
        ##     x = pleft + log(expm1(
        ##         pnorm(q = xV2 + V2steps, mean = V2mean, sd = V2sd,
        ##             log.p = TRUE, lower.tail = TRUE) - pleft
        ##     )),
        ##     na.rm = TRUE, dims = 1)
    }
    ##
    if(nVN) {
        out <- out + colSums(
            x = log(VNprobs[xVN, , , drop = FALSE]),
            na.rm = TRUE, dims = 1)
    }
    ##
    if(nVB) {
        ## VBprob is the probability that x = 1
        out <- out + colSums(
            x = log(1 - xVB - VBprobs + 2 * xVB * VBprobs),
            na.rm = TRUE, dims = 1)
    }

    if(is.null(temporarydir)){
        out
    }else{
        saveRDS(out,
            file.path(temporarydir,
                paste0(lab, ii, '__.rds'))
        )
    }
    })
}

#' Calculate probabilities, quantiles, etc, for all Y and X combinations
#'
#' @keywords internal
util_combineYX <- function(
    iyx,
    temporarydir, usememory = TRUE,
    doquantiles, quantiles,
    dosamples, nsamples,
    Qerror
) {

    if(usememory) {
        lprobX <- readRDS(file.path(temporarydir,
            paste0('__X', iyx['jx'], '__.rds')
        ))
        lprobY <- readRDS(file.path(temporarydir,
            paste0('__Y', iyx['jy'], '__.rds')
        ))
    }

    FF <- colSums(x = exp(lprobX + lprobY), na.rm = TRUE) /
        colSums(x = exp(lprobX), na.rm = TRUE)

    list(
        values = mean(x = FF, na.rm = TRUE),
        ##
        quantiles = if(doquantiles) {
            quantile(x = FF, probs = quantiles, type = 6,
                na.rm = TRUE, names = FALSE)
        },
        ##
        samples = if(dosamples) {
            FF <- FF[!is.na(FF)]
            FF[round(seq(1, length(FF), length.out = nsamples))]
        },
        ##
        values.MCaccuracy = funMCSELD(x = FF),
        ##
        quantiles.MCaccuracy = if(doquantiles) {
            temp <- funMCEQ(x = FF, prob = quantiles, Qpair = Qerror)
            (temp[2, ] - temp[1, ]) / 2
        }
        ##
        ## error = sd(FF, na.rm = TRUE)/sqrt(nmcsamples)
    )
}





## #' Calculate collection of log-probabilities for different components and samples
## #' @return Matrix with as many rows as components and as many cols as samples
## #' @keywords internal
## util_lprobssave <- function(xVs, params, logW, temporarydir, lab) {
## 
##     out <- util_lprobsbase(xVs = xVs, params = params, logW = logW)
## 
##     saveRDS(out,
##         file.path(temporarydir,
##             paste0(lab, xVs$ii, '__.rds'))
##     )
## }

#' Utility function to avoid finite-precision accuracys
#' 
#' @keywords internal
util_denorm <- function(lprob) {
    apply(X = lprob, MARGIN = 2, FUN = function(xx) {
        xx - max(xx[is.finite(xx)])
    }, simplify = TRUE)
}


#' Calculate pairs of log-probabilities for mutualinfo()
#' 
#' @keywords internal
util_lprobsmi <- function(xVs, params1, params2, lWnorm, lW) {

    lprobY1 <- util_lprobsbase(xVs = xVs[1:6], params = params1, logW = 0)
    lprobY2 <- util_lprobsbase(xVs = xVs[7:12], params = params2, logW = 0)

    celWnorm <- colSums(exp(lWnorm))

### Construct probabilities from lprobY1, lprobY2
    lpY1and2 <- log(mean(
        colSums(exp(lprobY1 + lprobY2 + lWnorm)) / celWnorm,
        na.rm = TRUE))

    lpY1 <- log(mean(
        colSums(exp(lprobY1 + lWnorm)) / celWnorm,
        na.rm = TRUE))

    lpY2 <- log(mean(
        colSums(exp(lprobY2 + lWnorm)) / celWnorm,
        na.rm = TRUE))

    lprobnorm <- util_denorm(lprobY2 + lW)
    lpY1given2 <- log(mean(
        colSums(exp(lprobY1 + lprobnorm)) / colSums(exp(lprobnorm)),
        na.rm = TRUE))

    lprobnorm <- util_denorm(lprobY1 + lW)
    lpY2given1 <- log(mean(
        colSums(exp(lprobY2 + lprobnorm)) / colSums(exp(lprobnorm)),
        na.rm = TRUE))

    mi <- lpY1and2 - lpY1 - lpY2
    c(
        MI = mi,
        CondEn12 = -lpY1given2,
        CondEn21 = -lpY2given1,
        En1 = -lpY1,
        En2 = -lpY2
        ## MIalt = (mi + lpY1given2 - lpY1 + lpY2given1 - lpY2) / 3,
    )
}



#' Calculate quantiles for continuous Y by bisection
#'
#' @keywords internal
util_qYXcont <- function(
    iyx,
    params1, params2,
    auxmetadata,
    temporarydir, usememory = TRUE,
    doquantiles, quantiles,
    dosamples, nsamples,
    Qerror,
    eps = .Machine$double.eps * 3
) {
    pY <- iyx['pY']

    if(usememory) {
        lprobX <- readRDS(file.path(temporarydir,
            paste0('__X', iyx['jx'], '__.rds')
        ))
    }
    sumlpX <- colSums(exp(lprobX), na.rm = TRUE)

    maxsamples <- ncol(params1)

    Yvals <- .Machine$double.xmax * c(-0.125, 0.125)

#### Calculate quantile for the posterior probability distribution

    values <- (Yvals[1] + Yvals[2]) / 2

    FF <- mean(colSums(exp(
        lprobX + pnorm(q = values,
            mean = params1, sd = params2,
            lower.tail = TRUE, log.p = TRUE)
    ), na.rm = TRUE) / sumlpX) - pY

    while(abs(FF) > eps && Yvals[2] - Yvals[1] > eps){
        Yvals[(FF > 0) + 1L] <- values
        values <- (Yvals[1] + Yvals[2]) / 2
        FF <- mean(colSums(exp(
            lprobX + pnorm(q = values,
                mean = params1, sd = params2,
                lower.tail = TRUE, log.p = TRUE)
        ), na.rm = TRUE) / sumlpX) - pY
    }
    values <- unname(unlist(vtransform(values,
        auxmetadata = auxmetadata,
        Rout = 'original',
        Cout = 'original',
        Dout = 'original',
        Oout = 'original',
        Nout = 'original',
        Bout = 'original',
        variates <- auxmetadata$name,
        logjacobianOr = NULL)))


#### Calculate quantile for the frequency samples
    if(doquantiles){
        selsamples <- TRUE
    } else if(dosamples) {
        nmaxsamples <- nsamples
        selsamples <- round(seq(1, ncol(params1), length.out = nsamples))
    }

    if(doquantiles || dosamples) {
        params1 <- t(params1[, selsamples])
        params2 <- t(params2[, selsamples])
        lprobX <- t(lprobX[, selsamples])

        Yvals <- rep(.Machine$double.xmax * c(-0.125, 0.125), each = nmaxsamples)
        dim(Yvals) <- c(nmaxsamples, 2)
        sampleseq <- 1:nmaxsamples

        samples <- (Yvals[, 1] + Yvals[, 2]) / 2
        FF <- rowSums(exp(
            lprobX + pnorm(q = samples,
                mean = params1, sd = params2,
                lower.tail = TRUE, log.p = TRUE)
        ), na.rm = TRUE) / sumlpX - pY

        tocheck <- abs(FF) > eps
        while(any(tocheck)) {
            choose <- c(sampleseq[tocheck], (FF[tocheck] > 0) + 1L)
            dim(choose) <- c(sum(tocheck), 2)
            Yvals[choose] <- samples[tocheck]
            samples[tocheck] <- (Yvals[tocheck, 1] + Yvals[tocheck, 2]) / 2
            FF <- rowSums(exp(
                lprobX + pnorm(q = samples,
                    mean = params1, sd = params2,
                    lower.tail = TRUE, log.p = TRUE)
            ), na.rm = TRUE) / sumlpX - pY
        tocheck <- abs(FF) > eps & Yvals[, 2] - Yvals[, 1] > eps
        }
        samples <- unname(unlist(vtransform(samples,
            auxmetadata = auxmetadata,
            Rout = 'original',
            Cout = 'original',
            Dout = 'original',
            Oout = 'original',
            Nout = 'original',
            Bout = 'original',
            variates <- auxmetadata$name,
            logjacobianOr = NULL)))
    }

    list(
        values = values,
        ##
        quantiles = if(doquantiles) {
            quantile(x = samples, probs = quantiles, type = 6,
                na.rm = TRUE, names = FALSE)
        },
        ##
        samples = if(dosamples) {
            samples[round(seq(1, length(samples), length.out = nsamples))]
        }
        ## values.MCaccuracy
        ## quantiles.MCaccuracy
)
}



#' Calculate quantiles for discrete Y by bisection
#'
#' @keywords internal
util_qYXdiscr <- function(
    iyx,
    params1, params2,
    auxmetadata,
    temporarydir, usememory = TRUE,
    doquantiles, quantiles,
    dosamples, nsamples,
    Qerror,
    eps = NULL
) {
    pY <- iyx['pY']

    if(usememory) {
        lprobX <- readRDS(file.path(temporarydir,
            paste0('__X', iyx['jx'], '__.rds')
        ))
    }
    sumlpX <- colSums(exp(lprobX), na.rm = TRUE)

    nmaxsamples <- dim(params1)[3]
    Nvalues <- auxmetadata$Nvalues

#### Calculate quantile for the posterior probability distribution

    values <- 1L

    FF <- mean(colSums(exp(
        lprobX + params1[values, ,]
    ), na.rm = TRUE) / sumlpX)

    while(FF < pY && values <= Nvalues){
        values <- values + 1L
        FF <- mean(colSums(exp(
            lprobX + params1[values, ,]
        ), na.rm = TRUE) / sumlpX)
    }
    values <- unname(unlist(vtransform(values,
        auxmetadata = auxmetadata,
        Rout = 'original',
        Cout = 'original',
        Dout = 'original',
        Oout = 'original',
        Nout = 'original',
        Bout = 'original',
        variates <- auxmetadata$name,
        logjacobianOr = NULL)))

#### Calculate quantile for the frequency samples
    if(doquantiles){
        selsamples <- TRUE
    } else if(dosamples) {
        nmaxsamples <- nsamples
        selsamples <- round(seq(1, ncol(params1), length.out = nsamples))
    }

    if(doquantiles || dosamples) {
        params1 <- aperm(a = params1[, , selsamples],
            perm = c(1, 3, 2), resize = TRUE)
        lprobX <- t(lprobX[, selsamples])

        samples <- rep(1L, nmaxsamples)

        i <- 1L

        FF <- rowSums(exp(
            lprobX + params1[i, ,]
        ), na.rm = TRUE) / sumlpX

        tocheck <- FF < pY
        while(any(tocheck)) {
            i <-  i + 1L
            samples[tocheck] <- i
            FF <- rowSums(exp(
            lprobX + params1[i, ,]
            ), na.rm = TRUE) / sumlpX
            tocheck <- FF < pY
        }
        samples <- unname(unlist(vtransform(samples,
            auxmetadata = auxmetadata,
            Rout = 'original',
            Cout = 'original',
            Dout = 'original',
            Oout = 'original',
            Nout = 'original',
            Bout = 'original',
            variates <- auxmetadata$name,
            logjacobianOr = NULL)))
    }

    list(
        values = values,
        ##
        quantiles = if(doquantiles) {
            quantile(x = samples, probs = quantiles, type = 6,
                na.rm = TRUE, names = FALSE)
        },
        ##
        samples = if(dosamples) {
            samples[round(seq(1, length(samples), length.out = nsamples))]
        }
        ## values.MCaccuracy
        ## quantiles.MCaccuracy
)
}
