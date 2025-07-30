#' Format datapoints for testing of MCMC progress
#'
#' @param x Datapoints to be used for checking MCMC progress
#' @param auxmetadata auxmetadata object
#' @param pointsid Id of datapoints
#'
#' @keywords internal
#'
#' @return some arguments to be repeatedly used in util_Pcheckpoints
util_prepPcheckpoints <- function(
    x, auxmetadata, pointsid = NULL
) {

    auxV0a <- auxV0b <- auxV1a <- auxV1b <- auxV1c <- auxV1d <- NULL
    auxV2 <- auxVN1 <- auxVN2 <- auxVB <- NULL

    nV0 <- nV1 <- nV2 <- nVN <- nVB <- FALSE

    xV0 <- xV1 <- xV2 <- xVN <- xVB <- matrix(data = NA_real_,
        nrow = 0, ncol = nrow(x), dimnames = NULL)

###
### point probability density
###

### R-variates not in 'cumul'
    toselect <- which(auxmetadata$mcmctype == 'R')
    if(length(toselect) > 0){
        nV0 <- TRUE
        aux <- auxmetadata[toselect, ]
        auxV0a <- aux$id
        xV0 <- rbind(xV0,
            t(as.matrix(vtransform(
                x[, aux$name, drop = FALSE],
                auxmetadata = auxmetadata,
                Rout = 'normalized',
                logjacobianOr = NULL
            )))
        )
    }

### C-variates not in 'cumul' and with some non-boundary value
    toselect <- which(auxmetadata$mcmctype == 'C')
    if(length(toselect) > 0) {
        toselect <- toselect[sapply(toselect, function(i){
            any(x[,auxmetadata$name[i]] > auxmetadata$domainmin[i] &
                    x[,auxmetadata$name[i]] < auxmetadata$domainmax[i],
                na.rm = TRUE)
        })]
    }
    if(length(toselect) > 0){
        nV0 <- TRUE
        aux <- auxmetadata[toselect, ]
        auxV0b <- aux$id
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

### C-variates not in 'cumul' and with left boundary values
    toselect <- which(auxmetadata$mcmctype == 'C')
    if(length(toselect) > 0) {
        toselect <- toselect[sapply(toselect, function(i){
            any(x[,auxmetadata$name[i]] <= auxmetadata$domainmin[i],
                na.rm = TRUE)
        })]
    }
    if(length(toselect) > 0){
        nV1 <- TRUE
        aux <- auxmetadata[toselect, ]
        auxV1a <- aux$id
        xV1 <- rbind(xV1,
            t(as.matrix(vtransform(
                x[, aux$name, drop = FALSE],
                auxmetadata = auxmetadata,
                Cout = 'leftbound',
                logjacobianOr = NULL
            )))
        )
    }

### C-variates not in 'cumul' and with right boundary values
    toselect <- which(auxmetadata$mcmctype == 'C')
    if(length(toselect) > 0) {
        toselect <- toselect[sapply(toselect, function(i){
            any(x[,auxmetadata$name[i]] >= auxmetadata$domainmax[i],
                na.rm = TRUE)
        })]
    }
    if(length(toselect) > 0){
        nV1 <- TRUE
        aux <- auxmetadata[toselect, ]
        auxV1b <- aux$id
        xV1 <- rbind(xV1,
            - t(as.matrix(vtransform( # minus sign
                x[, aux$name, drop = FALSE],
                auxmetadata = auxmetadata,
                Cout = 'rightbound',
                logjacobianOr = NULL
            )))
        )
    }

### D-variates not in 'cumul' and with left boundary values
    toselect <- which(auxmetadata$mcmctype == 'D')
    if(length(toselect) > 0) {
        toselect <- toselect[sapply(toselect, function(i){
            any(x[,auxmetadata$name[i]] <= auxmetadata$domainminplushs[i],
                na.rm = TRUE)
        })]
    }
    if(length(toselect) > 0){
        nV1 <- TRUE
        aux <- auxmetadata[toselect, ]
        auxV1c <- aux$id
        xV1 <- rbind(xV1,
            t(as.matrix(vtransform(
                x[, aux$name, drop = FALSE],
                auxmetadata = auxmetadata,
                Dout = 'leftbound',
                logjacobianOr = NULL
            )))
        )
    }

### D-variates not in 'cumul' and with right boundary values
    toselect <- which(auxmetadata$mcmctype == 'D')
    if(length(toselect) > 0) {
        toselect <- toselect[sapply(toselect, function(i){
            any(x[,auxmetadata$name[i]] >= auxmetadata$domainmaxminushs[i],
                na.rm = TRUE)
        })]
    }
    if(length(toselect) > 0){
        nV1 <- TRUE
        aux <- auxmetadata[toselect, ]
        auxV1d <- aux$id
        xV1 <- rbind(xV1,
            - t(as.matrix(vtransform( # minus sign
                x[, aux$name, drop = FALSE],
                auxmetadata = auxmetadata,
                Dout = 'rightbound',
                logjacobianOr = NULL
            )))
        )
    }

###
### interval probability
###

### D-variates not in 'cumul'
    toselect <- which(auxmetadata$mcmctype == 'D')
    if(length(toselect) > 0) {
        toselect <- toselect[sapply(toselect, function(i){
            any(x[,auxmetadata$name[i]] > auxmetadata$domainminplushs[i] &
                    x[,auxmetadata$name[i]] < auxmetadata$domainmaxminushs[i],
                na.rm = TRUE)
        })]
    }
    if(length(toselect) > 0){
        nV2 <- TRUE
        aux <- auxmetadata[toselect, ]
        auxV2 <- aux$id
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
    Nshift <- 0L
### O-variates not in 'cumul'
    toselect <- which(auxmetadata$mcmctype == 'O')
    if(length(toselect) > 0){
        nVN <- TRUE
        aux <- auxmetadata[toselect, ]
        auxVN1 <- aux$id
        Nindices <- unlist(mapply(FUN = function(i, n) {i + seq_len(n)},
            aux$indexpos, aux$Nvalues,
            SIMPLIFY = FALSE))
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
    toselect <- which(auxmetadata$mcmctype == 'N')
    if(length(toselect) > 0){
        nVB <- TRUE
        aux <- auxmetadata[toselect, ]
        auxVN2 <- aux$id
        Nindices <- unlist(mapply(FUN = function(i, n) {i + seq_len(n)},
            aux$indexpos, aux$Nvalues,
            SIMPLIFY = FALSE))
        xVN <- rbind(xVN,
            t(as.matrix(vtransform(
                x[, aux$name, drop = FALSE],
                auxmetadata = auxmetadata,
                Nout = 'numeric',
                logjacobianOr = NULL
            ))) +
                Nshift + c(0, cumsum(aux$Nvalues[-1]))
        )
    }

###
### binary case
###

### B-variates
    toselect <- which(auxmetadata$mcmctype == 'B')
    if(length(toselect) > 0){
        aux <- auxmetadata[toselect, ]
        auxVB <- aux$id
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
        nV0 = nV0, nV1 = nV1, nV2 = nV2, nVN = nVN, nVB = nVB,
        xV0 = xV0, xV1 = xV1, xV2 = xV2, xVN = xVN, xVB = xVB,
        auxV0a = auxV0a, auxV0b = auxV0b,
        auxV1a = auxV1a, auxV1b = auxV1b, auxV1c = auxV1c, auxV1d = auxV1d,
        auxV2 = auxV2,
        auxVN1 = auxVN1, auxVN2 = auxVN2,
        auxVB = auxVB,
        V2steps = V2steps,
        pointsid = pointsid
    )
}




#' Calculate joint frequencies for checkpoints in learn()
#'
#' @param testdata List of objects calculated with util_prepPcheckpoints
#' @param learnt mcsamples object
#' 
#' @keywords internal
#'
#' @return The joint frequencies of Y correspoinding to the Monte Carlo samples
util_Pcheckpoints <- function(
    testdata, learnt
) {

    nsamples <- ncol(learnt$W)
    ncomponents <- nrow(learnt$W)

    with(testdata, {

###
### point probability density
###
        V0mean <- V0sd <- NULL

### R-variates not in 'cumul'
        if(length(auxV0a) > 0){
            V0mean <- learnbind(V0mean,
                learnt$Rmean)
            V0sd <- learnbind(V0sd,
                sqrt(learnt$Rvar))
        }

### C-variates not in 'cumul' and with some non-boundary value
        if(length(auxV0b) > 0) {
            V0mean <- learnbind(V0mean,
                learnt$Cmean[auxV0b, , , drop = FALSE])
            V0sd <- learnbind(V0sd,
                sqrt(learnt$Cvar[auxV0b, , , drop = FALSE]))
        }

###
### tail probability
###
        V1mean <- V1sd <- NULL

### C-variates not in 'cumul' and with left boundary values
        if(length(auxV1a) > 0){
            V1mean <- learnbind(V1mean,
                learnt$Cmean[auxV1a, , , drop = FALSE])
            V1sd <- learnbind(V1sd,
                sqrt(learnt$Cvar[auxV1a, , , drop = FALSE]))
        }

### C-variates not in 'cumul' and with right boundary values
        if(length(auxV1b) > 0){
            V1mean <- learnbind(V1mean,
                - learnt$Cmean[auxV1b, , , drop = FALSE]) # minus sign
            V1sd <- learnbind(V1sd,
                sqrt(learnt$Cvar[auxV1b, , , drop = FALSE]))
        }

### D-variates not in 'cumul' and with left boundary values
        if(length(auxV1c) > 0){
            V1mean <- learnbind(V1mean,
                learnt$Dmean[auxV1c, , , drop = FALSE])
            V1sd <- learnbind(V1sd,
                sqrt(learnt$Dvar[auxV1c, , , drop = FALSE]))
        }

### D-variates not in 'cumul' and with right boundary values
        if(length(auxV1d) > 0){
            V1mean <- learnbind(V1mean,
                - learnt$Dmean[auxV1d, , , drop = FALSE]) # minus sign
            V1sd <- learnbind(V1sd,
                sqrt(learnt$Dvar[auxV1d, , , drop = FALSE]))
        }

###
### interval probability
###

### D-variates not in 'cumul'
        if(length(auxV2) > 0){
            V2mean <- learnt$Dmean
            V2sd <- sqrt(learnt$Dvar)
        }

###
### discrete case
###
        VNprobs <- NULL

### O-variates not in 'cumul'
        if(length(auxVN1) > 0){
            VNprobs <- learnbind(VNprobs,
                learnt$Oprob)
        }
### N-variates
        if(length(auxVN2) > 0){
            VNprobs <- learnbind(VNprobs,
                learnt$Nprob)
        }

###
### binary case
###

### B-variates
        if(length(auxVB) > 0){
            VBprobs <- learnt$Bprob
        }

        do.call(cbind, foreach(
            xV0 = xV0,
            xV1 = xV1,
            xV2 = xV2,
            xVN = xVN,
            xVB = xVB
        ) %do% {
            lprobY <- util_lprobs(
                nV0 = nV0,
                V0mean = V0mean,
                V0sd = V0sd,
                xV0 = xV0,
                nV1 = nV1,
                V1mean = V1mean,
                V1sd = V1sd,
                xV1 = xV1,
                nV2 = nV2,
                V2mean = V2mean,
                V2sd = V2sd,
                V2steps = V2steps,
                xV2 = xV2,
                nVN = nVN,
                VNprobs = VNprobs,
                xVN = xVN,
                nVB = nVB,
                VBprobs = VBprobs,
                xVB = c(xVB)
            )
#### Output: rows=components, columns=samples
            lprobX <- log(learnt$W)
            colSums(exp(lprobX + lprobY)) / colSums(exp(lprobX))
        })
    })
}
