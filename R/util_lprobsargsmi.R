#' Prepare arguments for util_lprobs from data
#'
#' @keywords internal
util_lprobsargsmi <- function(
    x, auxmetadata, learnt
) {
    Xv <- colnames(x)
    nX <- nrow(x)

    xV0 <- xV1 <- xV2 <- xVN <- xVB <- matrix(NA_real_, 0, nX)

###
### point probability density
###
    nV0 <- FALSE
    V0mean <- V0sd <- NULL

### R-variates not in 'tails'
    toselect <- which((auxmetadata$name %in% Xv) &
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
                Rout = 'mi',
                logjacobianOr = NULL
            )))
        )
    }

### C-variates not in 'tails' and with some non-boundary value
    toselect <- which((auxmetadata$name %in% Xv) &
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
                Cout = 'miboundisna',
                logjacobianOr = NULL
            )))
        )
    }

###
### tail probability
###
    nV1 <- FALSE
    V1mean <- V1sd <- NULL

### C-variates not in 'tails' and with left boundary values
    toselect <- which((auxmetadata$name %in% Xv) &
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
                    Cout = 'mileftbound',
                    logjacobianOr = NULL
                )))
        )
    }

### C-variates not in 'tails' and with right boundary values
    toselect <- which((auxmetadata$name %in% Xv) &
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
                    Cout = 'mirightbound',
                    logjacobianOr = NULL
                )))
        )
    }

### D-variates not in 'tails' and with left boundary values
    toselect <- which((auxmetadata$name %in% Xv) &
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
                    Dout = 'mileftbound',
                    logjacobianOr = NULL
                ))) +
                aux$halfstep / aux$tscale
        )
    }

### D-variates not in 'tails' and with right boundary values
    toselect <- which((auxmetadata$name %in% Xv) &
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
                    Dout = 'mirightbound',
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
                    Dout = 'miboundisna',
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
                    Oout = 'mi',
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
                    Nout = 'mi',
                    logjacobianOr = NULL
                ))) +
                    Nshift + c(0, cumsum(aux$Nvalues[-1]))
        )
        Nshift <- Nshift + length(Nindices)
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
                Bout = 'mi',
                logjacobianOr = NULL
            )))
        )
    }

    list(
        nV0 = nV0, xV0 = xV0, V0mean = V0mean, V0sd = V0sd,
        nV1 = nV1, xV1 = xV1, V1mean = V1mean, V1sd = V1sd,
        nV2 = nV2, xV2 = xV2, V2mean = V2mean, V2sd = V2sd,
        V2steps = V2steps,
        nVN = nVN, VNprobs = VNprobs, xVN = xVN,
        nVB = nVB, VBprobs = VBprobs, xVB = xVB
    )
}
