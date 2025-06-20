#' Calculate joint frequencies for checkpoints in learn()
#' @keywords internal
#' @param Y matrix or data.table: values of some already-transformed
#'   variates of which we want the joint probability; one variate per column
#' @param learnt object internal to `learn()`,
#'   containing partial Monte Carlo draws
#' @param auxmetadata object internal to `learn()`,
#'   containing processed metadata information
#'
#' @return The joint frequencies of Y correspoinding to the Monte Carlo samples
util_Pcheckpoints <- function(
    Y,
    learnt,
    auxmetadata
) {

    nsamples <- ncol(learnt$W)
    ncomponents <- nrow(learnt$W)

### Guide to indices:
    ## .i. = order in Y corresponding to appearance in vnames
    ## .t. = vnames present in Y, kept in their vnames-order

    Yv <- colnames(Y)

#### Type R
    vnames <- auxmetadata[auxmetadata$mcmctype == 'R', 'name']
    YiR <- match(vnames, Yv)
    YtR <- which(!is.na(YiR))
    YiR <- YiR[YtR]
    YnR <- length(YiR)
    if (YnR > 0) {
        learnt$Rvar <- sqrt(learnt$Rvar)
    }

#### Type C
    vnames <- auxmetadata[auxmetadata$mcmctype == 'C', 'name']
    YiC <- match(vnames, Yv)
    YtC <- which(!is.na(YiC))
    YiC <- YiC[YtC]
    YnC <- length(YiC)
    if (YnC > 0) {
        learnt$Cvar <- sqrt(learnt$Cvar)
        Clefts <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainmin']
        Crights <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainmax']
    }

#### Type D
    vnames <- auxmetadata[auxmetadata$mcmctype == 'D', 'name']
    YiD <- match(vnames, Yv)
    YtD <- which(!is.na(YiD))
    YiD <- YiD[YtD]
    YnD <- length(YiD)
    if (YnD > 0) {
        learnt$Dvar <- sqrt(learnt$Dvar)
        Dsteps <- auxmetadata[match(vnames, auxmetadata$name), 'halfstep'] /
            auxmetadata[match(vnames, auxmetadata$name), 'tscale']
        Dlefts <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainminplushs']
        Drights <- auxmetadata[match(vnames, auxmetadata$name), 'tdomainmaxminushs']
    }

#### Type O
    vnames <- auxmetadata[auxmetadata$mcmctype == 'O', 'name']
    YiO <- match(vnames, Yv)
    YtO <- which(!is.na(YiO))
    YiO <- YiO[YtO]
    YnO <- length(YiO)

#### Type N
    vnames <- auxmetadata[auxmetadata$mcmctype == 'N', 'name']
    YiN <- match(vnames, Yv)
    YtN <- which(!is.na(YiN))
    YiN <- YiN[YtN]
    YnN <- length(YiN)

#### Type B
    vnames <- auxmetadata[auxmetadata$mcmctype == 'B', 'name']
    YiB <- match(vnames, Yv)
    YtB <- which(!is.na(YiB))
    YiB <- YiB[YtB]
    YnB <- length(YiB)


    foreach(y = t(Y), .combine = `cbind`) %do% {

#### the loop is over the columns of y
#### each instance is a 1-column vector
        if (all(is.na(y))) {
            lprobY <- array(NA, dim = c(ncomponents, nsamples))
        } else {
                lprobY <- util_lprob(
                        x = y,
                        learnt = learnt,
                        nR = YnR, iR = YiR, tR = YtR,
                        nC = YnC, iC = YiC, tC = YtC,
                        Clefts = Clefts, Crights = Crights,
                        nD = YnD, iD = YiD, tD = YtD,
                        Dsteps = Dsteps, Dlefts = Dlefts, Drights = Drights,
                        nO = YnO, iO = YiO, tO = YtO,
                        nN = YnN, iN = YiN, tN = YtN,
                        nB = YnB, iB = YiB, tB = YtB
                    )
        }
#### Output: rows=components, columns=samples
        ## ## seems to lead to garbage for extreme values
        ## lprobX <- apply(log(learnt$W), 2, function(xx) {
        ##     xx - max(xx[is.finite(xx)])
        ## })
        lprobX <- log(learnt$W)
        colSums(exp(lprobX + lprobY)) / colSums(exp(lprobX))
    }
}
