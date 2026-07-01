#' Calculate mutual information between groups of joint variates
#'
#' @description This function calculates various entropic information measures between two grops of joint variates: the mutual information, the conditional entropies, and the entropies.
#'
#' @details If \eqn{Y_1} and \eqn{Y_2} are two variates, each of which can be a joint variate such as \eqn{Y_1 = (Y_{1,1}, Y_{1,2}, \dotsc)}, and \eqn{X} a third, also possibly join, variate, then the mutual information \eqn{\mathit{MI}} between \eqn{Y_1} and \eqn{Y_2}, conditional on \eqn{X = x}, is given by
#' \deqn{\mathit{MI}(Y_1, Y_2 \vert X = x) \mathrel{:=}
#' \sum_{y_1, y_2}
#' \mathrm{Pr}(Y_1 = y_1, Y_2 = y_2 \vert X = x, \text{data})
#' \log_2\frac{
#' \mathrm{Pr}(Y_1 = y_1, Y_2 = y_2 \vert X = x, \text{data})
#' }{
#' \mathrm{Pr}(Y_1 = y_1 \vert X = x, \text{data})
#' \cdot
#' \mathrm{Pr}(Y_2 = y_2 \vert X = x, \text{data})
#' } \, \mathrm{Sh}
#' }
#' an expression which can also be written in several other equivalent ways. It is an information-theoretic measure of association that is model-free, that is, does not depend on assumptions such as linearity, gaussianity, and similar. See `vignette('mutualinfo')` for discussion and example uses, and also the "References" section.  If \eqn{Y_1, Y_2} are *jointly gaussian variates*, then there is a mathematical correspondence between their mutual information and their Pearson correlation coefficient; see output `MI.rGauss` in the "Value" section.
#'
#' The conditional entropy of \eqn{Y_1} with respect to \eqn{Y_2}, conditional on \eqn{X = x}, is given by
#' \deqn{\mathit{CondEn12}(Y_1, Y_2 \vert X = x) \mathrel{:=}
#' -\sum_{y_1, y_2}
#' \mathrm{Pr}(Y_1 = y_1 \vert Y_2 = y_2, X = x, \text{data})
#' \log_2
#' \mathrm{Pr}(Y_1 = y_1 \vert Y_2 = y_2, X = x, \text{data})
#' \cdot
#' \mathrm{Pr}(Y_2 = y_2 \vert X = x, \text{data})
#' \, \mathrm{Sh}
#' }
#'
#' The (differential) entropy of \eqn{Y_1}, conditional on \eqn{X = x}, is given by
#' \deqn{\mathit{En1}(Y_1 \vert X = x) \mathrel{:=}
#' -\sum_{y_1}
#' \mathrm{Pr}(Y_1 = y_1 \vert X = x, \text{data})
#' \log_2
#' \mathrm{Pr}(Y_1 = y_1 \vert  X = x, \text{data})
#' \, \mathrm{Sh}
#' }
#'
#' see "References" section for discussions about entropy and conditional entropy.
#'
#' The function `mutualinfo()` calculates the quantities above for the joint variates specified in the arguments `Y1names` and `Y2names`, conditional on the values of the variates specified in the data frame `X`. If `X` is omitted or `NULL`, then the posterior probabilities \eqn{\mathrm{Pr}(Y_1 | \text{data})} etc. are used. Each variate in the argument `X` can be specified either as a point-value \eqn{X = x} or as a left-open interval \eqn{X \le x} or as a right-open interval \eqn{X \ge x}, through the argument `tails`.
#'
#' The computation of these quantities is done via Monte Carlo integration, using the samples produced by the [learn()] function. The present function also output the numerical error associated with this computation.
#'
#' @param Y1names Character vector: first group of joint variates
#' @param Y2names Character vector or `NULL`: second group of joint variates
#' @param X Matrix or data.frame or `NULL`: values of some variates conditional on which we want the probabilities.
#' @param learnt Either a character with the name of a directory or full path
#'   for an 'learnt.rds' object, or such an object itself.
#' @param tails Named vector or list, or `NULL` (default). The names must match some or all of the variates in arguments `X`. For variates in this list, the probability conditional is understood in an semi-open interval sense: \eqn{X \le x} or \eqn{X \ge x}, an so on. See analogous argument in [Pr()].
#' @param n Integer or `NULL` (default): number of samples from which to approximately calculate the mutual information. Default as many as Monte Carlo samples in `learnt`.
#' @param unit Either one of 'Sh' for *shannon* (default), 'Hart' for *hartley*, 'nat' for *natural unit*, or a positive real indicating the base of the logarithms to be used.
#' @param parallel Logical or positive integer or cluster object. `TRUE` (default): use roughly half of available cores; `FALSE`: use serial computation; integer: use this many cores. It can also be a cluster object previously created with [parallel::makeCluster()]; in this case the parallel computation will use this object.
#' @param verbose Logical, default `FALSE`: give messages about parallel processing?
#'
#' @return A list consisting of the following elements:
#'
#' - `MI`, a vector of `value` and `accuracy`: the mutual information between (joint) variates `Y1names` and (joint) variates `Y2names`.
#' - `CondEn12`, `CondEn21`, vectors of `value` and `accuracy`: the conditional entropy of the first variate given the second, and vice versa.
#' - `En1`, `En2`, vectors of `value` and `accuracy`: the (differential) entropies of the first and second variates.
#' - `MI.rGauss`, a vector of `value` and `accuracy`: the absolute value of the Pearson correlation coefficient \eqn{r} of a *multivariate Gaussian distribution* having mutual information `MI`; the two are related by \eqn{\mathrm{MI} = -\ln(1 - r^2)/2}. It may provide a vague intuition for the `MI` value for people more familiar with Pearson's correlation, but should be taken with a grain of salt.
#' - `unit`, `Y1names`, `Y1names`: same as the input arguments, included for the user's convenience.
#'
#' @seealso
#' [Pr()] to calculate probabilities and their variability.
#'
#' [learn()], which generates the `learnt` objects required by `mutualinfo()`.
#'
#' @examples
#' ## Load the example `learnt` object calculated from the "penguins" dataset;
#' ## variates: 'species' and 'bill_len'
#' learnt <- learntExample
#'
#' ## mutual information between variates 'species' and 'bill_len'
#' MI <- mutualinfo(Y1names = 'species', Y2names = 'bill_len',
#'   learnt = learnt, parallel = 1)
#'
#' paste0(MI$MI, ' ', MI$unit, collapse = ' +/- ')
#'
#' ## Shannon entropy of variate 'species'
#' paste0(MI$En1, ' ', MI$unit, collapse = ' +/- ')
#'
#'
#' \donttest{
#' ## Shannon entropy of variate 'species',
#' ## conditional on a bill length of 30 mm:
#' entr <- mutualinfo(
#'   Y1names = 'species',
#'   X = data.frame(bill_len = 30),
#'   learnt = learnt, parallel = 1
#' )
#'
#' paste0(entr$En1, ' ', entr$unit, collapse = ' +/- ')
#'
#' ## the entropy is now lower; indeed a penguin with a short bill length
#' ## is most probably of the 'Adelie' species:
#' probs <- Pr(
#'   Y = data.frame(species = c('Adelie', 'Gentoo', 'Chinstrap')),
#'   X = data.frame(bill_len = 30),
#'   learnt = learnt, parallel = 1
#' )
#'
#' print(probs)
#' }
#'
#' @importFrom extraDistr rcat
#' @importFrom extraDistr rbern
#' @import parallel
#' @import stats
#' @import utils
#'
#' @export
mutualinfo <- function(
    Y1names,
    Y2names = NULL,
    X = NULL,
    learnt,
    tails = NULL,
    n = NULL,
    unit = 'Sh',
    parallel = TRUE,
    verbose = FALSE
){
#### Mutual information and conditional entropy between Y2 and Y1
#### conditional on X, data, prior information
#### are calculated by Monte Carlo integration:
#### 0. adjusted component weights are calculated for conditioning on X:
####     all probabilities below are conditional on X & data & prior
#### 1. joint samples of Y1_i, Y2_i are drawn
#### 2. probabilities p(Y1|Y2) are calculated for each sample
####    the conditional entropy is Monte-Carlo approximated by
####    H(Y1|Y2) = - sum_{i} log p(Y1_i | Y2_i)
#### 3. probabilities p(Y1) are calculated for each sample
####    the entropy is Monte-Carlo approximated by
####    H(Y1) = - sum_{i} log p(Y1_i)
#### 4. the mutual info is Monte-Carlo approximated by
####    I(Y1|Y2) = sum_{i} [log p(Y1_i | Y2_i) - log p(Y1_i)]
####           = -H(Y1|Y2) + H(Y1)
####
#### For these computations it is not necessary to transform the Y1,Y2 variates
#### from the internal Monte Carlo representation to the original one

#### Requested parallel processing
    ## NB: doesn't make sense to have more cores than chains
    closeexit <- FALSE
    if ('cluster' %in% class(parallel)){
        ## user provides a cluster object
        cl <- parallel
    } else if (isTRUE(parallel)) {
        ## user wants us to register a parallel backend
        ## and to choose number of cores
        ncores <- max(1,
            floor(parallel::detectCores() / 2))
        cl <- parallel::makeCluster(ncores)
        closeexit <- TRUE
        if(verbose){message('Registered ', capture.output(print(cl)), '.')}
    } else if (isFALSE(parallel)) {
        ## user wants us not to use parallel cores
        ncores <- 1
        cl <- parallel::makeCluster(ncores)
        closeexit <- TRUE
    } else if (is.numeric(parallel) &&
                   is.finite(parallel) && parallel >= 1) {
        ## user wants us to register 'parallel' # of cores
        ncores <- parallel
        cl <- parallel::makeCluster(ncores)
        closeexit <- TRUE
        if(verbose){message('Registered ', capture.output(print(cl)), '.')}
    } else {
        stop("Unknown value of argument 'parallel'.")
    }

    ## Close parallel connections if any were opened
    if(closeexit) {
        closecoresonexit <- function(){
            if(verbose){message('Closing connections to cores.')}
            parallel::stopCluster(cl)
            ## parallel::setDefaultCluster(NULL)
        }
        on.exit(closecoresonexit())
    }

    ## Extract Monte Carlo output & aux-metadata
    ## If learnt is a string, check if it's a folder name or file name
    if (is.character(learnt)) {
        ## Check if 'learnt' is a folder containing learnt.rds
        if (file_test('-d', learnt) &&
                file.exists(file.path(learnt, 'learnt.rds'))) {
            learnt <- readRDS(file.path(learnt, 'learnt.rds'))
        } else {
            ## Assume 'learnt' the full path of learnt.rds
            ## possibly without the file extension '.rds'
            learnt <- paste0(sub('.rds$', '', learnt), '.rds')
            if (file.exists(learnt)) {
                learnt <- readRDS(learnt)
            } else {
                stop("The argument 'learnt' must be a folder containing learnt.rds, or the path to an rds-file containing the output from 'learn()'.")
            }
        }
    }
    ## Add check to see that learnt is correct type of object?
    auxmetadata <- learnt$auxmetadata
    learnt$auxmetadata <- NULL
    learnt$auxinfo <- NULL
    ncomponents <- nrow(learnt$W)
    nmcs <- ncol(learnt$W)

    if(is.null(n) || n == 0) {
        n <- 1 * nmcs
    } else if(n < 0) {
        n <- -n * nmcs
    }

    if(n <= nmcs) {
        sseq <- sort(sample.int(nmcs, n))
    } else {
        sseq <- c(rep(x = seq_len(nmcs), times = n %/% nmcs),
            seq_len(nmcs)[sort.int(sample.int(nmcs, n %% nmcs))])
    }

    if(all(is.na(X))){X <- NULL}
    if(!is.null(X)){
        X <- as.data.frame(X)
        if (nrow(X) > 1) {
            warning('Only the first row of X is considered')
            X <- X[1, , drop = FALSE]
        }
    }
    Xv <- colnames(X)

    if(!is.null(tails)){
        tails <- as.list(tails)
        if(is.null(names(tails))) {
            stop('Missing variate names in "tails"')
        }
    }
    tailsv <- names(tails)
    tailscentre <- list('==', 0, '0', NULL)
    tailsleft <- list('<=', -1, '-1', 'left')
    tailsright <- list('>=', 1, '+1', 'right')
    tailsvalues <- c(tailscentre, tailsleft, tailsright)

    ## Consistency checks
    if (unit == 'Sh') {
        lbase <- log(2)
    } else if (unit == 'Hart') {
        lbase <- log(10)
    } else if (unit == 'nat') {
        lbase <- 1
    } else if (is.numeric(unit) && unit > 0) {
        lbase <- log(unit)
    } else {
        stop("unit must be 'Sh', 'Hart', 'nat', or a positive real")
    }

    if(!is.character(Y1names) || any(is.na(Y1names))){
        stop('Y1names must be a vector of variate names')
    }
    if(!is.null(Y2names) && (!is.character(Y2names) || any(is.na(Y2names)))){
        stop('Y2names must be NULL or a vector of variate names')
    }

    ## More consistency checks
    if(!all(Y1names %in% auxmetadata$name)) {
        stop('unknown Y1 variates\n')
    }
    if(length(unique(Y1names)) != length(Y1names)) {
        stop('duplicate Y1 variates\n')
    }
    ##
    if(!is.null(Y2names) && !all(Y2names %in% auxmetadata$name)){
        stop('unknown Y2 variates\n')
    }
    if(length(unique(Y2names)) != length(Y2names)) {
        stop('duplicate Y2 variates\n')
    }


    if (!all(Xv %in% auxmetadata$name)) {
        stop('unknown X variates\n')
    }
    if (length(unique(Xv)) != length(Xv)) {
        stop('duplicate X variates\n')
    }
    ##
    if(any(Y1names %in% Y2names)) {
        stop('overlap in Y1 and Y2 variates\n')
    }
    if(any(Y1names %in% Xv)) {
        stop('overlap in Y1 and X variates\n')
    }
    if(any(Y2names %in% Xv)) {
        stop('overlap in Y2 and X variates\n')
    }

    if (!all(tailsv %in% Xv)) {
        warning('variate ',
            paste0(tailsv[!(tailsv %in% Xv)], collapse = ' '),
            ' not among X; ignored\n')
        tails <- tails[tailsv %in% Xv]
        tailsv <- names(tails)
    }
    if (length(unique(tailsv)) != length(tailsv)) {
        stop('duplicate "tails" variates\n')
    }
    if(!all(tails %in% tailsvalues)) {
        stop('"tails" values must be ',
            paste0(tailsvalues, collapse = ' '), '\n')
    }

    ## transform 'tails' to -1, +1
    ## +1: '<=',    -1: '>='
    ## this is opposite of the argument convention because
    ## interval probabilities are calculated with `lower.tail = TRUE`
    ## eg:
    ## pnorm(x, mean, sd, lower.tail = FALSE) ==
    ##     pnorm(-x, -mean, sd, lower.tail = TRUE)
    tails[tails %in% tailscentre] <- NULL
    cleft <- tails %in% tailsleft
    cright <- tails %in% tailsright
    tails[cleft] <- +1
    tails[cright] <- -1
    tails <- unlist(tails)


#### Step 0. Adjust component weights W for conditioning on X
    if(is.null(X)){
        lW <- log(learnt$W)
    } else {
        lpargs <- util_lprobsargsyx(
            x = X,
            auxmetadata = auxmetadata,
            learnt = learnt,
            tails = tails
        )

        lW <- util_lprobsbase(
            xVs = lpargs$xVs[[1]],
            params = lpargs$params,
            logW =  log(learnt$W)
        ) # rows=components, columns=samples

    } # end definition of lW if non-null X


#### Combine Y1,Y2 into single Y for speed
    Ynames <- c(Y1names, Y2names)

#### STEP 1. Draw samples of Ynames (that is, Y1names,Y2names)

    lWnorm <- util_denorm(lW)
    Ws <- extraDistr::rcat(n = n, prob = t(
        apply(X = lWnorm[, sseq, drop = FALSE], MARGIN = 2, FUN = function(xx){
            xx <- exp(xx)
            xx / sum(xx, na.rm = TRUE)
        }, simplify = TRUE)
    ))

    Yout <- NULL
    vYout <- NULL

    ## R
    toselect <- which((auxmetadata$name %in% Ynames) &
                          (auxmetadata$mcmctype == 'R'))
    nvrt <- length(toselect)
    if(nvrt > 0){
        aux <- auxmetadata[toselect, ]
        vYout <- c(vYout, aux$name)
        totake <- cbind(rep.int(x = aux$id, times = rep(n, nvrt)), Ws, sseq)
        Yout <- c(Yout,
            rnorm(n = n * nvrt,
                mean = learnt$Rmean[totake],
                sd = learnt$Rsd[totake] )
        )
    }

    ## C
    toselect <- which((auxmetadata$name %in% Ynames) &
                          (auxmetadata$mcmctype == 'C'))
    nvrt <- length(toselect)
    if(nvrt > 0){
        aux <- auxmetadata[toselect, ]
        vYout <- c(vYout, aux$name)
        totake <- cbind(rep.int(x = aux$id, times = rep(n, nvrt)), Ws, sseq)
        Yout <- c(Yout,
            rnorm(n = n * nvrt,
                mean = learnt$Cmean[totake],
                sd = learnt$Csd[totake] )
        )
    }

    ## D
    toselect <- which((auxmetadata$name %in% Ynames) &
                          (auxmetadata$mcmctype == 'D'))
    nvrt <- length(toselect)
    if(nvrt > 0){
        aux <- auxmetadata[toselect, ]
        vYout <- c(vYout, aux$name)
        totake <- cbind(rep.int(x = aux$id, times = rep(n, nvrt)), Ws, sseq)
        Yout <- c(Yout,
            rnorm(n = n * nvrt,
                mean = learnt$Dmean[totake],
                sd = learnt$Dsd[totake] )
        )
    }

    ## O
    toselect <- which((auxmetadata$name %in% Ynames) &
                          (auxmetadata$mcmctype == 'O'))
    nvrt <- length(toselect)
    if(nvrt > 0){
        vYout <- c(vYout, auxmetadata$name[toselect])
        for(i in toselect) {
            aux <- auxmetadata[i, ]
            totake <- cbind(Ws, sseq)
            Yout <- c(Yout,
                extraDistr::rcat(n = n,
                    prob = apply(
                        X = learnt$Oprob[aux$indexpos + seq_len(aux$Nvalues), ,],
                        MARGIN = 1, FUN = `[`, totake,
                        simplify = TRUE) )
            )
        }
    }

    ## N
    toselect <- which((auxmetadata$name %in% Ynames) &
                          (auxmetadata$mcmctype == 'N'))
    nvrt <- length(toselect)
    if(nvrt > 0){
        vYout <- c(vYout, auxmetadata$name[toselect])
        for(i in toselect) {
            aux <- auxmetadata[i, ]
            totake <- cbind(Ws, sseq)
            Yout <- c(Yout,
                extraDistr::rcat(n = n,
                    prob = apply(
                        X = learnt$Nprob[aux$indexpos + seq_len(aux$Nvalues), ,],
                        MARGIN = 1, FUN = `[`, totake,
                        simplify = TRUE) )
            )
        }
    }

    ## B
    toselect <- which((auxmetadata$name %in% Ynames) &
                          (auxmetadata$mcmctype == 'B'))
    nvrt <- length(toselect)
    if(nvrt > 0){
        aux <- auxmetadata[toselect, ]
        vYout <- c(vYout, aux$name)
        totake <- cbind(rep.int(x = aux$id, times = rep(n, nvrt)), Ws, sseq)
        Yout <- c(Yout,
            extraDistr::rbern(n = n * nvrt,
                prob = learnt$Bprob[totake])
        )
    }

    dim(Yout) <- c(n, length(Ynames))
    Yout <- Yout[, match(Ynames, vYout), drop = FALSE]
    colnames(Yout) <- Ynames

    Yout <- vtransform(Yout,
        auxmetadata = auxmetadata,
        Rout = 'original',
        Cout = 'original',
        Dout = 'original',
        Oout = 'original',
        Nout = 'original',
        Bout = 'original',
        logjacobianOr = NULL)

    Y1transf <- Yout[, Y1names, drop = FALSE]
    Y2transf <- Yout[, Y2names, drop = FALSE]
    rm(Yout)
    gc()


#### STEP 2. Calculate sum_i log2_p(Y1|Y2) for all samples
    lpargs1 <- util_lprobsargsyx(
        x = Y1transf,
        auxmetadata = auxmetadata,
        learnt = learnt,
        tails = NULL
    )

    lpargs2 <- util_lprobsargsyx(
        x = Y2transf,
        auxmetadata = auxmetadata,
        learnt = learnt,
        tails = NULL
    )

    out <- do.call(rbind,
        parallel::parLapply(cl = cl,
        X = mapply(c, lpargs1$xVs, lpargs2$xVs, SIMPLIFY = FALSE),
        fun = util_lprobsmi,
        params1 = lpargs1$params,
        params2 = lpargs2$params,
        lWnorm = lWnorm,
        lW = lW
        )
        )

    ## Jacobian factors
    logjacobians1 <- rowSums(
        as.matrix(vtransform(Y1transf,
            auxmetadata = auxmetadata,
            logjacobianOr = FALSE)),
        na.rm = TRUE)

    logjacobians2 <- rowSums(
        as.matrix(vtransform(Y2transf,
            auxmetadata = auxmetadata,
            logjacobianOr = FALSE)),
        na.rm = TRUE)

    out[, c('CondEn12', 'CondEn21', 'En1', 'En2')] <-
        out[, c('CondEn12', 'CondEn21', 'En1', 'En2')] -
        c(logjacobians1, logjacobians2)

    out <- unlist(apply(
        X = rbind(
            value = colMeans(x = out, na.rm = TRUE) / lbase,
            accuracy = signif(x = apply(
                X = out, MARGIN = 2, FUN = sd, na.rm = TRUE, simplify = TRUE
            ) / (sqrt(n) * lbase), digits = 2)
        ),
        MARGIN = 2, FUN = list, simplify = TRUE), recursive = FALSE)

    if(out$MI['value'] < 0){
        out$MI['accuracy'] <- out$MI['accuracy'] + out$MI['value']
        out$MI['value'] <- 0
    }


    ## ## generally there's no MI maximum for continous variates
    ## mmax <- paste0('En', which.min(c(out$En1['value'], out$En2['value'])) )
    rgauss <- sqrt(1 - exp(-2 * out$MI['value'] * lbase))

    c(out,
        list(#MImax = out[[mmax]],
            MI.rGauss = c(
                rgauss,
                signif(x = out$MI['accuracy'] *
                           exp(-2 * out$MI['value'] * lbase) / rgauss,
                    digits = 2)
            ),
            unit = unit, Y1names = Y1names, Y2names = Y2names
        )
    )
}
