#' Build preliminary metadata flie
#'
#' @param data data.frame object
#' @param metadata data.frame object
#'
#' @return an auxmetadata data.frame object
buildauxmetadata <- function(data, metadata) {

    ## In the internal, rescaled representation,
    ## with the SD of the means equal to 3 and
    ## the rate of the variances equal to 1,
    ## the median of Q3 of the possible F distribs is 1.95
    ## and its mean is 2.55, while
    ## the median of the IQR is 3.35
    ## and its mean is 5.10.
    ## So we could set the scale to
    ## ([(Q3-Q2)+(Q2-Q1)]/2) /2
    ## = (Q3-Q1)/4
    ## this way the transformed Q3 is at approx 2
    ## the factor "4" is instead "1" at the moment.
    ## This need to be studied some more
    iqrfactor <- 2
                                        #  iqrfactorLog <- 2 * 6
                                        #  iqrfactorQ <- 2 * 0.3

    idR <- idC <- idD <- idL <- idB <- idO <- idN <- 1L

    auxmetadata <- data.frame()

    for (name in metadata$name) {
        if(!is.null(data)) {
            x <- data[[name]]
            x <- x[!is.na(x)]
        }
        minfo <- as.list(metadata[metadata$name == name, ])
        ## make sure 'type' is lowercase
        minfo$type <- tolower(minfo$type)
        ordinal <- NA
        transf <- 'identity' # temporary
        datavalues <- minfo[grep('^V[0-9]+$', names(minfo))]
        if (minfo$type == 'binary') {
            mcmctype <- 'B'
            id <- idB
            idB <- idB + 1L
            Nvalues <- minfo$Nvalues
            step <- minfo$rounding / 2
            domainmin <- 0
            domainmax <- 1
            censormin <- -Inf
            censormax <- +Inf
            tcensormin <- -Inf
            tcensormax <- +Inf
            tlocation <- 0
            tscale <- 1
            plotmin <- NA
            plotmax <- NA
            ## mctest1 <- 1
            ## mctest2 <- 2
            ## mctest3 <- 2
            ## if(!is.null(data)) {
            ##   mctest1 <- match(names(which.min(table(x))), datavalues)
            ##   mctest2 <- match(names(which.max(table(x))), datavalues)
            ##   mctest3 <- match(names(which.min(table(x))), datavalues)
            ## } else {
            ##   mctest1 <- 1
            ##   mctest2 <- 2
            ##   mctest3 <- 2
            ## }
        } else if (minfo$type == 'ordinal') {
            ## nominal variate
            mcmctype <- 'O'
            id <- idO
            idO <- idO + 1L
            Nvalues <- minfo$Nvalues
            step <- 0.5
            domainmin <- 1 # Nimble index categorical from 1
            domainmax <- Nvalues
            censormin <- -Inf
            censormax <- +Inf
            tcensormin <- -Inf
            tcensormax <- +Inf
            tlocation <- 0
            tscale <- 1
            plotmin <- NA
            plotmax <- NA
            ## if(!is.null(data)) {
            ##   mctest1 <- match(names(which.min(table(x))), datavalues)
            ##   mctest2 <- match(names(which.max(table(x))), datavalues)
            ##   mctest3 <- match(names(which.min(table(x))), datavalues)
            ## } else {
            ## mctest1 <- 1
            ## mctest2 <- round(minfo$Nvalues/2)
            ## mctest3 <- minfo$Nvalues
            ## }
            ## mctest1 <- 1
            ## mctest2 <- round(Nvalues/2)
            ## mctest3 <- Nvalues
        } else if (minfo$type == 'nominal') {
            ## nominal variate
            mcmctype <- 'N'
            id <- idN
            idN <- idN + 1L
            Nvalues <- minfo$Nvalues
            step <- 0.5
            domainmin <- 1 # Nimble index categorical from 1
            domainmax <- Nvalues
            censormin <- -Inf
            censormax <- +Inf
            tcensormin <- -Inf
            tcensormax <- +Inf
            tlocation <- 0
            tscale <- 1
            plotmin <- NA
            plotmax <- NA
            ## if(!is.null(data)) {
            ##   mctest1 <- match(names(which.min(table(x))), datavalues)
            ##   mctest2 <- match(names(which.max(table(x))), datavalues)
            ##   mctest3 <- match(names(which.min(table(x))), datavalues)
            ## } else {
            ## mctest1 <- 1
            ## mctest2 <- round(minfo$Nvalues/2)
            ## mctest3 <- minfo$Nvalues
            ## }
            ## mctest1 <- 1
            ## mctest2 <- round(Nvalues/2)
            ## mctest3 <- Nvalues
        } else if (minfo$type == 'latent') {
            ## old treatment of ordinal variate
            ## to be deleted in the future, if not used
            mcmctype <- 'L'
            id <- idL
            idL <- idL + 1L
            transf <- 'Q'
            ordinal <- TRUE
            Nvalues <- minfo$Nvalues
            step <- 0.5
            domainmin <- minfo$domainmin
            domainmax <- minfo$domainmax
            censormin <- -Inf
            censormax <- +Inf
            tcensormin <- -Inf
            tcensormax <- +Inf
            ## datavalues <- as.vector(minfo[paste0('V',1:Nvalues)], mode='character')
            olocation <- (Nvalues * domainmin - domainmax) / (Nvalues - 1)
            oscale <- (domainmax - domainmin) / (Nvalues - 1)
            if(!is.null(data)) {
                tlocation <- median(util_Q(
                    round((x - olocation) / oscale) / Nvalues
                ))
                tscale <- mad(util_Q(
                    round((x - olocation) / oscale) / Nvalues
                )) / iqrfactor
            } else {
                tlocation <- 0
                tscale <- 1 / iqrfactor
            }
            ## tlocation <- Qf(round((minfo$centralvalue - olocation) / oscale) / Nvalues)
            ## tscale <- abs(
            ##   Qf(round((minfo$highvalue - olocation) / oscale) / Nvalues) -
            ##   Qf(round((minfo$lowvalue - olocation) / oscale) / Nvalues)
            ## ) / iqrfactorQ
            plotmin <- (if(is.finite(minfo$plotmin)){minfo$plotmin}else{minfo$domainmin})
            plotmax <- (if(is.finite(minfo$plotmax)){minfo$plotmax}else{minfo$domainmax})
            ## Q1 <- mctest1 <- quantile(x, probs = 0.25, type = 6)
            ## Q2 <- mctest2 <- quantile(x, probs = 0.5, type = 6)
            ## Q3 <- mctest3 <- quantile(x, probs = 0.75, type = 6)
            ## mctest1 <- minfo$lowvalue
            ## mctest2 <- minfo$centralvalue
            ## mctest3 <- minfo$highvalue
            ## if (mctest1 == mctest3) {
            ##   mctest1 <- if (sum(x < Q1) > 0) {
            ##     max(x[x < Q1])
            ##   } else {
            ##     max(x[x <= Q1])
            ##   }
            ##   mctest3 <- if (sum(x > Q3) > 0) {
            ##     min(x[x > Q3])
            ##   } else {
            ##     min(x[x >= Q3])
            ##   }
            ## }
        } else if (minfo$type == 'continuous') { # continuous variate (R,C,D)
            Nvalues <- +Inf

            ## Rounded variates
            if (is.null(minfo$rounding) || is.na(minfo$rounding)) {
                step <- 0
            } else {
                step <- minfo$rounding / 2
            }
            ## If the variate is rounded,
            ## no latent-variable representation is needed anyway
            ## if the datapoints are distinct in higher dimension
            if(step > 0 &&
               (is.null(data) || nrow(unique(data))/nrow(data) > 0.5)) {
                step <- 0
            }

            domainmin <- minfo$domainmin
            domainmax <- minfo$domainmax
            ## censormin <- max(domainmin, minfo$minincluded, na.rm=TRUE)
            ## censormax <- min(domainmax, minfo$maxincluded, na.rm=TRUE)
            ## closeddomain <- (censormin > domainmin) || (censormax < domainmax)
            plotmin <- minfo$plotmin
            plotmax <- minfo$plotmax
            ## Q1 <- mctest1 <- quantile(x, probs = 0.25, type = 6)
            ## Q2 <- mctest2 <- quantile(x, probs = 0.5, type = 6)
            ## Q3 <- mctest3 <- quantile(x, probs = 0.75, type = 6)
            ## mctest1 <- minfo$lowvalue
            ## mctest2 <- minfo$centralvalue
            ## mctest3 <- minfo$highvalue

#### variate transformation to real-line domain
            if (is.finite(minfo$domainmin) && is.finite(minfo$domainmax) &&
                !minfo$minincluded && !minfo$maxincluded) {
                ## doubly-bounded open domain
                transf <- 'Q'
                if(!is.null(data)) {
                    tlocation <- median(util_Q(0.5 +
                                               (x - (domainmax + domainmin)/2) /
                                               (domainmax - domainmin)
                    ))
                    tscale <- mad(util_Q(0.5 +
                                         (x - (domainmax + domainmin)/2) /
                                         (domainmax - domainmin)
                    )) / iqrfactor
                } else {
                    tlocation <- 0
                    tscale <- 1 / iqrfactor
                }
                ## tlocation <- Qf(0.5 +
                ##                 (minfo$centralvalue - (domainmax + domainmin)/2) /
                ##                 (domainmax - domainmin))
                ## tscale <- abs(
                ##   Qf(0.5 +
                ##      (minfo$highvalue - (domainmax + domainmin)/2) /
                ##      (domainmax - domainmin)) -
                ##   Qf(0.5 +
                ##      (minfo$lowvalue - (domainmax + domainmin)/2) /
                ##      (domainmax - domainmin))
                ## ) / iqrfactorQ
                closeddomain <- FALSE
                censormin <- domainmin
                censormax <- domainmax
                tcensormin <- (util_Q(0.5 +
                                      (censormin - (domainmax + domainmin)/2) /
                                      (domainmax - domainmin)
                ) - tlocation) / tscale
                tcensormax <- (util_Q(0.5 +
                                      (censormax - (domainmax + domainmin)/2) /
                                      (domainmax - domainmin)
                ) - tlocation) / tscale
            } else if (is.finite(minfo$domainmin) && is.finite(minfo$domainmax) &&
                       minfo$minincluded && !minfo$maxincluded) {
                ## doubly-bounded left-closed domain
                transf <- 'logminus'
                if(!is.null(data)) {
                    tlocation <- median(log(domainmax - x))
                    tscale <- mad(log(domainmax - x)) / iqrfactor
                } else {
                    tlocation <- 0
                    tscale <- 1 / iqrfactor
                }
                ## tlocation <- log(domainmax - minfo$centralvalue)
                ## tscale <- abs(
                ##   log(domainmax - minfo$highvalue) -
                ##   log(domainmax - minfo$lowvalue)
                ## ) / iqrfactorLog
                closeddomain <- TRUE
                censormin <- domainmin
                censormax <- domainmax
                tcensormin <- (log(domainmax - censormin) - tlocation) / tscale
                tcensormax <- (log(domainmax - censormax) - tlocation) / tscale
            } else if (is.finite(minfo$domainmin) && is.finite(minfo$domainmax) &&
                       !minfo$minincluded && minfo$maxincluded) {
                ## doubly-bounded right-closed domain
                transf <- 'log'
                if(!is.null(data)) {
                    tlocation <- median(log(x - domainmin))
                    tscale <- mad(log(x - domainmin)) / iqrfactor
                } else {
                    tlocation <- 0
                    tscale <- 1 / iqrfactor
                }
                ## tlocation <- log(minfo$centralvalue - domainmin)
                ## tscale <- abs(
                ##   log(minfo$highvalue - domainmin) -
                ##   log(minfo$lowvalue - domainmin)
                ## ) / iqrfactorLog
                closeddomain <- TRUE
                censormin <- domainmin
                censormax <- domainmax
                tcensormin <- (log(censormin - domainmin) - tlocation) / tscale
                tcensormax <- (log(censormax - domainmin) - tlocation) / tscale
            } else if (is.finite(minfo$domainmin) && is.finite(minfo$domainmax) &&
                       minfo$minincluded && minfo$maxincluded) {
                ## doubly-bounded closed domain
                transf <- 'identity'
                if(!is.null(data)) {
                    tlocation <- median(x)
                    tscale <- mad(x) / iqrfactor
                } else {
                    tlocation <- 0
                    tscale <- 1 / iqrfactor
                }
                ## tlocation <- minfo$centralvalue
                ## tscale <- abs(minfo$highvalue - minfo$lowvalue) / iqrfactor
                closeddomain <- TRUE
                censormin <- domainmin
                censormax <- domainmax
                tcensormin <- (censormin - tlocation) / tscale
                tcensormax <- (censormax - tlocation) / tscale
            } else if (is.finite(minfo$domainmin) && !minfo$minincluded) {
                ## left-bounded left-open domain
                transf <- 'log'
                if(!is.null(data)) {
                    tlocation <- median(log(x - domainmin))
                    tscale <- mad(log(x - domainmin)) / iqrfactor
                } else {
                    tlocation <- 0
                    tscale <- 1 / iqrfactor
                }
                ## tlocation <- log(minfo$centralvalue - domainmin)
                ## tscale <- abs(
                ##   log(minfo$highvalue - domainmin) -
                ##   log(minfo$lowvalue - domainmin)
                ## ) / iqrfactorLog
                closeddomain <- FALSE
                censormin <- domainmin
                censormax <- domainmax
                tcensormin <- (log(censormin - domainmin) - tlocation) / tscale
                tcensormax <- (log(censormax - domainmin) - tlocation) / tscale
            } else if (is.finite(minfo$domainmin) && minfo$minincluded) {
                ## left-bounded left-closed domain
                transf <- 'identity'
                if(!is.null(data)) {
                    tlocation <- median(x)
                    tscale <- mad(x) / iqrfactor
                } else {
                    tlocation <- 0
                    tscale <- 1 / iqrfactor
                }
                ## tlocation <- minfo$centralvalue
                ## tscale <- abs(minfo$highvalue - minfo$lowvalue) / iqrfactor
                closeddomain <- TRUE
                censormin <- domainmin
                censormax <- domainmax
                domainmin <- -Inf
                tcensormin <- (censormin - tlocation) / tscale
                tcensormax <- (censormax - tlocation) / tscale
            } else if (is.finite(minfo$domainmax) && !minfo$maxincluded) {
                ## right-bounded right-open domain
                transf <- 'logminus'
                if(!is.null(data)) {
                    tlocation <- median(log(domainmax - x))
                    tscale <- mad(log(domainmax - x)) / iqrfactor
                } else {
                    tlocation <- 0
                    tscale <- 1 / iqrfactor
                }
                ## tlocation <- log(domainmax - minfo$centralvalue)
                ## tscale <- abs(
                ##   log(domainmax - minfo$highvalue) -
                ##   log(domainmax - minfo$lowvalue)
                ## ) / iqrfactorLog
                closeddomain <- FALSE
                censormin <- domainmin
                censormax <- domainmax
                tcensormin <- (log(domainmax - censormin) - tlocation) / tscale
                tcensormax <- (log(domainmax - censormax) - tlocation) / tscale
            } else if (is.finite(minfo$domainmax) && minfo$maxincluded) {
                ## right-bounded right-closed domain
                transf <- 'identity'
                if(!is.null(data)) {
                    tlocation <- median(x)
                    tscale <- mad(x) / iqrfactor
                } else {
                    tlocation <- 0
                    tscale <- 1 / iqrfactor
                }
                ## tlocation <- minfo$centralvalue
                ## tscale <- abs(minfo$highvalue - minfo$lowvalue) / iqrfactor
                closeddomain <- TRUE
                censormin <- domainmin
                censormax <- domainmax
                domainmax <- +Inf
                tcensormin <- (censormin - tlocation) / tscale
                tcensormax <- (censormax - tlocation) / tscale
            } else {
                ## unbounded, open domain
                transf <- 'identity'
                if(!is.null(data)) {
                    tlocation <- median(x)
                    tscale <- mad(x) / iqrfactor
                } else {
                    tlocation <- 0
                    tscale <- 1 / iqrfactor
                }
                ## tlocation <- minfo$centralvalue
                ## tscale <- abs(minfo$highvalue - minfo$lowvalue) / iqrfactor
                closeddomain <- FALSE
                censormin <- domainmin
                censormax <- domainmax
                tcensormin <- (censormin - tlocation) / tscale
                tcensormax <- (censormax - tlocation) / tscale
            }
            ## # Old cases
            ## if (is.finite(minfo$domainmin) && is.finite(minfo$domainmax)) {
            ##   transf <- 'Q'
            ##   if (minfo$minincluded && !minfo$maxincluded) {
            ##     censormin <- domainmin
            ##     censormax <- +Inf
            ##     domainmin <- (8 * domainmin - domainmax) / 7
            ##   } else if (!minfo$minincluded && minfo$maxincluded) {
            ##     censormax <- domainmax
            ##     censormin <- -Inf
            ##     domainmax <- (8 * domainmax - domainmin) / 7
            ##   } else if (minfo$minincluded && minfo$maxincluded) {
            ##     censormin <- domainmin
            ##     censormax <- domainmax
            ##     domainmin <- (7 * domainmin - domainmax) / 6
            ##     domainmax <- (7 * domainmax - domainmin) / 6
            ##   }
            ##   tlocation <- Qf((tlocation - domainmin) / (domainmax - domainmin))
            ##   tscale <- abs(Qf((minfo$highvalue - domainmin) / (domainmax - domainmin)) -
            ##                  Qf((minfo$lowvalue - domainmin) / (domainmax - domainmin))) *
            ##     sdoveriqr
            ## } else if (is.finite(minfo$domainmin)) {
            ##   transf <- 'log'
            ##   tlocation <- log(tlocation - domainmin)
            ##   tscale <- abs(log(minfo$highvalue - domainmin) -
            ##                  log(minfo$lowvalue - domainmin)) * sdoveriqr
            ## } else if (is.finite(minfo$domainmax)) {
            ##   transf <- 'logminus'
            ##   tlocation <- log(domainmax - tlocation)
            ##   tscale <- abs(log(domainmax - minfo$highvalue) -
            ##                log(domainmax - minfo$lowvalue)) * sdoveriqr
            ## }


            if (step > 0) { # discretized
                mcmctype <- 'D'
                id <- idD
                idD <- idD + 1L
            } else if (closeddomain) { # censored
                mcmctype <- 'C'
                id <- idC
                idC <- idC + 1L
            } else { # continuous
                mcmctype <- 'R'
                id <- idR
                idR <- idR + 1L
            }
        } else { # end continuous case
            stop(paste0('ERROR: unknown variate type for ', name))
        }
        ## Debugging
        ## print(auxmetadata[nrow(auxmetadata)])
        ## print(as.data.table(c(list(name=name, type=mcmctype, transform=transf,
        ## Nvalues=Nvalues, step=step, domainmin=domainmin, domainmax=domainmax,
        ## censormin=censormin, censormax=censormax, tlocation=tlocation,
        ## tscale=tscale, plotmin=plotmin, plotmax=plotmax, Q1=Q1, Q2=Q2, Q3=Q3),
        ## datavalues
        ## )))

        auxmetadata <- rbind(auxmetadata,
            c(
                list(
                    name = name, mcmctype = mcmctype, id = id, # censored=cens,
                    transform = transf, Nvalues = Nvalues, step = step,
                    domainmin = domainmin, domainmax = domainmax,
                    censormin = censormin, censormax = censormax,
                    tcensormin = tcensormin, tcensormax = tcensormax,
                    tlocation = tlocation, tscale = tscale,
                    plotmin = plotmin, plotmax = plotmax
                    ## Q1 = Q1, Q2 = Q2, Q3 = Q3, # not used in other scripts, possibly remove
                    ## mctest1 = mctest1, mctest2 = mctest2, mctest3 = mctest3
                ),
                datavalues
            )
            ## ## previously, with data.table:
            ##      fill = FALSE
        )
    }

    auxmetadata
    ## if (!missing(file) && file != FALSE) { # must save to file
    ##   if (is.character(file)) {
    ##     file <- paste0(sub('.rds$', '', file), '.rds')
    ##   } else {
    ##     file <- paste0('auxmetadata_', datafile)
    ##     file <- paste0(sub('.csv$', '', file), '.rds')
    ##   }
    ##   if (file.exists(file)) {
    ##     file.rename(from = file,
    ##                 to = paste0(sub('.rds$', '', file), '_bak',
    ##                             format(Sys.time(), '%y%m%dT%H%M%S'), '.rds'))
    ##   }
    ##   saveRDS(auxmetadata, file)
    ##   cat('Saved proposal aux-metadata file as', file, '\n')
    ## } else {
    ##   auxmetadata
    ## }
}
