#' Build preliminary metadata flie
#'
#' @param data data.frame object
#' @param metadata data.frame object
#' @param Dthreshold positive number: threshold of fraction
#'   of unique datapoints to total datapoints, to decide
#'   whether to treat a rounded variate as continuous
#'
#' @return an auxmetadata data.frame object
buildauxmetadata <- function(data, metadata, Dthreshold = 1) {

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
    ##  iqrfactorLog <- 2 * 6
    ##  iqrfactorQ <- 2 * 0.3

    Othreshold <- 10

    idR <- idC <- idD <- idB <- idO <- idN <- 1L

#### Prepare empty metadata frame
    auxmetadata <- as.data.frame(
                list(
                    name = NA, mcmctype = NA, id = NA, # censored= NA,
                    transform = NA, Nvalues = NA, halfstep = NA,
                    domainmin = NA, domainmax = NA,
                    tdomainmin = NA, tdomainmax = NA,
                    leftbound = NA, rightbound = NA,
                    tleftbound = NA, trightbound = NA,
                    tlocation = NA, tscale = NA,
                    plotmin = NA, plotmax = NA
                )
    )[-1,]

    for (name in metadata$name) {
        if(!is.null(data)) {
            x <- data[[name]]
            x <- x[!is.na(x)]
        }
        minfo <- as.list(metadata[metadata$name == name, ])
        ## make sure 'type' is lowercase
        minfo$type <- tolower(minfo$type)
        transf <- 'identity' # temporary
        datavalues <- minfo[grep('^V[0-9]+$', names(minfo))]
        ndatavalues <- sum(!is.na(datavalues))
        domainmin <- minfo$domainmin
        domainmax <- minfo$domainmax
        leftbound <- NA
        tleftbound <- NA
        rightbound <- NA
        trightbound <- NA
        Nvalues <- +Inf
        halfstep <- minfo$datastep / 2
        ## Nvalues <- minfo$Nvalues
        plotmin <- minfo$plotmin
        plotmax <- minfo$plotmax

        ## Catch discrepancies in specification
        if(ndatavalues > 0 && (
            !is.na(domainmin) || !is.na(domainmax) || !is.na(halfstep)
        )) {
            stop('Variate "',
                name, '": Please specify the domain either as a list of values',
                '(fields "V1", "V2", ...)\nor as a combination of',
                '"domainmin", "domainmax", "datastep";',
                'not both.')
        }
        if(minfo$type %in% c('nominal', 'ordinal')) {
            if(ndatavalues == 0 && (
            is.na(domainmin) || is.na(domainmax) || is.na(halfstep)
            )) {
            stop('Variate "',
                name, '": Insufficient domain information.\n',
                'Please specify domain either as a list of values',
                '(fields "V1", "V2", ...)\nor as a combination of',
                '"domainmin", "domainmax", "datastep".')
            }
            if(ndatavalues > 0){
                Nvalues <- ndatavalues
            } else {
                Nvalues <- length(seq(domainmin, domainmax, by = 2 * halfstep))
            }
        }

        ## Identify binary variate
        if(minfo$type %in% c('nominal', 'ordinal') && Nvalues == 2) {
            mcmctype <- 'B'
            id <- idB
            idB <- idB + 1L
            if(ndatavalues != 2) {
                datavalues <- as.list(as.character(
                    c(domainmin, domainmax)
                ))
                names(datavalues) <- paste0('V', seq_len(2))
            }
            halfstep <- NA
            domainmin <- NA
            domainmax <- NA
            tdomainmin <- 0 # r/dbeta takes 0/1 inputs
            tdomainmax <- 1
            tlocation <- 0
            tscale <- 1
            plotmin <- NA
            plotmax <- NA

        } else if (minfo$type == 'nominal') {
            ## nominal variate
            mcmctype <- 'N'
            id <- idN
            idN <- idN + 1L
            halfstep <- NA
            domainmin <- NA
            domainmax <- NA
            tdomainmin <- 1  # Nimble indexes categorical from 1
            tdomainmax <- Nvalues
            tlocation <- 0
            tscale <- 1
            plotmin <- NA
            plotmax <- NA

        } else if (minfo$type == 'ordinal') {
            ## Ordinal variates can be specified in two different ways

### Consistency checks:
            if(ndatavalues > Othreshold) {
                stop('Variate "',
                    name, '": we cannot handle lists of values larger than ',
                    Othreshold,'.',
                    '\nPlease specifies its values as a numeric range')
            }

            ## Ordinal variates are internally handled in two different ways:
            ## if their values are less than Othreshold: type 'O', like nominal
            ## otherwise: type 'D', via a latent variable like rounded

### ordinal type 'O'
            if(Nvalues <= Othreshold) {
                mcmctype <- 'O'
                id <- idO
                idO <- idO + 1L
                ## convert domain specification to list
                if(ndatavalues == 0) {
                    datavalues <- as.list(as.character(
                        seq(domainmin, domainmax, length.out = Nvalues)
                    ))
                    names(datavalues) <- paste0('V', seq_len(Nvalues))
                }
                domainmin <- NA
                domainmax <- NA
                halfstep <- NA
                tdomainmin <- 1 # Nimble indexes categorical from 1
                tdomainmax <- Nvalues
                tlocation <- 0
                tscale <- 1
                plotmin <- NA
                plotmax <- NA

### ordinal type 'D'
            } else {
                mcmctype <- 'D'
                id <- idD
                idD <- idD + 1L
                transf <- 'identity'
                if(!is.null(data)) {
                    ## set the location to the nearest centre of datastep interval
                    tempvalue <- unname(quantile(x, 0.5, type = 6))
                    tlocation <- x[which.min(abs(x - tempvalue))]
                    tlocation <- tlocation +
                        round((tempvalue - tlocation) / (2 * halfstep)) *
                        2 * halfstep
                    tscale <- mad(x) / iqrfactor
                } else {
                    tlocation <- mean(domainmin, domainmax)
                    tscale <- ((domainmax - domainmin) / sqrt(12)) / iqrfactor
                }
                tdomainmin <-(domainmin - tlocation) / tscale
                tdomainmax <-(domainmax - tlocation) / tscale

                leftbound <- domainmin + halfstep
                rightbound <- domainmax - halfstep
                tleftbound <- (leftbound - tlocation) / tscale
                trightbound <- (rightbound - tlocation) / tscale

                if(is.finite(minfo$plotmin)){
                    plotmin <- minfo$plotmin
                } else {
                    plotmin <- domainmin
                }
                ## centre plotmin on a rounding bin
                tempvalue <- plotmin
                plotmin <- x[which.min(abs(x - tempvalue))]
                plotmin <- plotmin +
                    round((tempvalue - plotmin) / (2 * halfstep)) *
                    2 * halfstep

                if(is.finite(minfo$plotmax)){
                    plotmax <- minfo$plotmax
                } else {
                    plotmax <- domainmax
                }
                ## centre plotmax on a rounding bin
                tempvalue <- plotmax
                plotmax <- x[which.min(abs(x - tempvalue))]
                plotmax <- plotmax +
                    round((tempvalue - plotmax) / (2 * halfstep)) *
                    2 * halfstep

            }

        } else if (minfo$type == 'continuous') { # continuous variate (R,C,D)
            Nvalues <- +Inf
            ## Rounded variates
            if (is.null(minfo$datastep) || is.na(minfo$datastep)) {
                halfstep <- 0
            } else {
                halfstep <- minfo$datastep / 2
            }
            ## If the variate is rounded,
            ## we avoid a  latent-variable representation
            ## if the datapoints are "enough distinct" in higher dimension
            if(halfstep > 0 &&
               (is.null(data) ||
                nrow(unique(data))/nrow(data) > Dthreshold)) {
                halfstep <- 0
            }

### Rounded continuous case
            if (halfstep > 0) {

                mcmctype <- 'D'
                id <- idD
                idD <- idD + 1L
                ##
                transf <- 'identity'
                if(!is.null(data)) {
                    ## set the location to the nearest centre of datastep interval
                    tempvalue <- unname(quantile(x, 0.5, type = 6))
                    tlocation <- x[which.min(abs(x - tempvalue))]
                    tlocation <- tlocation +
                        round((tempvalue - tlocation) / (2 * halfstep)) *
                        2 * halfstep
                    tscale <- mad(x) / iqrfactor
                } else {
                    tlocation <- 0
                    tscale <- 1 / iqrfactor
                }
                tdomainmin <-(domainmin - tlocation) / tscale
                tdomainmax <-(domainmax - tlocation) / tscale

                ## For a rounded variate it does not matter whether
                ## a domain boundary is included or excluded
                ## hence 'domainmin/max' are ignored
                leftbound <- domainmin + halfstep
                rightbound <- domainmax - halfstep
                ## if(minfo$minincluded){
                ## } else {
                ##     leftbound <- domainmin + halfstep * 3
                ## }
                ## if(minfo$maxincluded){
                ## } else {
                ##     rightbound <- domainmin - halfstep * 3
                ## }
                tleftbound <- (leftbound - tlocation) / tscale
                trightbound <- (rightbound - tlocation) / tscale

            } else {
### Non-rounded continuous cases
                ## variate transformation to real-line domain
                if (is.finite(minfo$domainmin) && is.finite(minfo$domainmax) &&
                    minfo$minincluded && !minfo$maxincluded) {
                    ## doubly-bounded left-closed domain
                    mcmctype <- 'C'
                    id <- idC
                    idC <- idC + 1L
                    ##
                    transf <- 'logminus'
                    if(!is.null(data)) {
                        tlocation <- unname(quantile(-log(domainmax - x), 0.5, type = 6))
                        tscale <- mad(-log(domainmax - x)) / iqrfactor
                    } else {
                        tlocation <- 0
                        tscale <- 1 / iqrfactor
                    }
                    tdomainmin <- (-log(domainmax - domainmin) -
                                                 tlocation) / tscale
                    tdomainmax <- +Inf


                } else if (is.finite(minfo$domainmin) && is.finite(minfo$domainmax) &&
                           !minfo$minincluded && minfo$maxincluded) {
                    ## doubly-bounded right-closed domain
                    mcmctype <- 'C'
                    id <- idC
                    idC <- idC + 1L
                    ##
                    transf <- 'log'
                    if(!is.null(data)) {
                        tlocation <- unname(quantile(log(x - domainmin), 0.5, type = 6))
                        tscale <- mad(log(x - domainmin)) / iqrfactor
                    } else {
                        tlocation <- 0
                        tscale <- 1 / iqrfactor
                    }
                    tdomainmin <- -Inf
                    tdomainmax <- (log(domainmax - domainmin) - tlocation) / tscale

                } else if (is.finite(minfo$domainmin) && is.finite(minfo$domainmax) &&
                           minfo$minincluded && minfo$maxincluded) {
                    ## doubly-bounded closed domain
                    mcmctype <- 'C'
                    id <- idC
                    idC <- idC + 1L
                    ##
                    transf <- 'identity'
                    if(!is.null(data)) {
                        tlocation <- unname(quantile(x, 0.5, type = 6))
                        tscale <- mad(x) / iqrfactor
                    } else {
                        tlocation <- 0
                        tscale <- 1 / iqrfactor
                    }
                    tdomainmin <- (domainmin - tlocation) / tscale
                    tdomainmax <- (domainmax - tlocation) / tscale

                } else if (is.finite(minfo$domainmin) && minfo$minincluded) {
                    ## left-bounded left-closed domain
                    mcmctype <- 'C'
                    id <- idC
                    idC <- idC + 1L
                    ##
                    transf <- 'identity'
                    if(!is.null(data)) {
                        tlocation <- unname(quantile(x, 0.5, type = 6))
                        tscale <- mad(x) / iqrfactor
                    } else {
                        tlocation <- 0
                        tscale <- 1 / iqrfactor
                    }
                    tdomainmin <- (domainmin - tlocation) / tscale
                    tdomainmax <- +Inf

                } else if (is.finite(minfo$domainmax) && minfo$maxincluded) {
                    ## right-bounded right-closed domain
                    mcmctype <- 'C'
                    id <- idC
                    idC <- idC + 1L
                    ##
                    transf <- 'identity'
                    if(!is.null(data)) {
                        tlocation <- unname(quantile(x, 0.5, type = 6))
                        tscale <- mad(x) / iqrfactor
                    } else {
                        tlocation <- 0
                        tscale <- 1 / iqrfactor
                    }
                    tdomainmin <- -Inf
                    tdomainmax <- (domainmax - tlocation) / tscale

                } else if (is.finite(minfo$domainmin) && is.finite(minfo$domainmax) &&
                           !minfo$minincluded && !minfo$maxincluded) {
                    ## doubly-bounded open domain
                    mcmctype <- 'R'
                    id <- idR
                    idR <- idR + 1L
                    ##
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
                    tdomainmin <- -Inf
                    tdomainmax <- +Inf

                } else if (is.finite(minfo$domainmin) && !minfo$minincluded) {
                    ## left-bounded left-open domain
                    mcmctype <- 'R'
                    id <- idR
                    idR <- idR + 1L
                    ##
                    transf <- 'log'
                    if(!is.null(data)) {
                        tlocation <- unname(quantile(log(x - domainmin), 0.5, type = 6))
                        tscale <- mad(log(x - domainmin)) / iqrfactor
                    } else {
                        tlocation <- 0
                        tscale <- 1 / iqrfactor
                    }
                    tdomainmin <- -Inf
                    tdomainmax <- +Inf

                } else if (is.finite(minfo$domainmax) && !minfo$maxincluded) {
                    ## right-bounded right-open domain
                    mcmctype <- 'R'
                    id <- idR
                    idR <- idR + 1L
                    ##
                    transf <- 'logminus'
                    if(!is.null(data)) {
                        tlocation <- unname(quantile(-log(domainmax - x), 0.5, type = 6))
                        tscale <- mad(-log(domainmax - x)) / iqrfactor
                    } else {
                        tlocation <- 0
                        tscale <- 1 / iqrfactor
                    }
                    tdomainmin <- -Inf
                    tdomainmax <- +Inf

                } else {
                    ## unbounded, open domain
                    mcmctype <- 'R'
                    id <- idR
                    idR <- idR + 1L
                    ##
                    transf <- 'identity'
                    if(!is.null(data)) {
                        tlocation <- unname(quantile(x, 0.5, type = 6))
                        tscale <- mad(x) / iqrfactor
                    } else {
                        tlocation <- 0
                        tscale <- 1 / iqrfactor
                    }
                    tdomainmin <- -Inf
                    tdomainmax <- +Inf
                }

                leftbound <- domainmin
                rightbound <- domainmax
                tleftbound <- tdomainmin
                trightbound <- tdomainmax
            }
            ## end non-rounded case

            if(is.finite(minfo$plotmin)){
                plotmin <- minfo$plotmin
            } else {
                plotmin <- max(domainmin,
                    if(!is.null(data)){
                        min(x) - IQR(x, type = 6) / 2
                    } else {
                        -2 * iqrfactor
                    })
            }

            if(is.finite(minfo$plotmax)){
                plotmax <- minfo$plotmax
            } else {
                plotmax <- min(domainmax,
                    if(!is.null(data)){
                        max(x) + IQR(x, type = 6) / 2
                    } else {
                        2 * iqrfactor
                    })
            }

### end continuous case
        } else {
            stop('ERROR: unknown variate type for "', name, '"')
        }
        ## Debugging
        ## print(auxmetadata[nrow(auxmetadata)])
        ## print(as.data.table(c(list(name=name, type=mcmctype, transform=transf,
        ## Nvalues=Nvalues, halfstep=halfstep,
        ## domainmin=domainmin, domainmax=domainmax,
        ## censormin=censormin, censormax=censormax, tlocation=tlocation,
        ## tscale=tscale, plotmin=plotmin, plotmax=plotmax, Q1=Q1, Q2=Q2, Q3=Q3),
        ## datavalues
        ## )))
        auxmetadata <- merge(auxmetadata,
            c(
                list(
                    name = name, mcmctype = mcmctype, id = id, # censored=cens,
                    transform = transf, Nvalues = Nvalues, halfstep = halfstep,
                    domainmin = domainmin, domainmax = domainmax,
                    tdomainmin = tdomainmin, tdomainmax = tdomainmax,
                    leftbound = leftbound, rightbound = rightbound,
                    tleftbound = tleftbound, trightbound = trightbound,
                    tlocation = tlocation, tscale = tscale,
                    plotmin = plotmin, plotmax = plotmax
                    ## Q1 = Q1, Q2 = Q2, Q3 = Q3, # not used in other scripts, possibly remove
                    ## mctest1 = mctest1, mctest2 = mctest2, mctest3 = mctest3
                ),
                as.list(datavalues)
            ),
            sort = FALSE, all = TRUE)
    }

    ## print(auxmetadata) # for debugging
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
