#' Build preliminary metadata flie
#'
#' @param data data.frame object
#' @param metadata data.frame object
#' @param Dthreshold Positive number: threshold of fraction
#'   of unique datapoints to total datapoints, to decide
#'   whether to treat a rounded variate as continuous
#' @param tscalefactor Positive number: scaling factor for variate conversion
#'
#' @return an auxmetadata data.frame object
#' @keywords internal
buildauxmetadata <- function(data, metadata, Dthreshold = 1, tscalefactor = 4.266) {

    ## In the internal, rescaled representation,
    ## with the SD of the means equal to 3 and
    ## the rate of the variances equal to 1,
    ## the median of Q3 of the possible F distribs is 1.95
    ## and its mean is 2.55, while
    ## the median of the IQR is 3.35
    ## and its mean is 5.10.
    ## The MAD of the prior for the next unit is 3.67
    ## The IQR of the prior for the next unit is 4.96
    ## The prior median of the IQRs of the components is 1.349
    ## The prior quartiles of the IQRs of the components are 0.559 & 3.25
    ## The prior 0.89-quantile of the IQRs of the components is 7.72
    ## half-order of magnitude smaller than the IQR-median is 4.26585
    ## which is approx 80.5%-quantile of the component-IQR distribution

    Othreshold <- 10

    idR <- idC <- idD <- idB <- idO <- idN <- 1L
    posN <- posO <- 0L

#### Prepare empty metadata frame
    auxmetadata <- as.data.frame(
                list(
                    name = NA, type = NA, mcmctype = NA, id = NA, # censored= NA,
                    transform = NA,
                    Nvalues = NA_integer_, # keep this as integer
                    indexpos = NA_integer_, # keep this as integer
                    halfstep = NA,
                    domainmin = NA, domainmax = NA,
                    minincluded = NA, maxincluded = NA,
                    tdomainmin = NA, tdomainmax = NA,
                    domainminplushs = NA, domainmaxminushs = NA,
                    tdomainminplushs = NA, tdomainmaxminushs = NA,
                    tlocation = NA, tscale = NA,
                    plotmin = NA, plotmax = NA
                )
    )[-1, ]

    for (name in metadata$name) {
        minfo <- as.list(metadata[metadata$name == name, ])
        ## make sure 'type' is lowercase
        minfo$type <- tolower(minfo$type)
        transf <- 'identity' # temporary
        datavalues <- minfo[grep('^V[0-9]+$', names(minfo))]
        ndatavalues <- sum(!is.na(datavalues))
        domainmin <- as.numeric(minfo$domainmin)
        domainmax <- as.numeric(minfo$domainmax)
        domainminplushs <- NA
        tdomainminplushs <- NA
        domainmaxminushs <- NA
        tdomainmaxminushs <- NA
        Nvalues <- NA_integer_
        indexpos <- NA_integer_
        halfstep <- as.numeric(minfo$datastep) / 2
        if(!is.null(data)) {
            x <- data[[name]]
            x <- x[!is.na(x)]
        }
        ## Nvalues <- minfo$Nvalues
        ## plotmin <- minfo$plotmin
        ## plotmax <- minfo$plotmax

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
            minincluded <- NA
            maxincluded <- NA

        } else if (minfo$type == 'nominal') {
            ## nominal variate
            mcmctype <- 'N'
            id <- idN
            idN <- idN + 1L
            indexpos <- posN
            posN <- posN + Nvalues
            halfstep <- NA
            domainmin <- NA
            domainmax <- NA
            tdomainmin <- 1  # Nimble indexes categorical from 1
            tdomainmax <- Nvalues
            tlocation <- 0
            tscale <- 1
            plotmin <- NA
            plotmax <- NA
            minincluded <- NA
            maxincluded <- NA

        } else if (minfo$type == 'ordinal') {
            ## Ordinal variates can be specified in two different ways

### Consistency checks:
            if(ndatavalues > Othreshold) {
                stop('Variate "',
                    name, '": we cannot handle lists of values larger than ',
                    Othreshold, '.',
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
                indexpos <- posO
                posO <- posO + Nvalues
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
                minincluded <- TRUE
                maxincluded <- TRUE

### ordinal type 'D'
            } else {
                mcmctype <- 'D'
                id <- idD
                idD <- idD + 1L
                transf <- 'identity'
                if(!is.null(data)) {
                    ## set the location to the nearest centre of datastep interval
                    tempvalue <- quantile(x, probs = 0.5, type = 6, names = FALSE)
                    tlocation <- x[which.min(abs(x - tempvalue))]
                    tlocation <- tlocation +
                        round((tempvalue - tlocation) / (2 * halfstep)) *
                        2 * halfstep
                    tscale <- IQR(type = 6, na.rm = TRUE, x = x) / tscalefactor
                    if(tscale == 0){tscale <- 1 / tscalefactor}
                    plotmin <- min(x[is.finite(x)])
                    plotmin <- max(plotmin - IQR(x, type = 6) / 2,
                    (plotmin + domainmin) / 2, na.rm = TRUE)
                    plotmax <- max(x[is.finite(x)])
                    plotmax <- min(plotmax + IQR(x, type = 6) / 2,
                    (plotmax + domainmax) / 2, na.rm = TRUE)
                } else {
                    tlocation <- mean(domainmin, domainmax)
                    tscale <- ((domainmax - domainmin) / sqrt(12)) / tscalefactor
                    plotmin <- (31 * domainmin + domainmax) / 32
                    plotmax <- (31 * domainmax + domainmin) / 32
                }
                tdomainmin <- (domainmin - tlocation) / tscale
                tdomainmax <- (domainmax - tlocation) / tscale

                domainminplushs <- domainmin + halfstep
                domainmaxminushs <- domainmax - halfstep
                tdomainminplushs <- (domainminplushs - tlocation) / tscale
                tdomainmaxminushs <- (domainmaxminushs - tlocation) / tscale

                ## plotmin <- max(domainmin, plotmin)
                ## centre plotmin on a rounding bin
                tempvalue <- plotmin
                plotmin <- x[which.min(abs(x - tempvalue))]
                plotmin <- plotmin +
                    round((tempvalue - plotmin) / (2 * halfstep)) *
                    2 * halfstep

                ## plotmax <- min(domainmax, plotmax)
                ## centre plotmax on a rounding bin
                tempvalue <- plotmax
                plotmax <- x[which.min(abs(x - tempvalue))]
                plotmax <- plotmax +
                    round((tempvalue - plotmax) / (2 * halfstep)) *
                    2 * halfstep

                minincluded <- TRUE
                maxincluded <- TRUE
            }

#### Continuous variate (R,C,D)
        } else if (minfo$type == 'continuous') {
            Nvalues <- NA_integer_
            ## Rounded variates
            ## Empty datastep means 0
            if (is.null(halfstep) || is.na(halfstep)) {
                halfstep <- 0
            }
            ## Empty domainmin/max values mean -+Inf
            if(is.null(domainmin) || is.na(domainmin)) {
                domainmin <- -Inf
            }
            if(is.null(domainmax) || is.na(domainmax)) {
                domainmax <- +Inf
            }
            minincluded <- (tolower(minfo$minincluded) %in%
                                c('true', 't', 'yes', 'y', '1'))
            maxincluded <- (tolower(minfo$maxincluded) %in%
                                c('true', 't', 'yes', 'y', '1'))
            ## If the variate is rounded,
            ## we avoid a  latent-variable representation
            ## if the datapoints are "enough distinct" in higher dimension
            if(halfstep > 0 &&
               (is.null(data) ||
                nrow(unique(data)) / nrow(data) > Dthreshold)) {
                halfstep <- 0
            }

            if(!is.null(data)) {
                plotmin <- min(x[is.finite(x)])
                plotmax <- max(x[is.finite(x)])
            }

            ## Right now, some algorithms to set particular values
            ## are the same for all cases below,
            ## so they could be collected here
            ## But they are left in the individual cases
            ## in the event that they need a case-by-case tretment
            ## in the future

### Rounded continuous case
            if (halfstep > 0) {

                mcmctype <- 'D'
                id <- idD
                idD <- idD + 1L
                ##
                transf <- 'identity'
                if(!is.null(data)) {
                    ## set the location to the nearest centre of datastep interval
                    tempvalue <- quantile(x, probs = 0.5, type = 6, names = FALSE)
                    tlocation <- x[which.min(abs(x - tempvalue))]
                    tlocation <- tlocation +
                        round((tempvalue - tlocation) / (2 * halfstep)) *
                        2 * halfstep
                    tscale <- IQR(type = 6, na.rm = TRUE, x = x) / tscalefactor
                    if(tscale == 0){tscale <- 1 / tscalefactor}

                    plotmin <- max(plotmin - IQR(x, type = 6) / 2,
                    (plotmin + domainmin) / 2, na.rm = TRUE)
                    plotmax <- min(plotmax + IQR(x, type = 6) / 2,
                    (plotmax + domainmax) / 2, na.rm = TRUE)

                } else {
                    tlocation <- 0
                    tscale <- 1 / tscalefactor
                    plotmin <- max((31 * domainmin + domainmax) / 32,
                        -3 * tscale)
                    plotmax <- min((31 * domainmax + domainmin) / 32,
                        3 * tscale)
                }
                tdomainmin <- (domainmin - tlocation) / tscale
                tdomainmax <- (domainmax - tlocation) / tscale

                ## For a rounded variate it does not matter whether
                ## a domain boundary is included or excluded
                ## hence 'domainmin/max' are ignored
                domainminplushs <- domainmin + halfstep
                domainmaxminushs <- domainmax - halfstep
                ## if(minfo$minincluded){
                ## } else {
                ##     domainminplushs <- domainmin + halfstep * 3
                ## }
                ## if(minfo$maxincluded){
                ## } else {
                ##     domainmaxminushs <- domainmin - halfstep * 3
                ## }
                tdomainminplushs <- (domainminplushs - tlocation) / tscale
                tdomainmaxminushs <- (domainmaxminushs - tlocation) / tscale

                ## plotmin <- max(domainmin, plotmin)
                ## centre plotmin on a rounding bin
                tempvalue <- plotmin
                plotmin <- x[which.min(abs(x - tempvalue))]
                plotmin <- plotmin +
                    round((tempvalue - plotmin) / (2 * halfstep)) *
                    2 * halfstep

                ## plotmax <- min(domainmax, plotmax)
                ## centre plotmax on a rounding bin
                tempvalue <- plotmax
                plotmax <- x[which.min(abs(x - tempvalue))]
                plotmax <- plotmax +
                    round((tempvalue - plotmax) / (2 * halfstep)) *
                    2 * halfstep

            } else {
### Non-rounded continuous cases
                ## variate transformation to real-line domain
                if (is.finite(domainmin) && is.finite(domainmax) &&
                    minincluded && !maxincluded) {
                    ## doubly-bounded left-closed domain
                    mcmctype <- 'C'
                    id <- idC
                    idC <- idC + 1L
                    ##
                    transf <- 'logminus'
                    if(!is.null(data)) {
                        tlocation <- quantile(-log(domainmax - x), probs = 0.5,
                            type = 6, names = FALSE)
                        tscale <- IQR(type = 6, na.rm = TRUE,
                            x = -log(domainmax - x)) / tscalefactor
                        if(tscale == 0){tscale <- 1 / tscalefactor}

                        plotmin <- max(plotmin - IQR(x, type = 6) / 2,
                            domainmin, na.rm = TRUE)
                        plotmax <- min(plotmax + IQR(x, type = 6) / 2,
                        (plotmax + domainmax) / 2, na.rm = TRUE)

                    } else {
                        tlocation <- 0
                        tscale <- 1 / tscalefactor
                        plotmin <- domainmin
                        plotmax <- (31 * domainmax + domainmin) / 32
                    }
                    tdomainmin <- (-log(domainmax - domainmin) -
                                                 tlocation) / tscale
                    tdomainmax <- +Inf


                } else if (is.finite(domainmin) && is.finite(domainmax) &&
                           !minincluded && maxincluded) {
                    ## doubly-bounded right-closed domain
                    mcmctype <- 'C'
                    id <- idC
                    idC <- idC + 1L
                    ##
                    transf <- 'log'
                    if(!is.null(data)) {
                        tlocation <- quantile(log(x - domainmin), probs = 0.5,
                            type = 6, names = FALSE)
                        tscale <- IQR(type = 6, na.rm = TRUE, x = log(x - domainmin)) / tscalefactor
                        if(tscale == 0){tscale <- 1 / tscalefactor}

                        plotmin <- max(plotmin - IQR(x, type = 6) / 2,
                            (plotmin + domainmin) / 2, na.rm = TRUE)
                        plotmax <- min(plotmax + IQR(x, type = 6) / 2,
                            domainmax, na.rm = TRUE)

                    } else {
                        tlocation <- 0
                        tscale <- 1 / tscalefactor
                        plotmin <- (31 * domainmin + domainmax) / 32
                        plotmax <- domainmax
                    }
                    tdomainmin <- -Inf
                    tdomainmax <- (log(domainmax - domainmin) - tlocation) / tscale

                } else if (is.finite(domainmin) && is.finite(domainmax) &&
                           minincluded && maxincluded) {
                    ## doubly-bounded closed domain
                    mcmctype <- 'C'
                    id <- idC
                    idC <- idC + 1L
                    ##
                    transf <- 'identity'
                    if(!is.null(data)) {
                        tlocation <- quantile(x, probs = 0.5,
                            type = 6, names = FALSE)
                        tscale <- IQR(type = 6, na.rm = TRUE, x = x) / tscalefactor
                        if(tscale == 0){tscale <- 1 / tscalefactor}

                        plotmin <- max(plotmin - IQR(x, type = 6) / 2,
                            domainmin, na.rm = TRUE)
                        plotmax <- min(plotmax + IQR(x, type = 6) / 2,
                            domainmax, na.rm = TRUE)

                    } else {
                        tlocation <- 0
                        tscale <- 1 / tscalefactor
                        plotmin <- domainmin
                        plotmax <- domainmax
                    }
                    tdomainmin <- (domainmin - tlocation) / tscale
                    tdomainmax <- (domainmax - tlocation) / tscale

                } else if (is.finite(domainmin) && minincluded) {
                    ## left-bounded left-closed domain
                    mcmctype <- 'C'
                    id <- idC
                    idC <- idC + 1L
                    ##
                    transf <- 'identity'
                    if(!is.null(data)) {
                        tlocation <- quantile(x, probs = 0.5,
                            type = 6, names = FALSE)
                        tscale <- IQR(type = 6, na.rm = TRUE, x = x) / tscalefactor
                        if(tscale == 0){tscale <- 1 / tscalefactor}

                        plotmin <- max(plotmin - IQR(x, type = 6) / 2,
                            domainmin, na.rm = TRUE)
                        plotmax <- plotmax + IQR(x, type = 6) / 2

                    } else {
                        tlocation <- 0
                        tscale <- 1 / tscalefactor
                        plotmin <- domainmin
                        plotmax <- domainmin + 6
                    }
                    tdomainmin <- (domainmin - tlocation) / tscale
                    tdomainmax <- +Inf

                } else if (is.finite(domainmax) && maxincluded) {
                    ## right-bounded right-closed domain
                    mcmctype <- 'C'
                    id <- idC
                    idC <- idC + 1L
                    ##
                    transf <- 'identity'
                    if(!is.null(data)) {
                        tlocation <- quantile(x, probs = 0.5,
                            type = 6, names = FALSE)
                        tscale <- IQR(type = 6, na.rm = TRUE, x = x) / tscalefactor
                        if(tscale == 0){tscale <- 1 / tscalefactor}

                        plotmin <- plotmin - IQR(x, type = 6) / 2
                        plotmax <- min(plotmax + IQR(x, type = 6) / 2,
                            domainmax, na.rm = TRUE)

                    } else {
                        tlocation <- 0
                        tscale <- 1 / tscalefactor
                        plotmin <- domainmax - 6
                        plotmax <- domainmax
                    }
                    tdomainmin <- -Inf
                    tdomainmax <- (domainmax - tlocation) / tscale

                } else if (is.finite(domainmin) && is.finite(domainmax) &&
                           !minincluded && !maxincluded) {
                    ## doubly-bounded open domain
                    mcmctype <- 'R'
                    id <- idR
                    idR <- idR + 1L
                    ##
                    transf <- 'Q'
                    if(!is.null(data)) {
                        tlocation <- median(
                            util_Q(0.5 +
                                       (x - (domainmax + domainmin) / 2) /
                                       (domainmax - domainmin)
                            ))
                        tscale <- IQR(type = 6, na.rm = TRUE, x =
                            util_Q(0.5 +
                                       (x - (domainmax + domainmin) / 2) /
                                       (domainmax - domainmin)
                            )) / tscalefactor
                        if(tscale == 0){tscale <- 1 / tscalefactor}

                        plotmin <- max(plotmin - IQR(x, type = 6) / 2,
                        (plotmin + domainmin) / 2, na.rm = TRUE)
                        plotmax <- min(plotmax + IQR(x, type = 6) / 2,
                        (plotmax + domainmax) / 2, na.rm = TRUE)

                    } else {
                        tlocation <- 0
                        tscale <- 1 / tscalefactor
                        plotmin <- (31 * domainmin + domainmax) / 32
                        plotmax <- (31 * domainmax + domainmin) / 32

                    }
                    tdomainmin <- -Inf
                    tdomainmax <- +Inf

                } else if (is.finite(domainmin) && !minincluded) {
                    ## left-bounded left-open domain
                    mcmctype <- 'R'
                    id <- idR
                    idR <- idR + 1L
                    ##
                    transf <- 'log'
                    if(!is.null(data)) {
                        tlocation <- quantile(log(x - domainmin), probs = 0.5,
                            type = 6, names = FALSE)
                        tscale <- IQR(type = 6, na.rm = TRUE, x = log(x - domainmin)) / tscalefactor
                        if(tscale == 0){tscale <- 1 / tscalefactor}

                        plotmin <- max(plotmin - IQR(x, type = 6) / 2,
                        (plotmin + domainmin) / 2, na.rm = TRUE)
                        plotmax <- plotmax + IQR(x, type = 6) / 2

                    } else {
                        tlocation <- 0
                        tscale <- 1 / tscalefactor
                        plotmin <- domainmin + 1
                        plotmax <- domainmin + 7
                    }
                    tdomainmin <- -Inf
                    tdomainmax <- +Inf

                } else if (is.finite(domainmax) && !maxincluded) {
                    ## right-bounded right-open domain
                    mcmctype <- 'R'
                    id <- idR
                    idR <- idR + 1L
                    ##
                    transf <- 'logminus'
                    if(!is.null(data)) {
                        tlocation <- quantile(-log(domainmax - x), probs = 0.5,
                            type = 6, names = FALSE)
                        tscale <- IQR(type = 6, na.rm = TRUE, x = -log(domainmax - x)) / tscalefactor
                        if(tscale == 0){tscale <- 1 / tscalefactor}

                        plotmin <- plotmin - IQR(x, type = 6) / 2
                        plotmax <- min(plotmax + IQR(x, type = 6) / 2,
                        (plotmax + domainmax) / 2, na.rm = TRUE)

                    } else {
                        tlocation <- 0
                        tscale <- 1 / tscalefactor
                        plotmin <- domainmax - 7
                        plotmax <- domainmax - 1
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
                        tlocation <- quantile(x, probs = 0.5,
                            type = 6, names = FALSE)
                        tscale <- IQR(type = 6, na.rm = TRUE, x = x) / tscalefactor
                        if(tscale == 0){tscale <- 1 / tscalefactor}

                        plotmin <- plotmin - IQR(x, type = 6) / 2
                        plotmax <- plotmax + IQR(x, type = 6) / 2

                    } else {
                        tlocation <- 0
                        tscale <- 1 / tscalefactor
                        plotmin <- -3
                        plotmax <- 3
                    }
                    tdomainmin <- -Inf
                    tdomainmax <- +Inf
                }

                domainminplushs <- domainmin
                domainmaxminushs <- domainmax
                tdomainminplushs <- tdomainmin
                tdomainmaxminushs <- tdomainmax

            }
            ## end non-rounded case


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
                    name = name, type = minfo$type, mcmctype = mcmctype,
                    id = id, # censored=cens,
                    transform = transf,
                    Nvalues = Nvalues, indexpos = indexpos,
                    halfstep = halfstep,
                    domainmin = domainmin, domainmax = domainmax,
                    minincluded = minincluded, maxincluded = maxincluded,
                    tdomainmin = tdomainmin, tdomainmax = tdomainmax,
                    domainminplushs = domainminplushs,
                    domainmaxminushs = domainmaxminushs,
                    tdomainminplushs = tdomainminplushs,
                    tdomainmaxminushs = tdomainmaxminushs,
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
