#' Plot one-dimensional posterior probabilities
#'
#' @param file string: name of plot output file
#' @param learned Either a string with the name of a directory or full
#'   path for an 'learnt.rds' object, or such an object itself
#' @param data data.table object or filepath: datapoints
#' @param plotprobability logical: plot the resulting probability curve
#' @param plotvariability string, either 'samples' or 'quantiles':
#'   how to plot the variability of the probability distribution with new samples
#' @param nFsamples positive number: if plotvariability='samples', then
#'   number of samples of representative frequency distributions to display
#'   as variability;
#'   if plotvariability='quantiles', then
#'   the quantiles (in range 0 to 0.5) to show
#' @param datahistogram logical: plot the data as histogram?
#' @param datascatter logical: plot the data as scatterplot along the x-axis?
#' @param parallel Bool or numeric: whether to use pre-existing parallel
#'   workers, or how many to create and use
#' @param silent logical: give warnings or updates in the computation
#'
#' @return A list with the mutual information, its error, and its unit
#'
#' @export
plotFsamples <- function(
    file,
    learnt,
    data,
    plotprobability = TRUE,
    plotvariability = 'samples',
    nFsamples = NULL,
    datahistogram = !(missing(data) || is.null(data)),
    datascatter = !(missing(data) || is.null(data)),
    parallel = TRUE,
    silent = FALSE
) {

    fontfamily <- 'Palatino'
    ## Tol colour-blind-friendly palette
    black <- 'black'
    red <- '#EE6677'
    blue <- '#4477AA'
    green <- '#228833'
    yellow <- '#CCBB44'
    purple <- '#AA3377'
    cyan <- '#66CCEE'
    grey <- '#BBBBBB'
    midgrey <- '#888888'

    ## Extract Monte Carlo output & auxmetadata
    ## If 'learnt' is a string, check if it's a folder name or file name
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
                stop('The argument "learnt" must be a folder containing learnt.rds, or the path to an rds-file containing the output from "learn".')
            }
        }
    }

    auxmetadata <- learnt$auxmetadata
    ## learnt$auxmetadata <- NULL
    nsamples <- ncol(learnt$W)

    nodata <- missing(data) || is.null(data) || isFALSE(data)
    if (datahistogram && nodata) {
        datahistogram <- FALSE
        cat('\nNOTE: "datahistogram" is TRUE but there is no data\n')
    }
    if (datascatter && nodata) {
        datascatter <- FALSE
        cat('\nNOTE: "datascatter" is TRUE but there is no data\n')
    }

    addylab <- ''
    if (plotvariability == 'quantiles') {
        if(is.null(nFsamples)) {nFsamples <- c(5.5, 94.5)/100}
        ## plotprobability <- TRUE
        if (any(nFsamples <= 0 | nFsamples >= 1)) {
            nFsamples <- c(5.5, 94.5)/100
        }
        quants <- sort(unique(round(c(nFsamples, 1 - nFsamples), 6)))
        nmcsamples <- NULL
        addylab <- paste0(' & ', ceiling(diff(quants) * 100), '% variability')
    } else {
        if(is.null(nFsamples)) {nFsamples <- 100}
        if (nFsamples == 'all') {
            nFsamples <- nsamples
        }
        quants <- NULL
        nmcsamples <- abs(nFsamples)
    }

    addplot <- FALSE

    graphics.off()
    pdf(file = paste0(file, '.pdf'), family = fontfamily,
        height = 8.27, width = 11.69)
    par(mfrow = c(1, 1))

    for (name in auxmetadata[['name']]) {
        ## Check if we have proper data for this variate
        if (datahistogram || datascatter) {
            theresdata <- (sum(!is.na(data[[name]])) > 1)
        }

        ## save the possibly existing V-values for this variate
        datavalues <- as.list(auxmetadata[auxmetadata$name == name,
            grep('^V[0-9]+$', names(auxmetadata))
        ])
        datavalues <- unlist(datavalues[!is.na(datavalues)])
        nvaluelist <- length(datavalues)

        with(as.list(auxmetadata[auxmetadata$name == name, ]),
        {

            if(mcmctype == 'R') {
                Xgrid <- cbind(seq(plotmin, plotmax, length.out = 256))
                colnames(Xgrid) <- name

                probabilities <- Pr(Y = Xgrid, X = NULL,
                    learnt = learnt,
                    quantiles = quants,
                    nsamples = nmcsamples,
                    parallel = parallel,
                    silent = TRUE,
                    keepYX = FALSE)

                dim(probabilities$values) <- NULL

                ## Find appropriate plot height across plots
                if (plotvariability == 'samples') {
                    dim(probabilities$samples) <- dim(probabilities$samples)[-2]
                    ymax <- quantile(apply(probabilities$samples, 2,
                        function(x) {
                            quantile(c(x), probs = 31 / 32, type = 6,
                                na.rm = TRUE, names = FALSE)
                        }
                    ), probs = 31 / 32, type = 6, na.rm = TRUE, names = FALSE)
                } else {
                    dim(probabilities$quantiles) <- dim(probabilities$quantiles)[-2]
                    ymax <- max(probabilities$quantiles[
                        is.finite(probabilities$quantiles)])
                    ##     apply(probabilities[, , drop = FALSE], 1,
                    ##     function(x) {
                    ##         quantile(x, max(quants), type = 6, na.rm = TRUE)
                    ##     }
                    ## )
                    ## ymax <- max(ymax[is.finite(ymax)])
                }

                ## prepare info for histogram data plots if required
                if (datahistogram && theresdata) {
                    datum <- data[[name]]
                    datum <- datum[!is.na(datum)]
                    ## try to guess optimal bin size
                    nh <- max(32, round(length(datum) / 64))

                    ## histo <- thist(datum, n = nh, extendbreaks = FALSE)
                    histo <- hist(datum, breaks = nh, plot = FALSE)
                    hmax <- max(histo$density)

                    if(ymax/hmax < 0.3 | ymax/hmax > 3) {
                        histo$density <- histo$density/max(histo$density) * ymax
                    } else {
                        ymax <- max(ymax, histo$density)
                    }
                }

                ## If required, plot frequency samples
                if (plotvariability == 'samples') {
                    ## tplot(
                    ##     x = Xgrid, y = probabilities$samples,
                    ##     xlim = range(Xgrid), ylim = c(0, ymax),
                    ##     type = 'l', lty = 1, lwd = 2,
                    ##     col = , alpha = 1/8,
                    ##     xlab = name,
                    ##     ylab = paste0('probability density', addylab),
                    ##     family = fontfamily
                    ## )
                    flexiplot(
                        x = Xgrid, y = probabilities$samples,
                        xlim = range(Xgrid), ylim = c(0, ymax),
                        type = 'l', lty = 1, lwd = 2,
                        col = adjustcolor(cyan, 1/8),
                        xlab = name,
                        ylab = paste0('probability density', addylab),
                        family = fontfamily
                    )
                    addplot <- TRUE # new plots must keep this one
                } else if (plotvariability == 'quantiles') {
                    ## marguncertainty <- t(apply(probabilities, 1, function(x) {
                    ##     quantile(x, quants, type = 6, na.rm = TRUE)
                    ## }))

                    plotquantiles(
                        x = Xgrid,
                        y = probabilities$quantiles,
                        ## y = marguncertainty[, , drop = FALSE],
                        col = cyan, alpha.f = 0.25,
                        xlim = range(Xgrid), ylim = c(0, ymax),
                        xlab = name,
                        ylab = paste0('probability density', addylab),
                        family = fontfamily,
                        add = addplot
                    )
                    addplot <- TRUE # new plots must keep this one
                }

                ## If required, plot probability
                if (plotprobability) {
                    ## tplot(
                    ##     x = Xgrid,
                    ##     y = probabilities$values,
                    ##     ## y = rowMeans(probabilities[, , drop = FALSE], na.rm = TRUE),
                    ##     xlim = range(Xgrid), ylim = c(0, ymax),
                    ##     type = 'l', cex = 0.5, lty = 1, lwd = 4,
                    ##     col = 1, alpha = 0.75,
                    ##     xlab = name,
                    ##     ylab = paste0('probability density', addylab),
                    ##     family = fontfamily,
                    ##     add = addplot
                    ## )
                    flexiplot(
                        x = Xgrid,
                        y = probabilities$values,
                        ## y = rowMeans(probabilities[, , drop = FALSE], na.rm = TRUE),
                        xlim = range(Xgrid), ylim = c(0, ymax),
                        type = 'l', cex = 0.5, lty = 1, lwd = 4,
                        col = adjustcolor(blue, 0.75),
                        xlab = name,
                        ylab = paste0('probability density', addylab),
                        family = fontfamily,
                        add = addplot
                    )
                }

                ## If required and possible, plot data histogram
                if (datahistogram && theresdata) {
                    ## tplot(
                    ##     x = histo$mids, y = histo$density,
                    ##     xlim = range(Xgrid), ylim = c(0, ymax),
                    ##     type = 'l', cex = 0.5, lty = 1, lwd = 2,
                    ##     col = 4, alpha = 0.5,
                    ##     border = '#555555', border.alpha = 1 / 4,
                    ##     xlab = name,
                    ##     ylab = paste0('probability density', addylab),
                    ##     family = fontfamily,
                    ##     add = addplot
                    ## )
                    flexiplot(
                        x = histo$mids, y = histo$density,
                        xlim = range(Xgrid), ylim = c(0, ymax),
                        type = 'l', cex = 0.5, lty = 1, lwd = 2,
                        col = adjustcolor(yellow, 0.5),
                        xlab = name,
                        ylab = paste0('probability density', addylab),
                        family = fontfamily,
                        add = addplot
                    )
                }
                ## End R case

            } else if(mcmctype == 'C') {
                ## C-type variates are peculiar because
                ## their boundary values have a finite probability,
                ## not a probability density.
                ## These values need different handling and representation

                Xgrid <- cbind(seq(plotmin, plotmax, length.out = 256))
                ## }
                colnames(Xgrid) <- name

                xin <- Xgrid > domainminplushs & Xgrid < domainmaxminushs

                probabilities <- Pr(Y = Xgrid, X = NULL,
                    learnt = learnt,
                    quantiles = quants,
                    nsamples = nmcsamples,
                    parallel = parallel,
                    silent = TRUE,
                    keepYX = FALSE)

                dim(probabilities$values) <- NULL

                ## Find appropriate plot height across plots
                if (plotvariability == 'samples') {
                    dim(probabilities$samples) <- dim(probabilities$samples)[-2]
                    ymax <- quantile(apply(
                        probabilities$samples[xin, , drop = FALSE], 2,
                        function(x) {
                            quantile(c(x), probs = 31 / 32, type = 6,
                                na.rm = TRUE, names = FALSE)
                        }
                    ), probs = 31 / 32, type = 6, na.rm = TRUE, names = FALSE)
                    ## ymax <- quantile(apply(
                    ##     probabilities[xin, subsamples, drop = FALSE],
                    ##     2, function(x) {
                    ##         quantile(x, 31 / 32, type = 6, na.rm = TRUE)
                    ##     }
                    ## ), 31 / 32, type = 6, na.rm = TRUE)
                } else {
                    dim(probabilities$quantiles) <- dim(probabilities$quantiles)[-2]
                    temp <- probabilities$quantiles[xin,]
                    ymax <- max(temp[is.finite(temp)])
                    ## ymax <- apply(probabilities[xin, , drop = FALSE], 1,
                    ##     function(x) {
                    ##         quantile(x, max(quants), type = 6, na.rm = TRUE)
                    ##     }
                    ## )
                    ## ymax <- max(ymax[is.finite(ymax)])
                }

                ## prepare info for histogram data plots if required
                if (datahistogram && theresdata) {
                    datum <- data[[name]]
                    datum <- datum[!is.na(datum)]

                    din <- datum > domainminplushs & datum < domainmaxminushs
                    ## try to guess optimal bin size
                    nh <- max(32, round(length(datum[din]) / 64))

                    ## histo <- thist(datum[din], n = nh, extendbreaks = FALSE)
                    histo <- hist(datum, breaks = nh, plot = FALSE)
                    hmax <- max(histo$density)

                    hleft <- sum(datum <= domainminplushs)/length(datum)
                    hright <- sum(datum >= domainmaxminushs)/length(datum)

                    if(ymax/hmax < 0.3 | ymax/hmax > 3) {
                        histo$density <- histo$density/max(histo$density) * ymax
                    } else {
                        ymax <- max(ymax, histo$density)
                    }
                }

                ## If required, plot frequency samples
                if (plotvariability == 'samples') {
                    if (any(xin)) {
                        ## tplot(
                        ##     x = Xgrid[xin],
                        ##     y = probabilities$samples[xin, , drop = FALSE],
                        ##     xlim = range(Xgrid), ylim = c(0, ymax),
                        ##     type = 'l', lty = 1, lwd = 2,
                        ##     col = 5, alpha = 1/8,
                        ##     xlab = name,
                        ##     ylab = paste0('probability density', addylab),
                        ##     family = fontfamily
                        ## )
                        flexiplot(
                            x = Xgrid[xin],
                            y = probabilities$samples[xin, , drop = FALSE],
                            xlim = range(Xgrid), ylim = c(0, ymax),
                            type = 'l', lty = 1, lwd = 2,
                            col = adjustcolor(cyan, 1/8),
                            xlab = name,
                            ylab = paste0('probability density', addylab),
                            family = fontfamily
                        )
                        addplot <- TRUE
                    }

                    ## Boundary points
                    if (any(!xin)) {
                        ## tplot(
                        ##     x = Xgrid[!xin],
                        ##     y = probabilities$values[!xin] * ymax,
                        ##     type = 'p', pch = 2, cex = 2,
                        ##     col = 5, alpha = 1/8,
                        ##     family = fontfamily,
                        ##     add = addplot
                        ## )
                        flexiplot(
                            x = Xgrid[!xin],
                            y = probabilities$values[!xin] * ymax,
                            type = 'p', pch = 2, cex = 2,
                            col = adjustcolor(cyan, 1/8),
                            family = fontfamily,
                            add = addplot
                        )
                        addplot <- TRUE
                    }

                } else if (plotvariability == 'quantiles') {
                    ## marguncertainty <- t(apply(probabilities, 1, function(x) {
                    ##     quantile(x, quants, type = 6, na.rm = TRUE)
                    ## }))

                    if (any(xin)) {
                        plotquantiles(
                            x = Xgrid[xin],
                            y = probabilities$quantiles[xin, , drop = FALSE],
                            col = cyan, alpha.f = 0.25,
                            xlim = range(Xgrid), ylim = c(0, ymax),
                            xlab = name,
                            ylab = paste0('probability density', addylab),
                            family = fontfamily,
                            add = addplot
                        )
                        addplot <- TRUE
                    }

                    ## Boundary points
                    if (any(!xin)) {
                        ## tplot(
                        ##     x = matrix(Xgrid[!xin],
                        ##         nrow = 2, ncol = sum(!xin), byrow = TRUE),
                        ##     y = t(probabilities$quantiles[!xin, , drop = FALSE]) * ymax,
                        ##     type = 'l', pch = 2, cex = 2,
                        ##     col = 5, alpha = 0.25,
                        ##     lty = 1, lwd = 16,
                        ##     add = addplot
                        ## )
                        flexiplot(
                            x = matrix(Xgrid[!xin],
                                nrow = 2, ncol = sum(!xin), byrow = TRUE),
                            y = probabilities$quantiles[!xin, , drop = FALSE] *
                                ymax,
                            type = 'l', pch = 2, cex = 2,
                            col = adjustcolor(cyan, 0.25),
                            lty = 1, lwd = 16,
                            add = addplot
                        )
                        addplot <- TRUE
                    }
                }

                ## If required, plot probability
                if (plotprobability) {
                    if (any(xin)) {
                        ## tplot(
                        ##     x = Xgrid[xin],
                        ##     y = probabilities$values[xin],
                        ##     xlim = range(Xgrid), ylim = c(0, ymax),
                        ##     type = 'l', cex = 0.5, lty = 1, lwd = 4,
                        ##     col = 1, alpha = 0.75,
                        ##     xlab = name,
                        ##     ylab = paste0('probability density', addylab),
                        ##     family = fontfamily,
                        ##     add = addplot
                        ## )
                        flexiplot(
                            x = Xgrid[xin],
                            y = probabilities$values[xin],
                            xlim = range(Xgrid), ylim = c(0, ymax),
                            type = 'l', cex = 0.5, lty = 1, lwd = 4,
                            col = adjustcolor(blue, 0.75),
                            xlab = name,
                            ylab = paste0('probability density', addylab),
                            family = fontfamily,
                            add = addplot
                        )
                        addplot <- TRUE
                    }

                    ## Boundary points
                    if (any(!xin)) {
                        ## tplot(
                        ##     x = Xgrid[!xin],
                        ##     y = probabilities$values[!xin] * ymax,
                        ##     type = 'p', pch = 2, cex = 2,
                        ##     col = 1, alpha = 0.75,
                        ##     lty = 1, lwd = 3,
                        ##     add = addplot
                        ## )
                        flexiplot(
                            x = Xgrid[!xin],
                            y = probabilities$values[!xin] * ymax,
                            type = 'p', pch = 2, cex = 2,
                            col = adjustcolor(blue, 0.75),
                            lty = 1, lwd = 3,
                            add = addplot
                        )
                        addplot <- TRUE
                    }
                }

                ## If required and possible, plot data histogram
                if (datahistogram && theresdata) {
                    if(hleft + hright < 1) {
                        ## tplot(
                        ##     x = histo$mids, y = histo$density,
                        ##     xlim = range(Xgrid), ylim = c(0, ymax),
                        ##     type = 'l', col = 4, alpha = 0.5,
                        ##     border = '#555555', border.alpha = 1/4,
                        ##     xlab = name,
                        ##     ylab = paste0('probability density', addylab),
                        ##     family = fontfamily,
                        ##     add = addplot
                        ## )
                        flexiplot(
                            x = histo$mids, y = histo$density,
                            xlim = range(Xgrid), ylim = c(0, ymax),
                            type = 'l', cex = 0.5, lty = 1, lwd = 2,
                            col = adjustcolor(yellow, 0.5),
                            ## border = '#555555', border.alpha = 1/4,
                            xlab = name,
                            ylab = paste0('probability density', addylab),
                            family = fontfamily,
                            add = addplot
                        )
                        addplot <- TRUE
                    }

                    ## Boundary points
                    if (hleft + hright > 0) {
                        ## tplot(
                        ##     x = c(if(hleft > 0){domainminplushs},
                        ##         if(hright > 0){domainmaxminushs}),
                        ##     y = c(if(hleft > 0){hleft},
                        ##         if(hright > 0){hright}) * ymax,
                        ##     type = 'p', pch = 0, cex = 2,
                        ##     col = 4, alpha = 0.5,
                        ##     lty = 1, lwd = 5,
                        ##     family = fontfamily,
                        ##     add = addplot
                        ## )
                        flexiplot(
                            x = c(if(hleft > 0){domainminplushs},
                                if(hright > 0){domainmaxminushs}),
                            y = c(if(hleft > 0){hleft},
                                if(hright > 0){hright}) * ymax,
                            type = 'p', pch = 0, cex = 2,
                            col = adjustcolor(yellow, 0.5),
                            lty = 1, lwd = 5,
                            family = fontfamily,
                            add = addplot
                        )
                        addplot <- TRUE
                    }
                }
                ## End C case

            } else if (mcmctype %in% c('D', 'O', 'N', 'B')) {
                ## These variate types all have finite probabilities
                if(nvaluelist > 0) {
                    ## Xgrid <-  as.matrix(
                    ##     vtransform(x = datavalues,
                    ##         variates = name,
                    ##         auxmetadata = auxmetadata,
                    ##         Oout = 'numeric',
                    ##         Nout = 'numeric',
                    ##         Bout = 'numeric',
                    ##         logjacobianOr = NULL
                    ##     ))
                    Xgrid <- cbind(datavalues)
                    colnames(Xgrid) <- name
                    rownames(Xgrid) <- datavalues

                } else {
                    ## we must construct an X-grid
                    Xgrid <- cbind(seq(plotmin, plotmax, by = halfstep * 2))
                    colnames(Xgrid) <- name
                    rownames(Xgrid) <- NULL
                    Xticks <- NULL
                }

                probabilities <- Pr(Y = Xgrid, X = NULL,
                    learnt = learnt,
                    quantiles = quants,
                    nsamples = nmcsamples,
                    parallel = parallel,
                    silent = TRUE,
                    keepYX = FALSE)

                dim(probabilities$values) <- NULL


                ## Find appropriate plot height across plots
                if (plotvariability == 'samples') {
                    dim(probabilities$samples) <- dim(probabilities$samples)[-2]
                    ymax <- quantile(apply(probabilities$samples, 2,
                        function(x) {
                            quantile(x, probs = 31 / 32, type = 6,
                                na.rm = TRUE, names = FALSE)
                        }
                    ), probs = 31 / 32, type = 6, na.rm = TRUE, names = FALSE)
                    ## ymax <- quantile(apply(
                    ##     probabilities[, subsamples, drop = FALSE],
                    ##     2, function(x) {
                    ##         quantile(x, 31 / 32, type = 6, na.rm = TRUE)
                    ##     }
                    ## ), 31 / 32, type = 6, na.rm = TRUE)
                } else {
                    dim(probabilities$quantiles) <- dim(probabilities$quantiles)[-2]
                    temp <- probabilities$quantiles
                    ymax <- max(temp[is.finite(temp)])
                    ## ymax <- apply(probabilities[, , drop = FALSE], 1,
                    ##     function(x) {
                    ##         quantile(x, max(quants), type = 6, na.rm = TRUE)
                    ##     }
                    ## )
                    ## ymax <- max(ymax[is.finite(ymax)])
                }

                ## prepare info for histogram data plots if required
                if (datahistogram && theresdata) {
                    datum <- data[[name]]
                    datum <- datum[!is.na(datum)]

                    if(nvaluelist > 0) {
                        histo <- as.vector(table(factor(datum,
                            levels = rownames(Xgrid)))) / length(datum)
                    } else {
                        ## histo <- thist(datum,
                        ##     n = c(Xgrid - halfstep, max(Xgrid) + halfstep),
                        ##     extendbreaks = FALSE)$counts / length(datum)
                        histo <- hist(datum,
                            breaks = c(Xgrid - halfstep, max(Xgrid) + halfstep),
                            plot = FALSE)$counts / length(datum)
                    }

                    ymax <- max(ymax, histo)
                }

                ## If required, plot frequency samples
                if (plotvariability == 'samples') {
                    ## tplot(
                    ##     x = Xgrid,
                    ##     y = probabilities$samples,
                    ##     xlim = range(Xgrid), ylim = c(0, ymax),
                    ##     xticks = Xticks, xlabels = rownames(Xgrid),
                    ##     type = 'l', lty = 1, lwd = 2,
                    ##     col = 5, alpha = 1/8,
                    ##     xlab = name,
                    ##     ylab = paste0('probability', addylab),
                    ##     family = fontfamily
                    ## )
                    flexiplot(
                        x = Xgrid,
                        y = probabilities$samples,
                        ## xlim = range(Xgrid),
                        ylim = c(0, ymax),
                        ## xticks = Xticks, xlabels = rownames(Xgrid),
                        type = 'l', lty = 1, lwd = 2,
                        col = adjustcolor(cyan, 1/8),
                        xlab = name,
                        ylab = paste0('probability', addylab),
                        family = fontfamily
                    )
                    addplot <- TRUE

                } else if (plotvariability == 'quantiles') {
                    ## marguncertainty <- t(apply(probabilities, 1, function(x) {
                    ##     quantile(x, quants, type = 6, na.rm = TRUE)
                    ## }))
                    plotquantiles(
                        x = Xgrid,
                        y = probabilities$quantiles,
                        col = cyan, alpha.f = 0.25,
                        ## xlim = range(Xgrid),
                        ylim = c(0, ymax),
                        ## xticks = Xticks, xlabels = rownames(Xgrid),
                        xlab = name,
                        ylab = paste0('probability', addylab),
                        family = fontfamily,
                        add = addplot
                    )
                    addplot <- TRUE
                }

                ## If required, plot probability
                if (plotprobability) {
                    ## tplot(
                    ##     x = Xgrid,
                    ##     y = probabilities$values,
                    ##     xlim = range(Xgrid), ylim = c(0, ymax),
                    ##     xticks = Xticks, xlabels = rownames(Xgrid),
                    ##     type = 'b', cex = 0.5, lty = 1, lwd = 4,
                    ##     col = 1, alpha = 0.75,
                    ##     xlab = name,
                    ##     ylab = paste0('probability', addylab),
                    ##     family = fontfamily,
                    ##     add = addplot
                    ## )
                    flexiplot(
                        x = Xgrid,
                        y = probabilities$values,
                        xlim = range(Xgrid), ylim = c(0, ymax),
                        ## xticks = Xticks, xlabels = rownames(Xgrid),
                        type = 'b', pch = 16, cex = 0.5, lty = 1, lwd = 4,
                        col = adjustcolor(blue, 0.75),
                        xlab = name,
                        ylab = paste0('probability', addylab),
                        family = fontfamily,
                        add = addplot
                    )
                    addplot <- TRUE
                }

                ## If required and possible, plot data histogram
                if (datahistogram && theresdata) {
                    ## tplot(
                    ##     x = Xgrid, y = histo,
                    ##     xlim = range(Xgrid), ylim = c(0, ymax),
                    ##     xticks = Xticks, xlabels = rownames(Xgrid),
                    ##     type = 'b', col = 4, alpha = 0.5,
                    ##     border = '#555555', border.alpha = 1/4,
                    ##     xlab = name,
                    ##     ylab = paste0('probability', addylab),
                    ##     family = fontfamily,
                    ##     add = addplot
                    ## )
                    flexiplot(
                        x = Xgrid, y = histo,
                        xlim = range(Xgrid), ylim = c(0, ymax),
                        ## xticks = Xticks, xlabels = rownames(Xgrid),
                        type = 'b', pch = 1, cex = 0.5, lty = 1, lwd = 2,
                        col = adjustcolor(yellow, 0.5),
                        ## border = '#555555', border.alpha = 1/4,
                        xlab = name,
                        ylab = paste0('probability', addylab),
                        family = fontfamily,
                        add = addplot
                    )
                    addplot <- TRUE
                }
                ## End D, O, N, B case
            } else {
                message('Unknown MC type for variate ', name)
            }

            ## Scatter plot of data if required
            if (datascatter && theresdata) {
                datum <- data[[name]]
                datum <- datum[!is.na(datum)]
                if (mcmctype %in% c('O', 'N', 'B')) {
                    datum <- as.matrix(
                        vtransform(x = datum,
                            variates = name,
                            auxmetadata = auxmetadata,
                            Oout = 'numeric',
                            Nout = 'numeric',
                            Bout = 'numeric',
                            logjacobianOr = NULL
                        ))
                }
                rug(x = jitter(datum, amount = 0), side = 1,
                    col = adjustcolor(yellow,
                        alpha.f = exp((-length(datum) + 1)/128)),
                    quiet = TRUE)
                ## scatteraxis(
                ##     side = 1, n = NA, alpha = 0.75, ext = 5,
                ##     x = datum + runif(length(datum),
                ##         min = -min(diff(sort(c(par('usr')[1:2],
                ##             unique(datum))))) / 1.5,
                ##         max = min(diff(sort(c(par('usr')[1:2],
                ##             unique(datum))))) / 1.5
                ##     ),
                ##     col = 4
                ## )
                ## ## These lines plot quartiles
                ## if (vtype %in% c('R', 'D', 'C', 'L')) {
                ## fiven <- fivenum(datum)
                ## abline(v = fiven, col = paste0(palette()[c(2, 4, 5, 4, 2)], '44'),
                ##     lwd = 5, lty = 2)
                ## }
            }
        }) # End with
    }

    dev.off()
}
