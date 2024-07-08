#' Plot one-dimensional posterior probabilities
#'
#' @param file string: name of plot output file
#' @param mcoutput Either a string with the name of a directory or full
#'   path for a 'FDistribution.rds' object, or such an object itself
#' @param data data.table object or filepath: datapoints
#' @param plotvariability string, either 'samples' or 'quantiles': how to plot the variability of the probability distribution with new samples
#' @param nFsamples positive number: if plotvariability='samples', then number of samples of representative frequency distributions to display as variability; if plotvariability='quantiles', then the quantiles (in range 0 to 0.5) to show
#' @param datahistogram bool: plot the data as histogram?
#' @param datascatter bool: plot the data as scatterplot along the x-axis?
#' @param parallel Bool or numeric: whether to use pre-existing parallel
#'   workers, or how many to create and use
#' @param useLquantiles Bool: internal, use metadata quantiles for ordinal
#'   variates
#' @param silent Bool: give warnings or updates in the computation
#'
#' @return A list with the mutual information, its error, and its unit
#' @export
plotFsamples <- function(file, mcoutput, data, plotmeans = TRUE,
                         plotvariability = 'samples', nFsamples = 100,
                         datahistogram = TRUE, datascatter = TRUE,
                         useLquantiles = FALSE, parallel = TRUE,
                         silent = FALSE) {

  ## old utility functions
  ## source('tplotfunctions.R')
  ## source('vtransform.R')
  ## source('mcsubset.R')
  ## source('samplesFDistribution.R')

  fontfamily <- 'Palatino'

  ## Extract Monte Carlo output & aux-metadata
  # This function does not make sense with the input format?
  if (is.character(mcoutput)) {
    if (file_test('-d', mcoutput) &&
          file.exists(paste0(mcoutput, '/Fdistribution.rds'))) {
      mcoutput <- readRDS(paste0(mcoutput, '/Fdistribution.rds'))
    } else {
      mcoutput <- paste0(sub('.rds$', '', mcoutput), '.rds')
      if (file.exists(mcoutput)) {
        mcoutput <- readRDS(mcoutput, '/Fdistribution.rds')
      } else {
        stop('cannot find mcoutput file')
      }
    }
  }
  auxmetadata <- mcoutput$auxmetadata
  ## mcoutput$auxmetadata <- NULL
  nsamples <- ncol(mcoutput$W)

  nodata <- missing(data) || is.null(data) || (is.logical(data) && !data)
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
    plotmeans <- TRUE
    if (any(nFsamples <= 0 | nFsamples >= 1)) {
      nFsamples <- c(1, 7) / 8
    }
    quants <- sort(unique(round(c(nFsamples, 1 - nFsamples), 6)))
    mcsubsamples <- subsamples <- 1:nsamples
    addylab <- paste0(' (', ceiling(diff(quants) * 100), '% unc.)')
  } else {
    if (nsamples == 'all') {
      nFsamples <- nsamples
    }
    nFsamples <- abs(nFsamples)
    if (plotmeans) {
      mcsubsamples <- 1:nsamples
    } else {
      mcsubsamples <- round(seq(1, nsamples, length.out = abs(nFsamples)))
    }
    subsamples <- round(seq(1, length(mcsubsamples), length.out = nFsamples))
  }

  addplot <- FALSE

  graphics.off()
  pdff(file, apaper = 4)
  par(mfrow = c(1, 1))

  for (v in auxmetadata[['name']]) {
    ## Check if we have proper data for this variate
    if (datahistogram || datascatter) {
      theresdata <- (sum(!is.na(data[[v]])) > 1)
    }

    varinfo <- as.list(auxmetadata[auxmetadata$name == v, ])
    varinfo$Nvalues <- abs(varinfo$Nvalues)
    vtype <- varinfo[['mcmctype']]
    if (vtype %in% c('R', 'D', 'C', 'L')) { #
      if (vtype == 'L') {
        Xgrid <- seq(varinfo[['domainmin']], varinfo[['domainmax']],
          length.out = varinfo[['Nvalues']]
        )
        Xgrid <- cbind(Xgrid[Xgrid >= varinfo[['plotmin']]
                             & Xgrid <= varinfo[['plotmax']]])
      } else {
        Xgrid <- cbind(seq(
          varinfo[['plotmin']], varinfo[['plotmax']],
          length.out = 256
        ))
      }
      colnames(Xgrid) <- v
      xleft <- Xgrid > varinfo[['censormin']]
      xright <- Xgrid < varinfo[['censormax']]

      plotsamples <- samplesFDistribution(Y = Xgrid, X = NULL,
                                          mcoutput = mcoutput,
                                          subsamples = mcsubsamples,
                                          jacobian = TRUE,
                                          useLquantiles = useLquantiles,
                                          parallel = parallel,
                                          silent = TRUE)

      if (plotvariability == 'samples') {
        ymax <- tquant(apply(
          plotsamples[xleft & xright, subsamples, drop = FALSE],
          2, function(x) {
            tquant(x, 31 / 32)
          }
        ), 31 / 32, na.rm = TRUE)
      } else {
        ymax <- apply(
          plotsamples[xleft & xright, , drop = FALSE], 1,
          function(x) {
            tquant(x, max(quants))
          }
        )
        ymax <- max(ymax[is.finite(ymax)])
      }

      ## data plots if required
      if (datahistogram && theresdata) {
        datum <- data[[v]]
        datum <- datum[!is.na(datum)]
        dleft <- datum > varinfo[['censormin']]
        dright <- datum < varinfo[['censormax']]
        if (vtype == 'L') {
          dh <- (varinfo[['domainmax']] - varinfo[['domainmin']]) /
            (varinfo[['Nvalues']] - 1L) / 2
          nh <- seq(varinfo[['domainmin']] - dh, varinfo[['domainmax']] + dh,
            length.out = varinfo[['Nvalues']] + 1L
          )
          nh <- nh[nh >= min(datum) - dh & nh <= max(datum) + dh]
        } else {
          ## nh <- seq(min(datum[dleft & dright]), max(datum[dleft & dright]),
          ## length.out=max(16, round(length(datum[dleft & dright])/64)))
          nh <- max(32, round(length(datum[dleft & dright]) / 64))
        }

        histo <- thist(datum[dleft & dright],
          n = nh, extendbreaks = FALSE
        )
        hleft <- sum(!dleft) / length(datum)
        hright <- sum(!dright) / length(datum)

        ymax <- max(ymax, histo$density)
      }

      ## Plot frequency samples
      if (plotvariability == 'samples') {
        tplot(
          x = Xgrid[xleft & xright],
          y = plotsamples[xleft & xright, subsamples, drop = FALSE],
          xlim = range(Xgrid), ylim = c(0, ymax),
          type = 'l', lty = 1, lwd = 2,
          col = 5, alpha = 7 / 8,
          xlab = v,
          ylab = paste0('frequency', (if (vtype == 'L') {
            ''
          } else {
            ' density'
          }), addylab),
          family = fontfamily
        )
        if (any(!(xleft & xright))) {
          tplot(
            x = Xgrid[!(xleft & xright)],
            y = plotsamples[!(xleft & xright), subsamples, drop = FALSE] * ymax,
            type = 'p', pch = 2, cex = 2,
            col = 5, alpha = 7 / 8,
            family = fontfamily, add = TRUE
          )
        }
        addplot <- TRUE
      }

      ## Plot FALSE means if required
      if (plotmeans) {
        tplot(
          x = Xgrid[xleft & xright],
          y = rowMeans(plotsamples[xleft & xright, , drop = FALSE],
                       na.rm = TRUE),
          xlim = range(Xgrid), ylim = c(0, ymax),
          type = (if (vtype == 'L') {
            'b'
          } else {
            'l'
          }), cex = 0.5, lty = 1, lwd = 4,
          col = 1, alpha = 0.25,
          xlab = v,
          ylab = paste0('frequency', (if (vtype == 'L') {
            ''
          } else {
            ' density'
          }), addylab),
          family = fontfamily,
          add = addplot
        )
        if (any(!(xleft & xright))) {
          tplot(
            x = Xgrid[!(xleft & xright)],
            y = rowMeans(plotsamples, na.rm = TRUE)[!(xleft & xright)] * ymax,
            type = 'p', pch = 2, cex = 2,
            col = 1, alpha = 0.25,
            lty = 1, lwd = 3,
            add = TRUE
          )
        }
      }

      if (plotvariability == 'quantiles') {
        marguncertainty <- t(apply(plotsamples, 1, function(x) {
          tquant(x, quants)
        }))
        plotquantiles(x = Xgrid[xleft & xright],
                      y = marguncertainty[xleft & xright, , drop = FALSE],
                      col = 5, alpha = 0.75)
        if (any(!(xleft & xright))) {
          tplot(
            x = matrix(Xgrid[!(xleft & xright)],
                       nrow = 2, ncol = sum(!(xleft & xright)), byrow = TRUE),
            y = t(marguncertainty[!(xleft & xright), , drop = FALSE]) * ymax,
            type = 'l', pch = 2, cex = 2,
            col = 5, alpha = 0.75,
            lty = 1, lwd = 16,
            add = TRUE
          )
        }
      }

      if (datahistogram && theresdata) {
        histomax <- 1 # max(rowMeans(plotsamples))/max(histo$density)
        tplot(
          x = histo$mids, y = histo$density * histomax,
          xlim = range(Xgrid), ylim = c(0, ymax),
          type = (if (vtype == 'L') {
            'b'
          } else {
            'l'
          }), cex = 0.5, lty = 1, lwd = 2,
          col = 4, alpha = 0.5, border = '#555555', border.alpha = 3 / 4,
          xlab = v,
          ylab = paste0('frequency', (if (vtype == 'L') {
            ''
          } else {
            ' density'
          }), addylab),
          family = fontfamily, add = TRUE
        )

        if (any(!(dleft & dright))) {
          tplot(
            x = c(if (hleft > 0) {
              varinfo[['censormin']]
            }, if (hright > 0) {
              varinfo[['censormax']]
            }), y = c(if (hleft > 0) {
              hleft
            }, if (hright > 0) {
              hright
            }) * ymax,
            type = 'p', pch = 0, cex = 2,
            col = 4, alpha = 0.5,
            lty = 1, lwd = 5,
            family = fontfamily, add = TRUE
          )
        }
        ## fiven <- fivenum(datum)
        ## abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'),lwd=5,lty=2)
      }

      #####
      ## nominal or binary variate
    } else {
      Xgrid <- cbind(unlist(varinfo[paste0('V', 1:varinfo[['Nvalues']])]))
      colnames(Xgrid) <- v
      Ngrid <- vtransform(
        x = Xgrid, auxmetadata = auxmetadata,
        Nout = 'numeric', Bout = 'numeric', useLquantiles = useLquantiles
      )

      plotsamples <- samplesFDistribution(Y = Xgrid, X = NULL,
                                          mcoutput = mcoutput,
                                          subsamples = mcsubsamples,
                                          jacobian = TRUE,
                                          useLquantiles = useLquantiles,
                                          parallel = parallel,
                                          silent = TRUE)

      if (plotvariability == 'samples') {
        ymax <- tquant(apply(
          plotsamples[, subsamples, drop = FALSE],
          2, function(x) {
            tquant(x, 31 / 32)
          }
        ), 31 / 32, na.rm = TRUE)
      } else {
        ymax <- apply(
          plotsamples[, , drop = FALSE],
          1, function(x) {
            tquant(x, max(quants))
          }
        )
        ymax <- max(ymax[is.finite(ymax)])
      }

      ## data plots if required
      if (datahistogram && theresdata) {
        datum <- data[[v]]
        datum <- datum[!is.na(datum)]
        histo <- as.vector(table(factor(datum, levels = Xgrid))) / length(datum)

        ymax <- max(ymax, histo)
      }

      ## Plot FALSE samples
      if (plotvariability == 'samples') {
        tplot(
          x = Ngrid, y = plotsamples[, subsamples, drop = FALSE],
          xlim = range(Ngrid), ylim = c(0, ymax),
          xticks = Ngrid, xlabels = Xgrid,
          type = 'l', lty = 1, lwd = 2,
          col = 5, alpha = 7 / 8,
          xlab = v,
          ylab = paste0('frequency', addylab),
          family = fontfamily
        )
        addplot <- TRUE
      }
      ## Plot FALSE means if required
      if (plotmeans) {
        tplot(
          x = Ngrid, y = rowMeans(plotsamples, na.rm = TRUE),
          xlim = range(Ngrid), ylim = c(0, ymax),
          xticks = Ngrid, xlabels = Xgrid,
          type = 'b', cex = 0.5, lty = 1, lwd = 4,
          col = 1, alpha = 0.25,
          xlab = v,
          ylab = paste0('frequency', addylab),
          family = fontfamily,
          add = addplot
        )
      }
      if (plotvariability == 'quantiles') {
        marguncertainty <- t(apply(plotsamples, 1, function(x) {
          tquant(x, quants)
        }))
        plotquantiles(x = Ngrid, y = marguncertainty, col = 5, alpha = 0.75)
      }

      if (datahistogram && theresdata) {
        histomax <- 1 # max(rowMeans(plotsamples))/max(histo$density)
        tplot(
          x = Ngrid, y = histo * histomax,
          xlim = range(Ngrid), ylim = c(0, ymax),
          xticks = Ngrid, xlabels = Xgrid,
          type = 'b', lty = 1, lwd = 2,
          col = 4, alpha = 0.5, border = '#555555', border.alpha = 3 / 4,
          xlab = v,
          ylab = paste0('frequency', addylab),
          family = fontfamily, add = TRUE
        )
      }
    }

    if (datascatter && theresdata) {
      datum <- data[[v]]
      datum <- datum[!is.na(datum)]
      if (!(vtype %in% c('R', 'D', 'C', 'L'))) {
        datum <- vtransform(
          x = matrix(datum, ncol = 1, nrow = length(datum),
                     dimnames = list(NULL, v)), auxmetadata = auxmetadata,
          Nout = 'numeric', Bout = 'numeric', useLquantiles = useLquantiles
        )
      }
      scatteraxis(
        side = 1, n = NA, alpha = 0.75, ext = 5,
        x = datum + runif(length(datum),
          min = -min(diff(sort(c(par('usr')[1:2], unique(datum))))) / 1.5,
          max = min(diff(sort(c(par('usr')[1:2], unique(datum))))) / 1.5
        ),
        col = 4
      )
      if (vtype %in% c('R', 'D', 'C', 'L')) {
        fiven <- fivenum(datum)
        abline(v = fiven, col = paste0(palette()[c(2, 4, 5, 4, 2)], '44'),
               lwd = 5, lty = 2)
      }
    }
  }
  dev.off()
}
