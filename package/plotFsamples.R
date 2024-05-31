# DESCRIPTION
#' @param mcoutput
#' @param silent, Boolean # could set samplesDFistribution to TRUE
plotFsamples <- function(file, mcoutput, data, plotmeans = TRUE,
                         plotuncertainty = 'samples', uncertainty = 100,
                         datahistogram = TRUE, datascatter = TRUE,
                         useOquantiles = TRUE, parallel = TRUE,
                         silent = FALSE) {
  #imports
  source('tplotfunctions.R')
  source('vtransform.R')
  source('mcsubset.R')
  source('samplesFDistribution.R')
  
  family <- 'Palatino'

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
  if (plotuncertainty == 'quantiles') {
    plotmeans <- TRUE
    if (any(uncertainty <= 0 | uncertainty >= 1)) {
      uncertainty <- c(1, 7) / 8
    }
    quants <- sort(unique(round(c(uncertainty, 1 - uncertainty), 6)))
    mcsubsamples <- subsamples <- 1:nsamples
    addylab <- paste0(' (', ceiling(diff(quants) * 100), '% unc.)')
  } else {
    if (uncertainty == 'all') {
      uncertainty <- nsamples
    }
    uncertainty <- abs(uncertainty)
    if (plotmeans) {
      mcsubsamples <- 1:nsamples
    } else {
      mcsubsamples <- round(seq(1, nsamples, length.out = abs(uncertainty)))
    }
    subsamples <- round(seq(1, length(mcsubsamples), length.out = uncertainty))
  }

  addplot <- FALSE

  graphics.off()
  pdff(file, apaper = 4)
  par(mfrow = c(1, 1))

  for (v in auxmetadata[['name']]) {
    ## Check if we have proper data for this variate
    if (datahistogram || datascatter) {
      theresdata <- (sum(!is.na(data)) > 1)
    }

    varinfo <- as.list(auxmetadata[name == v])
    vtype <- varinfo[['mcmctype']]
    if (vtype %in% c('R', 'D', 'C', 'O')) { #
      if (vtype == 'O') {
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
                                          useOquantiles = useOquantiles,
                                          parallel = parallel,
                                          silent = TRUE)

      if (plotuncertainty == 'samples') {
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
        if (vtype == 'O') {
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
      if (plotuncertainty == 'samples') {
        tplot(
          x = Xgrid[xleft & xright],
          y = plotsamples[xleft & xright, subsamples, drop = FALSE],
          xlim = range(Xgrid), ylim = c(0, ymax),
          type = 'l', lty = 1, lwd = 2,
          col = 5, alpha = 7 / 8,
          xlab = v,
          ylab = paste0('frequency', (if (vtype == 'O') {
            ''
          } else {
            ' density'
          }), addylab),
          family = family
        )
        if (any(!(xleft & xright))) {
          tplot(
            x = Xgrid[!(xleft & xright)],
            y = plotsamples[!(xleft & xright), subsamples, drop = FALSE] * ymax,
            type = 'p', pch = 2, cex = 2,
            col = 5, alpha = 7 / 8,
            family = family, add = TRUE
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
          type = (if (vtype == 'O') {
            'b'
          } else {
            'l'
          }), cex = 0.5, lty = 1, lwd = 4,
          col = 1, alpha = 0.25,
          xlab = v,
          ylab = paste0('frequency', (if (vtype == 'O') {
            ''
          } else {
            ' density'
          }), addylab),
          family = family,
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

      if (plotuncertainty == 'quantiles') {
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
          type = (if (vtype == 'O') {
            'b'
          } else {
            'l'
          }), cex = 0.5, lty = 1, lwd = 2,
          col = yellow, alpha = 0.5, border = darkgrey, border.alpha = 3 / 4,
          xlab = v,
          ylab = paste0('frequency', (if (vtype == 'O') {
            ''
          } else {
            ' density'
          }), addylab),
          family = family, add = TRUE
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
            family = family, add = TRUE
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
        Nout = 'numeric', Bout = 'numeric', useOquantiles = useOquantiles
      )

      plotsamples <- samplesFDistribution(Y = Xgrid, X = NULL,
                                          mcoutput = mcoutput,
                                          subsamples = mcsubsamples,
                                          jacobian = TRUE,
                                          useOquantiles = useOquantiles,
                                          parallel = parallel,
                                          silent = TRUE)

      if (plotuncertainty == 'samples') {
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
      if (plotuncertainty == 'samples') {
        tplot(
          x = Ngrid, y = plotsamples[, subsamples, drop = FALSE],
          xlim = range(Ngrid), ylim = c(0, ymax),
          xticks = Ngrid, xlabels = Xgrid,
          type = 'l', lty = 1, lwd = 2,
          col = 5, alpha = 7 / 8,
          xlab = v,
          ylab = paste0('frequency', addylab),
          family = family
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
          family = family,
          add = addplot
        )
      }
      if (plotuncertainty == 'quantiles') {
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
          col = yellow, alpha = 0.5, border = darkgrey, border.alpha = 3 / 4,
          xlab = v,
          ylab = paste0('frequency', addylab),
          family = family, add = TRUE
        )
      }
    }

    if (datascatter && theresdata) {
      datum <- data[[v]]
      datum <- datum[!is.na(datum)]
      if (!(vtype %in% c('R', 'D', 'C', 'O'))) {
        datum <- vtransform(
          x = matrix(datum, ncol = 1, nrow = length(datum),
                     dimnames = list(NULL, v)), auxmetadata = auxmetadata,
          Nout = 'numeric', Bout = 'numeric', useOquantiles = useOquantiles
        )
      }
      scatteraxis(
        side = 1, n = NA, alpha = 0.75, ext = 5,
        x = datum + runif(length(datum),
          min = -min(diff(sort(c(par('usr')[1:2], unique(datum))))) / 1.5,
          max = min(diff(sort(c(par('usr')[1:2], unique(datum))))) / 1.5
        ),
        col = yellow
      )
      if (vtype %in% c('R', 'D', 'C', 'O')) {
        fiven <- fivenum(datum)
        abline(v = fiven, col = paste0(palette()[c(2, 4, 5, 4, 2)], '44'),
               lwd = 5, lty = 2)
      }
    }
  }
  dev.off()
}
