#' Build preliminary metadata flie
#'
#' @param data data.frame object or filepath
#' @param file string: name of output metadata file; NULL: output metadata as data.frame
#' @param diagnosticvalues Bool: also output some diagnostic statistics?
#' @param backupfiles Bool: rename previous metadata file if it exists?
#'
#' @return nothing or data.table object
#' @export
buildmetadata <- function(data, file = NULL,
                          diagnosticvalues = FALSE,
                          backupfiles = FALSE) {
  ## require('data.table')

  gcd <- function(...) {
    suppressWarnings(Reduce(function(a, b) {
      if (b == 0) a else Recall(b, a %% b)
    }, c(...)))
  }
  ##
  datafile <- NULL
  # Read datafile if it exists
  if (is.character(data) && file.exists(data)) {
    datafile <- data
    data <- read.csv(datafile, na.strings = '')
  }
  data <- as.data.frame(data)
  metadata <- as.data.frame(c(list(name = NA, type = NA),
        if (diagnosticvalues) {
          list(datamin = NA, datamax = NA, datamaxrep = NA,
               dataNvalues = NA)
        },
        list(Nvalues = NA, rounding = NA, domainmin = NA,
             domainmax = NA, minincluded = NA,
             maxincluded = NA, centralvalue = NA, lowvalue = NA,
             highvalue = NA, plotmin = NA, plotmax = NA)
        ))[-1,]
  # Loop over columns in data
  for (xn in colnames(data)) {
    ## print(xn)
    x <- data[[xn]]
    x <- x[!is.na(x)]
    # Checks if the data contains useful information
    if (all(is.na(x)) || min(x) == max(x)) {
      message('\nWARNING: variate ', xn,
              ' does not have at least two distinct values!',
              '\nDiscarded because non-informative\n')
      next
    }
    # If the data is numeric, calculate the variation in the data
    # (first and second quartile, min and max values)
    if (is.numeric(x)) {
      ## print('numeric')
      Q1 <- loval <- quantile(x, probs = 0.25, type = 6)
      Q2 <- meval <- quantile(x, probs = 0.5, type = 6) # Median. Not used for anything
      Q3 <- hival <- quantile(x, probs = 0.75, type = 6)
      dmin <- min(x)
      dmax <- max(x)
      # Edge case if the first and second quartile have the same value
      if (loval == hival) {
        loval <- if (sum(x < Q1) > 0) {
          max(x[x < Q1])
        } else {
          max(x[x <= Q1])
        }
        hival <- if (sum(x > Q3) > 0) {
          min(x[x > Q3])
        } else {
          min(x[x >= Q3])
        }
      }
    } else {
      loval <- meval <- hival <- dmin <- dmax <- NA
    }
    # Making an educated guess for the type of variates
    # Binary variate (only has two unique values)
    if (length(unique(x)) == 2) {
      vtype <- 'binary'
      vn <- 2
      vd <- NA
      domainmin <- NA
      domainmax <- NA
      censormin <- TRUE
      censormax <- TRUE
      vval <- sort(as.character(unique(x)))
      names(vval) <- paste0('V', 1:2)
      loval <- meval <- hival <- NA
      plotmin <- NA
      plotmax <- NA
      # Nominal variate (non numeric values)
    } else if (!is.numeric(x)) {
      vtype <- 'nominal'
      vn <- length(unique(x))
      vd <- NA
      domainmin <- NA # Nimble index categorical from 1
      domainmax <- NA
      censormin <- TRUE
      censormax <- TRUE
      vval <- sort(as.character(unique(x)))
      names(vval) <- paste0('V', 1:vn)
      plotmin <- NA
      plotmax <- NA
    } else if (length(unique(x)) <= 10) {
      vtype <- 'ordinal'
      vn <- length(unique(x))
      vd <- NA
      domainmin <- NA # Nimble index categorical from 1
      domainmax <- NA
      censormin <- TRUE
      censormax <- TRUE
      vval <- sort(as.character(unique(x)))
      names(vval) <- paste0('V', 1:vn)
      plotmin <- NA
      plotmax <- NA
    } else {
      ## Ordinal, continuous, censored, or discretized variate
      ## (numeric with more than 2 different values)
      ix <- x[!(x %in% range(x))] # exclude boundary values
      maxrep <- max(table(ix)) # average of repeated inner values
      ud <- unique(signif(diff(sort(unique(x))), 3)) # differences
      rx <- diff(range(x))
      multi <- 10^(-min(floor(log10(ud))))
      max(table(x))
      dd <- gcd(round(ud * multi)) / multi # greatest common difference
      ##
      if (dd / rx < 1e-5 && maxrep > 100) {
        cat('Warning: variate', xn,
            'seems continuous but with singular values.\n')
      }
      # Continuous variate (the difference between values can be
      # very small, but there is variation)
      if (dd / rx < 1e-5 || (dd / rx < 1e-3 && maxrep <= 100)) {
        ## temporary values
        vtype <- 'continuous'
        vn <- Inf
        vd <- 0
        domainmin <- signif(min(x) - 3 * diff(range(x)), 1)
        domainmax <- signif(max(x) + 3 * diff(range(x)), 1)
        censormin <- FALSE
        censormax <- FALSE
        plotmin <- min(x) - (Q3 - Q1) / 2
        plotmax <- max(x) + (Q3 - Q1) / 2
        ##
        if (all(x >= 0)) { # seems to be strictly positive
          domainmin <- 0
          plotmin <- max((domainmin + min(x)) / 2, plotmin)
        }
        ix <- x[!(x %in% range(x))] # exclude boundary values
        repindex <- mean(table(ix)) # average of repeated inner values
        ## contindex <- length(unique(diff(sort(unique(ix)))))/length(ix) # check for repeated values
        if (sum(x == min(x)) > repindex) { # seems to be left-singular
          domainmin <- plotmin <- min(x)
          censormin <- TRUE
        }
        if (sum(x == max(x)) > repindex) { # seems to be right-singular
          domainmax <- plotmax <- max(x)
          censormax <- TRUE
        }
      } else {
        ## # Seems to be an ordinal variate
        ## vtype <- 'ordinal'
        ## if (dd >= 1) { # seems originally integer
        ##   domainmin <- min(1, x)
        ##   domainmax <- max(x)
        ##   vn <- domainmax - domainmin + 1
        ##   vd <- NA
        ##   censormin <- TRUE
        ##   censormax <- TRUE
        ##   ## location <- NA # (vn*domainmin-domainmax)/(vn-1)
        ##   ## scale <- NA # (domainmax-domainmin)/(vn-1)
        ##   plotmin <- max(domainmin, min(x) - (Q3 - Q1) / 2)
        ##   plotmax <- max(x)
        ## } else { # seems a rounded continuous variate
          vtype <- 'continuous'
          vn <- Inf
          vd <- dd
          if (diff(range(x)) / vd > 256) {
            message('\nNOTE: variate ', xn, ' is reported as "rounded",\n',
                    'but consider the possibility of treating it as',
                    'continuous, \n by setting its "rounding" to 0 in the',
                    ' metadata file.\n')
          }
          domainmin <- signif(min(x) - 4 * diff(range(x)), 1)
          domainmax <- signif(max(x) + 4 * diff(range(x)), 1)
          censormin <- FALSE
          censormax <- FALSE
          ## location <- Q2
          ## scale <- (Q3-Q1)/2
          plotmin <- min(x) - (Q3 - Q1) / 2
          plotmax <- max(x) + (Q3 - Q1) / 2
          ix <- x[!(x %in% range(x))] # exclude boundary values
          repindex <- mean(table(ix)) # average of repeated inner values
          ## contindex <- length(unique(diff(sort(unique(ix)))))/length(ix) # check for repeated values
          if (all(x >= 0)) { # seems to be strictly positive
            domainmin <- 0
            plotmin <- max((domainmin + min(x)) / 2, plotmin)
          }
          ix <- x[!(x %in% range(x))] # exclude boundary values
          repindex <- mean(table(ix)) # average of repeated inner values
          ## contindex <- length(unique(diff(sort(unique(ix)))))/length(ix) # check for repeated values
          if (sum(x == min(x)) > repindex) { # seems to be left-singular
            domainmin <- plotmin <- min(x)
            censormin <- TRUE
          }
          if (sum(x == max(x)) > repindex) { # seems to be right-singular
            domainmax <- plotmax <- max(x)
            censormax <- TRUE
          }
        ## } # end rounded
      } # end integer
      vval <- NULL
    } # end numeric
    ##
    ## Create metadata object
    cat('\n***test***\n')
    str(metadata)
    metadata <- merge(metadata,
      c(
        list(name = xn, type = vtype),
        if (diagnosticvalues) {
          list(datamin = dmin, datamax = dmax, datamaxrep = max(table(x)),
               dataNvalues = length(unique(x)))
        },
        list(Nvalues = vn, rounding = vd, domainmin = domainmin,
             domainmax = domainmax, minincluded = censormin,
             maxincluded = censormax, centralvalue = meval, lowvalue = loval,
             highvalue = hival, plotmin = plotmin, plotmax = plotmax),
        as.list(vval)
      ),
      sort = FALSE, all = TRUE)
    cat('\n***test2***\n')
    str(metadata)
  } # End loop over columns
  ## metadata <- cbind(name=names(data), metadata)

  # Save to file if the file parameter is set
  if (!is.null(file)) {
    # If file is a string, name the metadata file that
    if (is.character(file)) {
      file <- paste0(sub('.csv$', '', file), '.csv')
    }
    # If file is not a string, and the data is read from a file,
    # save the metadata as metadata_<datafile>.csv
    else if (!is.null(datafile)) {
      file <- paste0(sub('.csv$', '', file), '_metadata.csv')
    }
    # If the file already exists, rename the old file as a backup file
    # and name this metadata file as intended
    if (backupfiles) {
      if (file.exists(file)) {
        file.rename(from = file, to = paste0(sub('.csv$', '', file), '_bak', format(Sys.time(), '%y%m%dT%H%M%S'),'.csv'))
      }
    }
    # Save the file
    write.csv(metadata, file, row.names = FALSE, quote = FALSE, na = '')
    cat('Saved proposal metadata file as', file, '\n')

  } else {
    # Else just print to console
    metadata
  }
}
