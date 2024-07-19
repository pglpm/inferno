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
  iqrfactor <- 2 * 2
  iqrfactorLog <- 2 * 2 # 6
  iqrfactorQ <- 2 * 0.25

  idR <- idC <- idD <- idL <- idB <- idO <- idN <- 1L

  auxmetadata <- data.frame()

  for (name in metadata$name) {
    if(!is.null(data)) {
      x <- data[[name]]
      x <- x[!is.na(x)]
    }
    xinfo <- as.list(metadata[metadata$name == name, ])
    ## make sure 'type' is lowercase
    xinfo$type <- tolower(xinfo$type)
    ordinal <- NA
    transf <- 'identity' # temporary
    datavalues <- xinfo[grep('^V[0-9]+$', names(xinfo))]
    if (xinfo$type == 'binary') {
      mcmctype <- 'B'
      id <- idB
      idB <- idB + 1L
      Nvalues <- xinfo$Nvalues
      step <- xinfo$rounding / 2
      domainmin <- 0
      domainmax <- 1
      censormin <- -Inf
      censormax <- +Inf
      tlocation <- 0
      tscale <- 1
      plotmin <- NA
      plotmax <- NA
      mctest1 <- 1
      mctest2 <- 2
      mctest3 <- 2
      ## if(!is.null(data)) {
      ##   mctest1 <- match(names(which.min(table(x))), datavalues)
      ##   mctest2 <- match(names(which.max(table(x))), datavalues)
      ##   mctest3 <- match(names(which.min(table(x))), datavalues)
      ## } else {
      ##   mctest1 <- 1
      ##   mctest2 <- 2
      ##   mctest3 <- 2
      ## }
    } else if (xinfo$type == 'ordinal') {
      ## nominal variate
      mcmctype <- 'O'
      id <- idO
      idO <- idO + 1L
      Nvalues <- xinfo$Nvalues
      step <- 0.5
      domainmin <- 1 # Nimble index categorical from 1
      domainmax <- Nvalues
      censormin <- -Inf
      censormax <- +Inf
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
            ## mctest2 <- round(xinfo$Nvalues/2)
            ## mctest3 <- xinfo$Nvalues
      ## }
            mctest1 <- 1
            mctest2 <- round(Nvalues/2)
            mctest3 <- Nvalues
    } else if (xinfo$type == 'nominal') {
      ## nominal variate
      mcmctype <- 'N'
      id <- idN
      idN <- idN + 1L
      Nvalues <- xinfo$Nvalues
      step <- 0.5
      domainmin <- 1 # Nimble index categorical from 1
      domainmax <- Nvalues
      censormin <- -Inf
      censormax <- +Inf
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
            ## mctest2 <- round(xinfo$Nvalues/2)
            ## mctest3 <- xinfo$Nvalues
      ## }
            mctest1 <- 1
            mctest2 <- round(Nvalues/2)
            mctest3 <- Nvalues
    } else if (xinfo$type == 'latent') {
      ## old treatment of ordinal variate
      ## to be deleted in the future, if not used
      mcmctype <- 'L'
      id <- idL
      idL <- idL + 1L
      transf <- 'Q'
      ordinal <- TRUE
      Nvalues <- xinfo$Nvalues
      step <- 0.5
      domainmin <- xinfo$domainmin
      domainmax <- xinfo$domainmax
      censormin <- -Inf
      censormax <- +Inf
      ## datavalues <- as.vector(xinfo[paste0('V',1:Nvalues)], mode='character')
      olocation <- (Nvalues * domainmin - domainmax) / (Nvalues - 1)
      oscale <- (domainmax - domainmin) / (Nvalues - 1)
      tlocation <- Qf(round((xinfo$centralvalue - olocation) / oscale) / Nvalues)
      tscale <- abs(
        Qf(round((xinfo$highvalue - olocation) / oscale) / Nvalues) -
        Qf(round((xinfo$lowvalue - olocation) / oscale) / Nvalues)
      ) / iqrfactorQ
      plotmin <- (if(is.finite(xinfo$plotmin)){xinfo$plotmin}else{xinfo$domainmin})
      plotmax <- (if(is.finite(xinfo$plotmax)){xinfo$plotmax}else{xinfo$domainmax})
      ## Q1 <- mctest1 <- quantile(x, probs = 0.25, type = 6)
      ## Q2 <- mctest2 <- quantile(x, probs = 0.5, type = 6)
      ## Q3 <- mctest3 <- quantile(x, probs = 0.75, type = 6)
      mctest1 <- xinfo$lowvalue
      mctest2 <- xinfo$centralvalue
      mctest3 <- xinfo$highvalue
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
    } else if (xinfo$type == 'continuous') { # continuous variate (R,C,D)
      Nvalues <- +Inf

      ## Rounded variates
      if (is.null(xinfo$rounding) || is.na(xinfo$rounding)) {
        step <- 0
      } else {
        step <- xinfo$rounding / 2
      }
      ## If the variate is rounded,
      ## no latent-variable representation is needed anyway
      ## if the datapoints are distinct in higher dimension
      if(step > 0 &&
         (is.null(data) || nrow(unique(data))/nrow(data) > 0.5)) {
        step <- 0
      }

      domainmin <- xinfo$domainmin
      domainmax <- xinfo$domainmax
      ## censormin <- max(domainmin, xinfo$minincluded, na.rm=TRUE)
      ## censormax <- min(domainmax, xinfo$maxincluded, na.rm=TRUE)
      ## closeddomain <- (censormin > domainmin) || (censormax < domainmax)
      plotmin <- xinfo$plotmin
      plotmax <- xinfo$plotmax
      ## Q1 <- mctest1 <- quantile(x, probs = 0.25, type = 6)
      ## Q2 <- mctest2 <- quantile(x, probs = 0.5, type = 6)
      ## Q3 <- mctest3 <- quantile(x, probs = 0.75, type = 6)
      mctest1 <- xinfo$lowvalue
      mctest2 <- xinfo$centralvalue
      mctest3 <- xinfo$highvalue

#### variate transformation to real-line domain
      if (is.finite(xinfo$domainmin) && is.finite(xinfo$domainmax) &&
          !xinfo$minincluded && !xinfo$maxincluded) {
        ## doubly-bounded open domain
        transf <- 'Q'
        tlocation <- Qf(0.5 +
                        (xinfo$centralvalue - (domainmax + domainmin)/2) /
                        (domainmax - domainmin))
        tscale <- abs(
          Qf(0.5 +
             (xinfo$highvalue - (domainmax + domainmin)/2) /
             (domainmax - domainmin)) -
          Qf(0.5 +
             (xinfo$lowvalue - (domainmax + domainmin)/2) /
             (domainmax - domainmin))
        ) / iqrfactorQ
        closeddomain <- FALSE
        censormin <- -Inf
        censormax <- +Inf
      } else if (is.finite(xinfo$domainmin) && is.finite(xinfo$domainmax) &&
                 xinfo$minincluded && !xinfo$maxincluded) {
        ## doubly-bounded left-closed domain
        transf <- 'logminus'
        tlocation <- log(domainmax - xinfo$centralvalue)
        tscale <- abs(
          log(domainmax - xinfo$highvalue) -
          log(domainmax - xinfo$lowvalue)
        ) / iqrfactorLog
        closeddomain <- TRUE
        censormin <- domainmin
        censormax <- domainmax
      } else if (is.finite(xinfo$domainmin) && is.finite(xinfo$domainmax) &&
                 !xinfo$minincluded && xinfo$maxincluded) {
        ## doubly-bounded right-closed domain
        transf <- 'log'
        tlocation <- log(xinfo$centralvalue - domainmin)
        tscale <- abs(
          log(xinfo$highvalue - domainmin) -
          log(xinfo$lowvalue - domainmin)
        ) / iqrfactorLog
        closeddomain <- TRUE
        censormin <- domainmin
        censormax <- domainmax
      } else if (is.finite(xinfo$domainmin) && is.finite(xinfo$domainmax) &&
                 xinfo$minincluded && xinfo$maxincluded) {
        ## doubly-bounded closed domain
        transf <- 'identity'
        tlocation <- xinfo$centralvalue
        tscale <- abs(xinfo$highvalue - xinfo$lowvalue) / iqrfactor
        closeddomain <- TRUE
        censormin <- domainmin
        censormax <- domainmax
      } else if (is.finite(xinfo$domainmin) && !xinfo$minincluded) {
        ## left-bounded left-open domain
        transf <- 'log'
        tlocation <- log(xinfo$centralvalue - domainmin)
        tscale <- abs(
          log(xinfo$highvalue - domainmin) -
          log(xinfo$lowvalue - domainmin)
        ) / iqrfactorLog
        closeddomain <- FALSE
        censormin <- domainmin
        censormax <- +Inf
      } else if (is.finite(xinfo$domainmin) && xinfo$minincluded) {
        ## left-bounded left-closed domain
        transf <- 'identity'
        tlocation <- xinfo$centralvalue
        tscale <- abs(xinfo$highvalue - xinfo$lowvalue) / iqrfactor
        closeddomain <- TRUE
        censormin <- domainmin
        censormax <- +Inf
        domainmin <- -Inf
      } else if (is.finite(xinfo$domainmax) && !xinfo$maxincluded) {
        ## right-bounded right-open domain
        transf <- 'logminus'
        tlocation <- log(domainmax - xinfo$centralvalue)
        tscale <- abs(
          log(domainmax - xinfo$highvalue) -
          log(domainmax - xinfo$lowvalue)
        ) / iqrfactorLog
        closeddomain <- FALSE
        censormin <- -Inf
        censormax <- domainmax
      } else if (is.finite(xinfo$domainmax) && xinfo$maxincluded) {
        ## right-bounded right-closed domain
        transf <- 'identity'
        tlocation <- xinfo$centralvalue
        tscale <- abs(xinfo$highvalue - xinfo$lowvalue) / iqrfactor
        closeddomain <- TRUE
        censormin <- -Inf
        censormax <- domainmax
        domainmax <- +Inf
      } else {
        ## unbounded, open domain
        transf <- 'identity'
        tlocation <- xinfo$centralvalue
        tscale <- abs(xinfo$highvalue - xinfo$lowvalue) / iqrfactor
        closeddomain <- FALSE
        censormin <- -Inf
        censormax <- +Inf
      }
    ## # Old cases
    ## if (is.finite(xinfo$domainmin) && is.finite(xinfo$domainmax)) {
    ##   transf <- 'Q'
    ##   if (xinfo$minincluded && !xinfo$maxincluded) {
    ##     censormin <- domainmin
    ##     censormax <- +Inf
    ##     domainmin <- (8 * domainmin - domainmax) / 7
    ##   } else if (!xinfo$minincluded && xinfo$maxincluded) {
    ##     censormax <- domainmax
    ##     censormin <- -Inf
    ##     domainmax <- (8 * domainmax - domainmin) / 7
    ##   } else if (xinfo$minincluded && xinfo$maxincluded) {
    ##     censormin <- domainmin
    ##     censormax <- domainmax
    ##     domainmin <- (7 * domainmin - domainmax) / 6
    ##     domainmax <- (7 * domainmax - domainmin) / 6
    ##   }
    ##   tlocation <- Qf((tlocation - domainmin) / (domainmax - domainmin))
    ##   tscale <- abs(Qf((xinfo$highvalue - domainmin) / (domainmax - domainmin)) -
    ##                  Qf((xinfo$lowvalue - domainmin) / (domainmax - domainmin))) *
    ##     sdoveriqr
    ## } else if (is.finite(xinfo$domainmin)) {
    ##   transf <- 'log'
    ##   tlocation <- log(tlocation - domainmin)
    ##   tscale <- abs(log(xinfo$highvalue - domainmin) -
    ##                  log(xinfo$lowvalue - domainmin)) * sdoveriqr
    ## } else if (is.finite(xinfo$domainmax)) {
    ##   transf <- 'logminus'
    ##   tlocation <- log(domainmax - tlocation)
    ##   tscale <- abs(log(domainmax - xinfo$highvalue) -
    ##                log(domainmax - xinfo$lowvalue)) * sdoveriqr
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
          tlocation = tlocation, tscale = tscale,
          plotmin = plotmin, plotmax = plotmax,
          ## Q1 = Q1, Q2 = Q2, Q3 = Q3, # not used in other scripts, possibly remove
          mctest1 = mctest1, mctest2 = mctest2, mctest3 = mctest3
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
