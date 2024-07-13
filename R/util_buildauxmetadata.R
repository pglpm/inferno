#' Build preliminary metadata flie
#'
#' @param data data.frame object
#' @param metadata data.frame object
#'
#' @return an auxmetadata data.frame object
buildauxmetadata <- function(data, metadata) {
  ## Elements used in other scripts:
  ## name, mcmctype, Nvalues, mctest, transform,
  ## domainmin, domainmax, tscale, tlocation
  ## censormax, censormin, step,
  ## mctest1, mctest2, mctest3,
  ## plotmin, plotmax, V...

  sdoveriqr <- 0.5 / qnorm(0.75)

  ## Qf <- readRDS('Qfunction3600_3.rds')

  ## Q <- readRDS('Qfunction512.rds')
  ##
  idR <- idC <- idD <- idL <- idB <- idO <- idN <- 1L

  auxmetadata <- data.frame()

  for (xn in metadata$name) {
    if(!is.null(data)) {
      x <- data[[xn]]
      x <- x[!is.na(x)]
    }
    xinfo <- as.list(metadata[metadata$name == xn, ])
    ## make sure 'type' is lowercase
    xinfo$type <- tolower(xinfo$type)
    ordinal <- NA
    cens <- any(c(xinfo$minincluded, xinfo$maxincluded), na.rm = TRUE)
    rounded <- NA
    transf <- 'identity' # temporary
    vval <- xinfo[grep('^V[0-9]+$', names(xinfo))]
    ## print(xn)
    ## str(vval)
    ## Q1 <- NA
    ## Q2 <- NA
    ## Q3 <- NA
    if (xinfo$type == 'binary') { # seems binary variate
      ## if (length(unique(x)) != 2) {
      ##   cat('Warning: inconsistencies with variate', xn, '\n')
      ## }
      vtype <- 'B'
      vid <- idB
      idB <- idB + 1L
      vn <- xinfo$Nvalues
      vd <- xinfo$rounding / 2
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
      ##   mctest1 <- match(names(which.min(table(x))), vval)
      ##   mctest2 <- match(names(which.max(table(x))), vval)
      ##   mctest3 <- match(names(which.min(table(x))), vval)
      ## } else {
      ##   mctest1 <- 1
      ##   mctest2 <- 2
      ##   mctest3 <- 2
      ## }
    } else if (xinfo$type == 'ordinal') {
      ## nominal variate
      vtype <- 'O'
      vid <- idO
      idO <- idO + 1L
      vn <- xinfo$Nvalues
      vd <- 0.5
      domainmin <- 1 # Nimble index categorical from 1
      domainmax <- vn
      censormin <- -Inf
      censormax <- +Inf
      tlocation <- 0
      tscale <- 1
      plotmin <- NA
      plotmax <- NA
      ## if(!is.null(data)) {
      ##   mctest1 <- match(names(which.min(table(x))), vval)
      ##   mctest2 <- match(names(which.max(table(x))), vval)
      ##   mctest3 <- match(names(which.min(table(x))), vval)
      ## } else {
            ## mctest1 <- 1
            ## mctest2 <- round(xinfo$Nvalues/2)
            ## mctest3 <- xinfo$Nvalues
      ## }
            mctest1 <- 1
            mctest2 <- round(vn/2)
            mctest3 <- vn
    } else if (xinfo$type == 'nominal') {
      ## nominal variate
      vtype <- 'N'
      vid <- idN
      idN <- idN + 1L
      vn <- xinfo$Nvalues
      vd <- 0.5
      domainmin <- 1 # Nimble index categorical from 1
      domainmax <- vn
      censormin <- -Inf
      censormax <- +Inf
      tlocation <- 0
      tscale <- 1
      plotmin <- NA
      plotmax <- NA
      ## if(!is.null(data)) {
      ##   mctest1 <- match(names(which.min(table(x))), vval)
      ##   mctest2 <- match(names(which.max(table(x))), vval)
      ##   mctest3 <- match(names(which.min(table(x))), vval)
      ## } else {
            ## mctest1 <- 1
            ## mctest2 <- round(xinfo$Nvalues/2)
            ## mctest3 <- xinfo$Nvalues
      ## }
            mctest1 <- 1
            mctest2 <- round(vn/2)
            mctest3 <- vn
    } else if (xinfo$type == 'latent') {
      ## old treatment of ordinal variate
      ## to be deleted in the future, if not used
      vtype <- 'L'
      vid <- idL
      idL <- idL + 1L
      transf <- 'Q'
      ordinal <- TRUE
      vn <- xinfo$Nvalues
      vd <- 0.5
      domainmin <- xinfo$domainmin
      domainmax <- xinfo$domainmax
      if ((domainmax - domainmin) != vn - 1) {
        cat('Warning: "Nvalues" for variate', xn, 'could be incorrect.\n')
      }
      censormin <- -Inf
      censormax <- +Inf
      ## vval <- as.vector(xinfo[paste0('V',1:vn)], mode='character')
      olocation <- (vn * domainmin - domainmax) / (vn - 1)
      oscale <- (domainmax - domainmin) / (vn - 1)
      tlocation <- Qf(round((xinfo$centralvalue - olocation) / oscale) / vn)
      tscale <- abs(Qf(round((xinfo$highvalue - olocation) / oscale) / vn) -
        Qf(round((xinfo$lowvalue - olocation) / oscale) / vn)) * sdoveriqr
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
      vn <- +Inf
      if (is.null(xinfo$rounding) || is.na(xinfo$rounding)) {
        xinfo$rounding <- 0
      }
      vd <- xinfo$rounding / 2
      rounded <- (vd > 0)
      domainmin <- censormin <- xinfo$domainmin
      domainmax <- censormax <- xinfo$domainmax
      ## censormin <- max(domainmin, xinfo$minincluded, na.rm=TRUE)
      ## censormax <- min(domainmax, xinfo$maxincluded, na.rm=TRUE)
      ## cens <- (censormin > domainmin) || (censormax < domainmax)
      tlocation <- xinfo$centralvalue
      ## with the sd of the means equal to 3 and the rate of the vrns equal
      ## to 1, the median of Q3 of the possible F distribs is 1.97
      ## so we set the scale to
      ## ([(Q3-Q2)+(Q2-Q1)]/2) /2
      ## this way the transformed Q3 is at approx 2
      tscale <- abs(xinfo$highvalue - xinfo$lowvalue)/4
      plotmin <- xinfo$plotmin
      plotmax <- xinfo$plotmax
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
      # needs transformation
      if (is.finite(xinfo$domainmin) && is.finite(xinfo$domainmax)) {
        transf <- 'Q'
        if (xinfo$minincluded && !xinfo$maxincluded) {
          censormin <- domainmin
          censormax <- +Inf
          domainmin <- (8 * domainmin - domainmax) / 7
        } else if (!xinfo$minincluded && xinfo$maxincluded) {
          censormax <- domainmax
          censormin <- -Inf
          domainmax <- (8 * domainmax - domainmin) / 7
        } else if (xinfo$minincluded && xinfo$maxincluded) {
          censormin <- domainmin
          censormax <- domainmax
          domainmin <- (7 * domainmin - domainmax) / 6
          domainmax <- (7 * domainmax - domainmin) / 6
        }
        tlocation <- Qf((tlocation - domainmin) / (domainmax - domainmin))
        tscale <- abs(Qf((xinfo$highvalue - domainmin) / (domainmax - domainmin)) -
                       Qf((xinfo$lowvalue - domainmin) / (domainmax - domainmin))) *
          sdoveriqr
      } else if (is.finite(xinfo$domainmin)) {
        transf <- 'log'
        tlocation <- log(tlocation - domainmin)
        tscale <- abs(log(xinfo$highvalue - domainmin) -
                       log(xinfo$lowvalue - domainmin)) * sdoveriqr
      } else if (is.finite(xinfo$domainmax)) {
        transf <- 'logminus'
        tlocation <- log(domainmax - tlocation)
        tscale <- abs(log(domainmax - xinfo$highvalue) -
                     log(domainmax - xinfo$lowvalue)) * sdoveriqr
      }
      if (xinfo$rounding > 0) { # discretized
      ## if(diff(range(x))/xinfo$rounding > 256){
      ##     Message('\nVariate ',xn,' is reported as "rounded".\n
      ##            Consider the possibility of treating it as continuous 
      ##            setting "rounding" to 0.\n')
      ## }
        vtype <- 'D'
        vid <- idD
        idD <- idD + 1L
      } else if (cens) { # censored
        vtype <- 'C'
        vid <- idC
        idC <- idC + 1L
      } else { # continuous
        vtype <- 'R'
        vid <- idR
        idR <- idR + 1L
      }
    } else { # end continuous case
      stop(paste0('ERROR: unknown variate type for ', xn))
    }
    ## Debugging
    ## print(auxmetadata[nrow(auxmetadata)])
    ## print(as.data.table(c(list(name=xn, type=vtype, transform=transf,
    ## Nvalues=vn, step=vd, domainmin=domainmin, domainmax=domainmax,
    ## censormin=censormin, censormax=censormax, tlocation=tlocation,
    ## tscale=tscale, plotmin=plotmin, plotmax=plotmax, Q1=Q1, Q2=Q2, Q3=Q3),
    ## vval
    ## )))

    auxmetadata <- rbind(auxmetadata,
      c(
        list(
          name = xn, mcmctype = vtype, id = vid, # censored=cens,
          ## rounded = rounded, # not used in other scripts, possibly remove
          transform = transf, Nvalues = vn, step = vd,
          domainmin = domainmin, domainmax = domainmax,
          censormin = censormin, censormax = censormax,
          tlocation = tlocation, tscale = tscale,
          plotmin = plotmin, plotmax = plotmax,
          ## Q1 = Q1, Q2 = Q2, Q3 = Q3, # not used in other scripts, possibly remove
          mctest1 = mctest1, mctest2 = mctest2, mctest3 = mctest3
        ),
        vval
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
