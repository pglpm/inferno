# DESCRIPTION
# @param data is either a file or a data.table object.
# @param metadata is either a file or a data.table object.
# @param file Boolean, whether to save the aux metadatafile or not.

buildauxmetadata <- function(data, metadata) {
  ## Elements used in other scripts:
  ## name, mcmctype, Nvalues, mctest, transform,
  ## domainmin, domainmax, tscale, tlocation
  ## censormax, censormin, step,
  ## mctest1, mctest2, mctest3,
  ## plotmin, plotmax, V...

  sdoveriqr <- 0.5 / qnorm(0.75)

  Qf <- readRDS('Qfunction3600_3.rds')

  datafile <- NULL
  if (is.character(data) && file.exists(data)) {
    datafile <- data
    data <- data.table::fread(datafile, na.strings = '')
  data <- data.table::as.data.table(data)
  }
  ##
  ## Q <- readRDS('Qfunction512.rds')
  ##
  idR <- idC <- idD <- idO <- idB <- idN <- 1L
  auxmetadata <- data.table()
  for (xn in metadata$name) {
    if(!is.null(data)) {
      x <- data[[xn]]
      x <- x[!is.na(x)]
    }
    xinfo <- as.list(metadata[name == xn])
    xinfo$type <- tolower(xinfo$type)
    ordinal <- NA
    cens <- any(c(xinfo$minincluded, xinfo$maxincluded), na.rm = TRUE)
    rounded <- NA
    transf <- 'identity' # temporary
    vval <- xinfo[grep('^V[0-9]+$', names(xinfo))]
    ## print(xn)
    ## str(vval)
    Q1 <- NA
    Q2 <- NA
    Q3 <- NA
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
      location <- 0
      scale <- 1
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
    } else if (xinfo$type == 'nominal') { # nominal variate
      vtype <- 'N'
      vid <- idN
      idN <- idN + 1L
      vn <- xinfo$Nvalues
      vd <- 0.5
      domainmin <- 1 # Nimble index categorical from 1
      domainmax <- vn
      censormin <- -Inf
      censormax <- +Inf
      location <- 0
      scale <- 1
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
            mctest2 <- round(xinfo$Nvalues/2)
            mctest3 <- xinfo$Nvalues
    } else if (xinfo$type == 'ordinal') { # ordinal variate
      vtype <- 'O'
      vid <- idO
      idO <- idO + 1L
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
      location <- Qf(round((xinfo$centralvalue - olocation) / oscale) / vn)
      scale <- abs(Qf(round((xinfo$highvalue - olocation) / oscale) / vn) -
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
      location <- xinfo$centralvalue
      scale <- abs(xinfo$highvalue - xinfo$lowvalue)
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
        location <- Qf((location - domainmin) / (domainmax - domainmin))
        scale <- abs(Qf((xinfo$highvalue - domainmin) / (domainmax - domainmin)) -
                       Qf((xinfo$lowvalue - domainmin) / (domainmax - domainmin))) *
          sdoveriqr
      } else if (is.finite(xinfo$domainmin)) {
        transf <- 'log'
        location <- log(location - domainmin)
        scale <- abs(log(xinfo$highvalue - domainmin) -
                       log(xinfo$lowvalue - domainmin)) * sdoveriqr
      } else if (is.finite(xinfo$domainmax)) {
        transf <- 'logminus'
        location <- log(domainmax - location)
        scale <- abs(log(domainmax - xinfo$highvalue) -
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
    ## censormin=censormin, censormax=censormax, tlocation=location,
    ## tscale=scale, plotmin=plotmin, plotmax=plotmax, Q1=Q1, Q2=Q2, Q3=Q3),
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
          tlocation = location, tscale = scale,
          plotmin = plotmin, plotmax = plotmax,
          ## Q1 = Q1, Q2 = Q2, Q3 = Q3, # not used in other scripts, possibly remove
          mctest1 = mctest1, mctest2 = mctest2, mctest3 = mctest3
        ),
        vval
      ),
      fill = FALSE
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
