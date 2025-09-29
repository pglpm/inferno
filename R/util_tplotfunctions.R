#' Plot numeric or character values
#'
#' @description
#' Plot function that modifies and expands the **graphics** package's [graphics::matplot()] function in several ways.
#'
#' @details
#' Some of the additional features provided by `flexiplot` are the following. First, either or both `x` and `y` arguments can be of class [`base::character`]. In this case, axes labels corresponding to the unique values are used (see arguments `xdomain` and `ydomain`). A jitter can also be added to the generated points, via the `xjitter` and `yjitter` switches. Second, it allows for the specification of only a lower or upper limit in the `xlim` and `ylim` arguments. Third, it uses a cleaner plotting style and a default argument `type = 'l'` (line plot) rather than `type = 'p'` (point plot).
#'
#' @param x Numeric or character: vector of x-coordinates. If missing, a numeric vector `1:...` is created having as many values as the rows of `y`.
#' @param y Numeric or character: vector of y coordinates. If missing, a numeric vector `1:...` is created having as many values as the rows of `x`.
#' @param xdomain,ydomain Character or numeric or `NULL` (default): vector of possible values of the variables represented in the `x`- and `y`-axes, in case the `x` or `y` argument is a character vector. The ordering of the values is respected. If `NULL`, then `unique(x)` or `unique(y)` is used.
#' @param xlim,ylim `NULL` (default) or a vector of two values. In the latter case, if any of the two values is not finite (including `NA` or `NULL`), then the `min` or `max` `x`- or `y`-coordinates of the plotted points are used.
#' @param grid Logical: whether to plot a light grid. Default `TRUE`.
#' @param alpha.f Numeric, default 1: opacity of the colours, `0` being completely invisible and `1` completely opaque.
#' @param xjitter,yjitter Logical or `NULL` (default): add [base::jitter()] to `x`- or `y`-values? Useful when plotting discrete variates. If `NULL`, jitter is added if the values are of character class.
#' @param ... Other parameters to be passed to [graphics::matplot()].
#'
#' @export
flexiplot <- function(
    x, y,
    type = NULL,
    lty = c(1, 2, 4, 3, 6, 5),
    lwd = 2,
    pch = c(1, 2, 0, 5, 6, 3), #, 4,
    col = palette(),
    xlab = NULL, ylab = NULL,
    xlim = NULL, ylim = NULL,
    add = FALSE,
    xdomain = NULL, ydomain = NULL,
    alpha.f = 1,
    xjitter = NULL,
    yjitter = NULL,
    ## c( ## Tol's colour-blind-safe scheme
    ##     '#4477AA',
    ##     '#EE6677',
    ##     '#228833',
    ##     '#CCBB44',
    ##     '#66CCEE',
    ##     '#AA3377' #, '#BBBBBB'
    ## ),
    grid = TRUE,
    cex.main = 1,
    ...
){
    xat <- yat <- xaxp <- yaxp <- NULL

    if(missing('x') && !missing('y')){
        x <- y
        x[] <- rep(seq_len(NCOL(y)), each = NROW(y))
        if(is.null(ylab)){ ylab <- deparse1(substitute(y)) }
        if(is.null(yjitter)){ yjitter <- FALSE }
        if(is.null(xdomain) && is.null(xlim)){
            xat <- seq_len(NCOL(y))
            xdomain <- NA
            ## if(!is.null(xjitter)){
            ##     xlim <- range(x) + c(-0.04, 0.04)
            ## }
            if(is.null(xlab)){ xlab <- NA }
            if(is.null(type)){ type <- 'p' }
        }
    } else if(!missing('x') && missing('y')){
        y <- x
        y[] <- rep(seq_len(NCOL(x)), each = NROW(x))
        if(is.null(xlab)){ xlab <- deparse1(substitute(x)) }
        if(is.null(xjitter)){ xjitter <- FALSE }
        if(is.null(ydomain) && is.null(ylim)){
            yat <- seq_len(NCOL(x))
            ydomain <- NA
            ## if(!is.null(yjitter)){
            ##     ylim <- range(y) + c(-0.04, 0.04)
            ## }
            if(is.null(ylab)){ ylab <- NA }
            if(is.null(type)){ type <- 'p' }
        }
    } else if(!missing('x') && !missing('y')){
        if(is.null(xlab)){ xlab <- deparse1(substitute(x)) }
        if(is.null(ylab)){ ylab <- deparse1(substitute(y)) }
    } else {
        stop('Arguments "x" and "y" cannot both be missing')
    }

    if(NROW(y) == 1 && NCOL(y) == NCOL(x)){
        y <- rep(y, each = NROW(x))
        dim(y) <- dim(x)
        if(is.null(type)){ type <- 'p' }
    }
    if(NROW(x) == 1 && NCOL(x) == NCOL(y)){
        x <- rep(x, each = NROW(y))
        dim(x) <- dim(y)
        if(is.null(type)){ type <- 'p' }
    }

    if(is.character(x) && is.character(y)) {
        if(is.null(xjitter)){xjitter <- TRUE}
        if(is.null(yjitter)){yjitter <- TRUE}
    }
    ## if x is character, convert to numeric
    if(is.character(x)){
        if(is.null(xdomain)){ xdomain <- unique(x) }
        ## we assume the user has sorted the vaules in a meaningful order
        ## because the lexical order may not be correct
        ## (think of values like 'low', 'medium', 'high')
        . <- dim(x)
        x <- as.numeric(factor(x, levels = xdomain))
        dim(x) <- .
        xat <- seq_along(xdomain)
        xaxp <- c(range(xat), length(xat) - 1)
        if(is.null(type)){ type <- 'p' }
    }
    if(isTRUE(xjitter)){
        xaxp <- c(range(xat) + c(-0.5, 0.5), length(xat))
        ## xaxp <- c(range(xat), length(xat) - 1)
        x <- jitter(x, factor = 5/3)
    }

    ## if y is character, convert to numeric
    if(is.character(y)){
        if(is.null(ydomain)){ ydomain <- unique(y) }
        ## we assume the user has sorted the vaules in a meaningful order
        ## because the lexical order may not be correct
        ## (think of values like 'low', 'medium', 'high')
        . <- dim(y)
        y <- as.numeric(factor(y, levels = ydomain))
        dim(y) <- .
        yat <- seq_along(ydomain)
        yaxp <- c(range(yat), length(yat) - 1)
        if(is.null(type)){ type <- 'p' }
    }
    if(isTRUE(yjitter)){
        yaxp <- c(range(yat) + c(-0.5, 0.5), length(yat))
        ## yaxp <- c(range(yat), length(yat) - 1)
        y <- jitter(y, factor = 5/3)
    }

    ## Syntax of xlim and ylim that allows
    ## for the specification of only upper- or lower-bound
    if(length(xlim) == 2){
        if(is.null(xlim[1]) || !is.finite(xlim[1])){ xlim[1] <- min(x[is.finite(x)]) }
        if(is.null(xlim[2]) || !is.finite(xlim[2])){ xlim[2] <- max(x[is.finite(x)]) }
    }
    if(length(ylim) == 2){
        if(is.null(ylim[1]) || !is.finite(ylim[1])){ ylim[1] <- min(y[is.finite(y)]) }
        if(is.null(ylim[2]) || !is.finite(ylim[2])){ ylim[2] <- max(y[is.finite(y)]) }
    }

    if(is.null(type)){ type <- 'l' }

    if(is.na(alpha.f)){alpha.f <- 1}
    col <- adjustcolor(col, alpha.f = alpha.f)
    graphics::matplot(x, y, xlim = xlim, ylim = ylim, type = type, axes = FALSE,
        col = col, lty = lty, lwd = lwd, pch = pch, cex.main = cex.main, add = add, xlab = xlab, ylab = ylab, ...)
    if(!add){
        graphics::axis(1, at = xat, labels = xdomain, tick = !grid,
            col = 'black', lwd = 1, lty = 1, ...)
        graphics::axis(2, at = yat, labels = ydomain, tick = !grid,
            col = 'black', lwd = 1, lty = 1, ...)
        if(grid){
            if(exists('xaxp')){ par(xaxp = xaxp) }
            if(exists('yaxp')){ par(yaxp = yaxp) }
            graphics::grid(nx = NULL, ny = NULL, lty = 1, col = '#BBBBBB80')
        }
    }
}
## flexiplot <- function(
##     x, y,
##     type = 'l',
##     lty = c(1, 2, 4, 3, 6, 5),
##     lwd = 2,
##     pch = c(1, 2, 0, 5, 6, 3), #, 4,
##     col = palette(),
##     xlab = NULL, ylab = NULL,
##     xlim = NULL, ylim = NULL,
##     add = FALSE,
##     xdomain = NULL, ydomain = NULL,
##     alpha.f = 1,
##     xjitter = NULL,
##     yjitter = NULL,
##     ## c( ## Tol's colour-blind-safe scheme
##     ##     '#4477AA',
##     ##     '#EE6677',
##     ##     '#228833',
##     ##     '#CCBB44',
##     ##     '#66CCEE',
##     ##     '#AA3377' #, '#BBBBBB'
##     ## ),
##     grid = TRUE,
##     cex.main = 1,
##     ...
## ){
##     xat <- yat <- NULL
## 
##     if(missing('x') && !missing('y')){
##         x <- numeric(NROW(y))
##         if(is.null(xdomain) && is.null(xlim)){
##             xat <- 0
##             xdomain <- NA
##             if(!is.null(xjitter)){
##                 xlim <- c(-0.04, 0.04)
##             }
##             if(is.null(xlab)){ xlab <- NA }
##             if(is.null(ylab)){ ylab <- deparse1(substitute(y)) }
##         }
##     } else if(!missing('x') && missing('y')){
##         y <- numeric(NROW(x))
##         if(is.null(ydomain) && is.null(ylim)){
##             yat <- 0
##             ydomain <- NA
##             if(!is.null(yjitter)){
##                 ylim <- c(-0.04, 0.04)
##             }
##             if(is.null(ylab)){ ylab <- NA }
##             if(is.null(xlab)){ xlab <- deparse1(substitute(x)) }
##         }
##     } else if(!missing('x') && !missing('y')){
##         if(is.null(xlab)){ xlab <- deparse1(substitute(x)) }
##         if(is.null(ylab)){ ylab <- deparse1(substitute(y)) }
##     } else {
##         stop('Arguments "x" and "y" cannot both be missing')
##     }
## 
##     ## if x is character, convert to numeric
##     if(is.character(x)){
##         if(is.null(xdomain)){ xdomain <- unique(x) }
##         ## we assume the user has sorted the vaules in a meaningful order
##         ## because the lexical order may not be correct
##         ## (think of values like 'low', 'medium', 'high')
##         x <- as.numeric(factor(x, levels = xdomain))
##         if(is.null(xjitter)){xjitter <- TRUE}
##         xat <- seq_along(xdomain)
##     }
##     if(isTRUE(xjitter)){x <- jitter(x)}
## 
##     ## if y is character, convert to numeric
##     if(is.character(y)){
##         if(is.null(ydomain)){ ydomain <- unique(y) }
##         ## we assume the user has sorted the vaules in a meaningful order
##         ## because the lexical order may not be correct
##         ## (think of values like 'low', 'medium', 'high')
##         y <- as.numeric(factor(y, levels = ydomain))
##         if(is.null(yjitter)){yjitter <- TRUE}
##         yat <- seq_along(ydomain)
##     }
##     if(isTRUE(yjitter)){y <- jitter(y)}
## 
##     ## Syntax of xlim and ylim that allows
##     ## for the specification of only upper- or lower-bound
##     if(length(xlim) == 2){
##         if(is.null(xlim[1]) || !is.finite(xlim[1])){ xlim[1] <- min(x[is.finite(x)]) }
##         if(is.null(xlim[2]) || !is.finite(xlim[2])){ xlim[2] <- max(x[is.finite(x)]) }
##     }
##     if(length(ylim) == 2){
##         if(is.null(ylim[1]) || !is.finite(ylim[1])){ ylim[1] <- min(y[is.finite(y)]) }
##         if(is.null(ylim[2]) || !is.finite(ylim[2])){ ylim[2] <- max(y[is.finite(y)]) }
##     }
## 
##     if(is.null(xlab) && !missing(x)) {
##         xlab <- deparse1(substitute(x))
##     }
##     if(is.na(alpha.f)){alpha.f <- 1}
##     col <- adjustcolor(col, alpha.f = alpha.f)
## 
##     graphics::matplot(x, y, xlim = xlim, ylim = ylim, type = type, axes = F,
##         col = col, lty = lty, lwd = lwd, pch = pch, cex.main = cex.main, add = add, xlab = xlab, ylab = ylab, ...)
##     if(!add){
##         graphics::axis(1, at = xat, labels = xdomain, tick = !grid,
##             col = 'black', lwd = 1, lty = 1, ...)
##         graphics::axis(2, at = yat, labels = ydomain, tick = !grid,
##             col = 'black', lwd = 1, lty = 1, ...)
##         if(grid){
##             graphics::grid(nx = NULL, ny = NULL, lty = 1, col = '#BBBBBB80')
##         }
##     }
## }

#' Plot pairs of quantiles
#'
#' @description
#' Utility function to plot pair of quantiles obtained with [Pr()].
#'
#' @param x Numeric or character: vector of x-coordinates. See [flexiplot()].
#' @param y Numeric: a matrix having as many rows as `x` and an even number of columns, with one column per quantile. Typically these quantiles have been obtained with [Pr()], as their `$quantiles` value. This value is a three-dimensional array, and one of its columns (corresponding to the possible values of the `X` argument of [Pr()]) or one of its rows (corresponding to the possible values of the `Y` argument of [Pr()]) should be selected before being used as `y` input.
#' @param xdomain Character or numeric or `NULL` (default): vector of possible values of the variable represented in the x-axis, if the `x` argument is a character vector. The ordering of the values is respected. If `NULL`, then `unique(x)` is used.
#' @param alpha.f Numeric, default 0.25: opacity of the quantile bands, `0` being completely invisible and `1` completely opaque.
#' @param col Fill colour of the quantile bands. Can be specified in any of the usual ways, see for instance [grDevices::col2rgb()]. Default `#4477AA`.
#' @param border Fill colour of the quantile bands. Can be specified in any of the usual ways, see for instance [grDevices::col2rgb()]. If `NA` (default), no border is drawn.
#' @param ... Other parameters to be passed to [flexiplot()].
#'
#' @export
plotquantiles <- function(
    x, y,
    xdomain = NULL,
    alpha.f = 0.25,
    col = palette(),
    ##     c( ## Tol's colour-blind-safe scheme
    ##     '#4477AA',
    ##     '#EE6677',
    ##     '#228833',
    ##     '#CCBB44',
    ##     '#66CCEE',
    ##     '#AA3377' #, '#BBBBBB'
    ## ),
    border = NA,
    type = 'n',
    ...
){
    ## ## TODO: modify so that a vertical plot is also possible
    if(!is.matrix(y) || ncol(y) %% 2 != 0) {
        stop('"y" must be a matrix with an even number of columns.')
    }
    nquant <- ncol(y)

    isfin <- ( (is.numeric(x) & is.finite(x)) | !is.na(x)) &
        apply(y, 1, function(xx){all(is.finite(xx))})
    x <- unname(x[isfin])
    y <- unname(y[isfin, , drop = FALSE])

    ##
    ## col[!grepl('^#',col)] <- palette()[as.numeric(col[!grepl('^#',col)])]
    ## if(is.na(alpha.f)){alpha.f <- 1}
    ## col <- adjustcolor(col, alpha.f = alpha.f)
    ## if(is.na(alpha)){alpha <- ''}
    ## else if(!is.character(alpha)){alpha <- alpha2hex(alpha)}
    ## if(!(is.na(col) | nchar(col)>7)){col <- paste0(col, alpha)}
    ##
    flexiplot(x = x, y = y, xdomain = xdomain, type = 'n',
        xjitter = FALSE, yjitter = FALSE, ...)

    ## if x is character, convert to numeric
    if(is.character(x)){
        if(is.null(xdomain)){ xdomain <- unique(x) }
        ## we assume the user has sorted the vaules in a meaningful order
        ## because the lexical order may not be correct
        ## (think of values like 'low', 'medium', 'high')
        x <- as.numeric(factor(x, levels = xdomain))
    }
    if(is.na(alpha.f)){alpha.f <- 1}
    col <- adjustcolor(col, alpha.f = alpha.f)
    for(ii in seq_len(nquant/2)) {
        graphics::polygon(x=c(x, rev(x)), y=c(y[,ii], rev(y[, nquant + 1 - ii])),
            col = col, border = border)
    }
}

#' Plot an object of class "probability"
#'
#' @description
#' This [base::plot()] method is a utility to plot probabilities obtained with [Pr()], as well as their variabilities. The probabilities are plotted either against `Y`, with one curve for each value of `X`, or vice versa.
#'
#' @param p Object of class "probability", obtained with [Pr()].
#' @param variability One of the values `'quantiles'`, `'samples'`, `'none'` (equivalent to `NA` or `FALSE`), or `NULL` (default), in which case the variability available in `p` is used. This argument chooses how to represent the variability of the probability; see [Pr()]. If the requested variability is not available in the object `p`, then a warning is issued and no variability is plotted.
#' @param PvsY Logical or `NULL`: should probabilities be plotted against their `Y` argument? If `NULL`, the argument between `Y` and `X` having larger number of values is chosen. As many probability curves will be plotted as the number of values of the other argument.
#' @param legend One of the values `'bottomright'`, `'bottom'`, `'bottomleft'`, `'left'`, `'topleft'`, `'top'`, `'topright'`, `'right'`, `'center'` (see [graphics::legend()]): plot a legend at that position. A value `FALSE` or any other does not plot any legend. Default `'top'`.
#' @param alpha.f Numeric, default 0.25: opacity of the colours, `0` being completely invisible and `1` completely opaque.
#' @param var.alpha.f Numeric: opacity of the quantile bands or of the samples, `0` being completely invisible and `1` completely opaque.
#' @param ... Other parameters to be passed to [flexiplot()].
#'
#' @export
plot.probability <- function(
    p,
    variability = NULL,
    PvsY = NULL,
    legend = 'top',
    lty = c(1, 2, 4, 3, 6, 5),
    lwd = 2,
    col = palette(),
    type = NULL,
    ##     c( ## Tol's colour-blind-safe scheme, or palette()
    ##     '#4477AA',
    ##     '#EE6677',
    ##     '#228833',
    ##     '#CCBB44',
    ##     '#66CCEE',
    ##     '#AA3377' #, '#BBBBBB'
    ## ),
    alpha.f = 1,
    var.alpha.f = NULL,
    xlab = NULL,
    ylab = NULL,
    main = NULL,
    ylim = c(0, NA),
    grid = TRUE,
    add = FALSE,
    ...
){
    ## Check how we should represent the variability
    ## The user can choose among three options
    ## provided that option is available in argument 'p'

    if(is.null(variability)) { # User is not choosing
        ## We choose 'quantiles' or what's available
        if(!is.null(p$quantiles)) {
            variability <- 'quantiles'
        } else if(!is.null(p$samples)){
            variability <- 'samples'
        } else {
            variability <- 'none'
        }
    } else { # User is choosing
        if(is.na(variability) || isFALSE(variability)){ variability <- 'none'}

        ## handle shortenings
        variability <- match.arg(variability, c('quantiles', 'samples', 'none'))

        ## handle impossible requests
        if(
        (variability == 'quantiles' && is.null(p$quantiles)) ||
            (variability == 'samples' && is.null(p$samples))
        ) {
            message('Requested variability not available. Omitting its plot.')
            variability <- 'none'
        }
    }

    Ylen <- nrow(p$values)
    Xlen <- ncol(p$values)

    ## We rename the variability object so as to avoid if-else below
    if(variability == 'quantiles'){
        mainpercentiles <- c(5.5, 94.5) # By default we choose an 89% band
        pvar <- p$quantiles
        ## if we are only plotting more than one curve, just keep the 89% band
        if(Xlen > 1 && Ylen > 1){
            qnames <- as.numeric(sub('%', '', dimnames(pvar)[[3]]))
            choosepercentiles <- sapply(mainpercentiles,
                function(xx){which.min(abs(qnames - xx))})
            pvar <- pvar[, , choosepercentiles, drop = FALSE]
            qnames <- as.numeric(sub('%', '', dimnames(pvar)[[3]]))
        }
        qnames <- as.numeric(sub('%', '', dimnames(pvar)[[3]]))
        if(is.null(var.alpha.f)){var.alpha.f <- 0.25}
    } else if(variability == 'samples'){
        pvar <- p$samples
        if(is.null(var.alpha.f)){var.alpha.f <- 1/ceiling(sqrt(dim(pvar)[3]))}
    } else {
        pvar <- NULL
    }

    ## Handle the case of missing Y and X items in 'p'
    if(is.null(p$Y)){
        p$Y <- data.frame(Y = paste0('Y', seq_len(Ylen)))
        if(Xlen > 1){
            p$X <- data.frame(X = paste0('X', seq_len(Xlen)))
        }
    }

    ## If there's only one probability it doesn't make sense to plot anything
    if(length(p$values) == 1){
        print(sort(c(Pr = p$values, pvar[1, 1, ])))
    }

    ## If 'PvsY' is NULL, then we guess that the longest between Y and X
    ## is meant to be abscissa
    if(is.null(PvsY)){ PvsY <- (Ylen >= Xlen) }

    if(isTRUE(PvsY)){
        x <- p$Y
        leg <- p$X
        tempxlab <- 'Y'
    } else {
        x <- p$X
        leg <- p$Y
        tempxlab <- 'X'
        p$values <- t(p$values)
        if(!is.null(pvar)){ pvar <- aperm(pvar, c(2, 1, 3)) }
    }

    ## If the abscissa has more than one variate,
    ## then it becomes tricky to understand which of these we must plot against
    ## Heuristic: if there's one variate with as many unique elements as x,
    ## then use that one. Otherwise use a generic 'Y...'
    if(ncol(x) == 1){
        tempxlab <- colnames(x)
        x <- unlist(x)
    } else {
        uniquevrts <- apply(x, 2, function(xx){length(unique(xx))})
        toselect <- which(uniquevrts == nrow(x))[1]
        if(is.na(toselect)){
            x <- seq_len(nrow(x))
        } else {
            tempxlab <- colnames(x)[toselect]
            x <- x[, toselect]
        }
    }

    if(is.null(xlab)){xlab <- tempxlab}
    if(missing(main)){
        main <- paste0('prob. & freq.  ',
            paste0(names(p$Y), collapse = ', '),
            if(!is.null(p$X)){paste0(' | ', paste0(names(p$X), collapse = ', '))}
            )
        if(variability == 'quantiles'){
            main <- paste0(main, '  [',
                paste0(round(qnames, 1), '%', collapse = ', '),
                ']')
        }
    }
    if(is.null(ylab)){
        ylab <- 'probability & rel.frequency'
    }

    if(is.null(type)){
        if(is.character(x)){type <- 'b'} else {type <- 'l'}
    }

    ## Plot the variability first
    ## find maximum and minimum y-value first, if needed
    if(is.na(ylim[2])){
        ylim[2] <- max(pvar, p$values)
    }
    if(is.na(ylim[1])){
        ylim[1] <- min(pvar, p$values)
    }
    if(variability == 'quantiles' && length(p$values) > 1){
        for(i in seq_len(dim(pvar)[2])){
            plotquantiles(x = unlist(x), y = pvar[, i, ],
                col = col[(i - 1) %% length(col) + 1],
                alpha.f = var.alpha.f,
                lty =  lty[(i - 1) %% length(lty) + 1],
                xlab = xlab,
                ylab = ylab,
                ylim = ylim,
                main = main,
                grid = grid,
                add = (add || i > 1),
                ...)
            add <- TRUE
        }

    } else if(variability == 'samples' && length(p$values) > 1){
        ## the samples are plotted alternating between the different subgroups,
        ## rather than one group at a time, in order to avoid that
        ## the samples of the last subgroup cover the previous ones
        nx <- dim(pvar)[2]
        dim(pvar) <- c(dim(pvar)[1], prod(dim(pvar)[-1]))
        flexiplot(x = x, y = pvar,
            type = type,
            col = col[(seq_len(nx) - 1) %% length(col) + 1],
            alpha.f = var.alpha.f,
            lty =  1, #lty[(seq_len(nx) - 1) %% length(lty) + 1],
            lwd = 0.5, #lwd[(seq_len(nx) - 1) %% length(lwd) + 1] / 4,
            xlab = xlab,
            ylab = ylab,
            ylim = ylim,
            main = main,
            grid = grid,
            xjitter = FALSE,
            yjitter = FALSE,
            add = add,
            ...)
        add <- TRUE
    }

    ## Plot the probabilities
    if(length(p$values) > 1){
        flexiplot(x = x, y = p$values,
            type = type,
            col = col,
            alpha.f = alpha.f,
            lty = lty,
            lwd = lwd,
            xlab = xlab,
            ylab = ylab,
            ylim = ylim,
            main = main,
            grid = grid,
            xjitter = FALSE,
            yjitter = FALSE,
            add = add,
            ...)

        ## Plot legends
    if(!is.null(leg) && is.character(legend) &&
           (legend %in%
                c("bottomright", "bottom", "bottomleft", "left", "topleft",
                    "top", "topright", "right", "center"))){
        graphics::legend(x = legend,
            legend = apply(leg, 1, function(xxx){
                paste0(paste0(names(xxx), ' = ', xxx), collapse = ', ')
            }),
            bty = 'n',
            col = col,
            lty = lty,
            lwd = lwd,
            ...)
    }
    }
}


#' Plot the variability of an object of class "probability" as a histogram
#'
#' @description
#' This [graphics::hist()]ogram method is a utility to visualize the variability of the probabilities obtained with [Pr()], which can also be interpreted as the probability density for the whole-population frequencies.
#'
#' @param p Object of class "probability", obtained with [Pr()].
#' @param breaks `NULL` or as in function [graphics::hist()]. If `NULL` (default), an optimal number of breaks for each probability distribution is computed.
#' @param fill.alpha.f Numeric, default 0.125: opacity of the histogram filling. `0` means no filling.
#' @param legend One of the values `"bottomright"`, `"bottom"`, `"bottomleft"`, `"left"`, `"topleft"`, `"top"`, `"topright"`, `"right"`, `"center"` (see [graphics::legend()]): plot a legend at that position. A value `FALSE` or any other does not plot any legend. Default `"top"`.
#' @param showmean Logical, default `TRUE`: show the means of the probability distributions? The means correspond to the probabilities about the next observed unit.
#' @param ... Other parameters to be passed to [flexiplot()].
#'
#' @export
hist.probability <- function(
    p,
    breaks = NULL,
    legend = 'top',
    lty = c(1, 2, 4, 3, 6, 5),
    lwd = 2,
    col = palette(),
    alpha.f = 1,
    fill.alpha.f = 0.125,
    showmean = TRUE,
    ##     c( ## Tol's colour-blind-safe scheme, or palette()
    ##     '#4477AA',
    ##     '#EE6677',
    ##     '#228833',
    ##     '#CCBB44',
    ##     '#66CCEE',
    ##     '#AA3377' #, '#BBBBBB'
    ## ),
    xlab = NULL,
    ylab = NULL,
    xlim = NULL,
    ylim = c(0, NA),
    main = NULL,
    grid = TRUE,
    add = FALSE,
    ...
){
    ## Check that samples are available in the probability object

    if(is.null(p$samples)) {
        stop('The probability object does not contain any frequency samples')
        }
    pvar <- p$samples
    Ylen <- nrow(p$values)
    Xlen <- ncol(p$values)

    if(is.null(breaks)){n <- ceiling(sqrt(dim(pvar)[3])/2)} else {n <- NULL}

    ## Precompute histograms, to determine maximum y-value
    midslist <- densitylist <- list()
    i <- 0L
    for(xx in seq_len(Xlen)){ for(yy in seq_len(Ylen)){
        i <- i + 1L
        ff <- pvar[yy, xx, ]
        rg <- range(ff)
        if(diff(rg)==0){rg <- c(0, 1)}
        if(!is.null(n)){ breaks <- seq(rg[1], rg[2], length.out = n + 1) }
        hd <- graphics::hist(x = ff, breaks = breaks, plot = FALSE)
        midslist[[i]] <- hd$mids
        densitylist[[i]] <- hd$density
    } }

    if(is.null(xlab)){xlab <- 'relative frequency'}
    if(is.null(ylab)){ylab <- 'probability density'}
    if(isFALSE(fill.alpha.f) || !is.numeric(fill.alpha.f)){fill.alpha.f <- 0}

    if(missing(xlim)){xlim <- range(unlist(midslist))}
    if(is.na(ylim)[2]){ylim[2] <- max(unlist(densitylist))}

    i <- 0L
    for(xx in seq_len(Xlen)){ for(yy in seq_len(Ylen)){
        i <- i + 1L
        x <- midslist[[i]]
        y <- densitylist[[i]]
        thiscol <- col[(i - 1) %% length(col) + 1]
        thislty <- lty[(i - 1) %% length(lty) + 1]
        if(alpha.f > 0){
            plotquantiles(x = x, y = cbind(rep(0, length(y)), y),
                col = thiscol,
                alpha.f = fill.alpha.f,
                xlab = xlab, ylab = ylab,
                xlim = xlim, ylim = ylim,
                main = main,
                grid = grid,
                lty = 0,
                add = (add || i > 1),
                ...)
        }

        flexiplot(x = midslist[[i]], y = densitylist[[i]],
            xlab = xlab, ylab = ylab,
            xlim = xlim, ylim = ylim,
            main = main,
            col = thiscol,
            alpha.f = alpha.f,
            lty = thislty,
            lwd = lwd,
            grid = grid,
            xjitter = FALSE,
            yjitter = FALSE,
            add = (add || alpha.f >0 || i > 1),
            ...
        )
        if(isTRUE(showmean)){
            graphics::abline(v = p$values[yy, xx],
                col = adjustcolor(thiscol, alpha.f * 0.75),
                lty = thislty,
                lwd = lwd * 0.75)
        }
    } }

    ## Plot legends
    if(is.character(legend) &&
           (legend %in%
                 c("bottomright", "bottom", "bottomleft", "left", "topleft",
                     "top", "topright", "right", "center"))){
        legs <- paste0(apply(p$Y, 1, function(xxx){
                paste0(paste0(names(xxx), ' = ', xxx), collapse = ', ')
        }))
        if(!is.null(p$X)){
            legs <- c(outer( legs,
                paste0(' | ',
                    apply(p$X, 1, function(xxx){
                        paste0(paste0(names(xxx), ' = ', xxx), collapse = ', ')
                    })
                ),
                paste0))
        }

        graphics::legend(x = legend,
            legend = legs,
            bty = 'n',
            col = col,
            lty = lty,
            lwd = lwd,
            ...)
    }

}
