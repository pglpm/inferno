#' Plot numeric or character values
#'
#' Plot function that modifies and expands the `graphics` \code{\link[graphics]{matplot}} in several ways: First, either or both `x` and `y` arguments can be of class \code{\link[base]{character}}. In this case, axes labels corresponding to the unique values are used (see arguments `xdomain` and 'ydomain'). Second, it allows for the specification of only a lower or upper limit in `xlim` and `ylim`. Third, it uses a cleaner plotting style and a default argument `type = 'l'` (line plot) rather than `type = 'p'` (point plot).
#'
#' @param x Numeric or character: vector of x-coordinates. If missing, a numeric vector `1:...` is created having as many values as the rows of `y`.
#' @param y Numeric or character: vector of y coordinates. If missing, a numeric vector `1:...` is created having as many values as the rows of `x`.
#' @param xdomain Character or numeric or `NULL` (default): vector of possible values of the variable represented in the x-axis, if the `x` argument is a character vector. The ordering of the values is respected. If `NULL`, then `unique(x)` is used.
#' @param ydomain Character or numeric or `NULL` (default): like `xdomain` but for the y-coordinate.
#' @param xlim `NULL` (default) or a vector of two values. In the latter case, if any of the two values is not finite (including `NA` or `NULL`), then the `min` or `max` x-coordinate of the plotted points is used.
#' @param ylim `NULL` (default) or a vector of two values. Like argument `xlim`, but for the y-coordinates.
#' @param grid Logical: whether to plot a light grid. Default `TRUE`.
#' @param ... Other parameters to be passed to \code{\link[base]{matplot}}.
#'
#' @export
flexiplot <- function(
    x, y,
    xdomain = NULL, ydomain = NULL,
    xlim = NULL, ylim = NULL,
    type = 'l',
    pch = c(1, 0, 2, 5, 6, 3, 4),
    grid = TRUE,
    add = FALSE,
    lwd = 1,
    ...
){
    xat <- yat <- NULL

    if(missing('x') && !missing('y')){
        x <- seq_len(NROW(y))
    } else if(!missing('x') && missing('y')){
        y <- seq_len(NROW(x))
    } else if(missing('x') && missing('y')){
        stop('Arguments "x" and "y" cannot both be missing')
    }

    ## if x is character, convert to numeric
    if(is.character(x)){
        if(is.null(xdomain)){ xdomain <- unique(x) }
        ## we assume the user has sorted the vaules in a meaningful order
        ## because the lexical order may not be correct
        ## (think of values like 'low', 'medium', 'high')
        x <- as.numeric(factor(x, levels = xdomain))
        xat <- seq_along(xdomain)
    }

    ## if y is character, convert to numeric
    if(is.character(y)){
        if(is.null(ydomain)){ ydomain <- unique(y) }
        ## we assume the user has sorted the vaules in a meaningful order
        ## because the lexical order may not be correct
        ## (think of values like 'low', 'medium', 'high')
        y <- as.numeric(factor(y, levels = ydomain))
        yat <- seq_along(ydomain)
    }

    ## Syntax of xlim and ylim that allows
    ## for the specification of only upper- or lower-bond
    if(length(xlim) == 2){
        if(is.null(xlim[1]) || !is.finite(xlim[1])){ xlim[1] <- min(x[is.finite(x)]) }
        if(is.null(xlim[2]) || !is.finite(xlim[2])){ xlim[2] <- max(x[is.finite(x)]) }
    }
    if(length(ylim) == 2){
        if(is.null(ylim[1]) || !is.finite(ylim[1])){ ylim[1] <- min(y[is.finite(y)]) }
        if(is.null(ylim[2]) || !is.finite(ylim[2])){ ylim[2] <- max(y[is.finite(y)]) }
    }

    graphics::matplot(x, y, xlim = xlim, ylim = ylim, type = type, pch = pch, axes = F, add = add, lwd = lwd, ...)
    if(!add){
        graphics::axis(1, at = xat, labels = xdomain, lwd = 0, ...)
        graphics::axis(2, at = yat, labels = ydomain, lwd = 0, ...)
        if(grid){
            graphics::grid(nx = NULL, ny = NULL, lty = 1, col = '#BBBBBB80')
        }
    }
}

#' Plot pairs of quantiles
#'
#' Utility function to plot pair of quantiles obtained with \code{\link{Pr}} or \code{\link{tailPr}}.
#'
#' @param x Numeric or character: vector of x-coordinates. See \code{\link{flexiplot}}.
#' @param y Numeric: a matrix having as many rows as `x` and an even number of columns, with one column per quantile. Typically these quantiles have been obtained with \code{\link{Pr}} or \code{\link{tailPr}}, as their `$quantiles` value. This value is a three-dimensional array, and one of its columns (corresponding to the possible values of the `X` argument of \code{\link{Pr}}) or one of its rows (corresponding to the possible values of the `Y` argument of \code{\link{Pr}}) should be selected before being used as `y` input.
#' @param xdomain Character or numeric or `NULL` (default): vector of possible values of the variable represented in the x-axis, if the `x` argument is a character vector. The ordering of the values is respected. If `NULL`, then `unique(x)` is used.
#' @param alpha.f Numeric, default 0.25: opacity of the quantile bands, `0` being completely invisible and `1` completely opaque.
#' @param col Fill colour of the quantile bands. Can be specified in any of the usual ways, see for instance \code{\link[grDevices]{col2rgb}}. Default `9`.
#' @param border Fill colour of the quantile bands. Can be specified in any of the usual ways, see for instance \code{\link[grDevices]{col2rgb}}. If `NA` (default), no border is drawn.
#' @param ... Other parameters to be passed to \code{\link{flexiplot}}.
#'
#' @export
plotquantiles <- function(
    x, y,
    xdomain = NULL,
    alpha.f = 0.25,
    col = 9,
    border = NA,
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
    if(is.na(alpha.f)){alpha.f <- 1}
    col <- adjustcolor(col, alpha.f = alpha.f)
    ## if(is.na(alpha)){alpha <- ''}
    ## else if(!is.character(alpha)){alpha <- alpha2hex(alpha)}
    ## if(!(is.na(col) | nchar(col)>7)){col <- paste0(col, alpha)}
    ##
    flexiplot(x = x, y = y, xdomain = xdomain, type = 'n', ...)

    ## if x is character, convert to numeric
    if(is.character(x)){
        if(is.null(xdomain)){ xdomain <- unique(x) }
        ## we assume the user has sorted the vaules in a meaningful order
        ## because the lexical order may not be correct
        ## (think of values like 'low', 'medium', 'high')
        x <- as.numeric(factor(x, levels = xdomain))
    }

    for(ii in seq_len(nquant/2)) {
        graphics::polygon(x=c(x, rev(x)), y=c(y[,ii], rev(y[, nquant + 1 - ii])),
            col = col, border = border)
    }
}

#' Plot an object of class "probability"
#'
#' This plot method is a utility to plot probabilities obtained with \code{\link{Pr}} or \code{\link{tailPr}}, as well as their uncertainties. The probabilities are plotted either against `Y`, with one curve for each value of `X`, or vice versa.
#'
#' @param p Object of class "probability", obtained with \code{\link{Pr}} or \code{\link{tailPr}}.
#' @param variability One of the values `"quantiles"`, `"samples"`, `"none"` (equivalent to `NA` or `FALSE`), or `NULL` (default), in which case the variability available in `p` is used. This argument chooses how to represent the variability of the probability; see \code{\link{Pr}}. If the requested variability is not available in the object `p`, then a warning is issued and no variability is plotted.
#' @param PvsY Logical or `NULL`: should probabilities be plotted against their `Y` argument? If `NULL`, the argument between `Y` and `X` having larger number of values is chosen. As many probability curves will be plotted as the number of values of the other argument.
#' @param legend Logical: plot a legend of the different curves?
#' @param ... Other parameters to be passed to \code{\link[base]{matplot}}.
#'
#' @export
plot.probability <- function(
    p,
    variability = NULL,
    PvsY = NULL,
    legend = TRUE,
    lwd = 3,
    lty = 1:5,
    col = palette(),
    xlab = NULL,
    ylab = NULL,
    ylim = c(0, NA),
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
            warning('Requested variability not available. Omitting its plot.')
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
    } else if(variability == 'samples'){
        pvar <- p$samples
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
        pvar <- aperm(pvar, c(2, 1, 3))
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
    if(is.null(ylab)){
        ylab <- 'probability'
        if(variability == 'quantiles'){
            ylab <- paste0(ylab, ' (variab. ',
                paste0(round(qnames, 1), '%', collapse = ', '),
                ')')
        }
    }

    ## Plot the variability first
    ## find maximum and minimum y-value first, if needed
    if(is.na(ylim[2])){
        ylim[2] <- max(pvar)
    }
    if(is.na(ylim[1])){
        ylim[1] <- min(pvar)
    }
    if(variability == 'quantiles' && length(p$values) > 1){
        for(i in seq_len(dim(pvar)[2])){
            plotquantiles(x = unlist(x), y = pvar[, i, ],
                col = col[(i - 1) %% length(col) + 1],
                lty =  lty[(i - 1) %% length(lty) + 1],
                xlab = xlab,
                ylab = ylab,
                ylim = ylim,
                add = (add || i > 1),
                ...)
            add <- TRUE
        }

    } else if(variability == 'samples' && length(p$values) > 1){
        nx <- dim(pvar)[2]
        dim(pvar) <- c(dim(pvar)[1], prod(dim(pvar)[-1]))
        flexiplot(x = x, y = pvar,
            col = adjustcolor(col[(seq_len(nx) - 1) %% length(col) + 1],
                alpha.f = 0.25),
            lty =  lty[(seq_len(nx) - 1) %% length(lty) + 1],
            lwd = lwd[(seq_len(nx) - 1) %% length(lwd) + 1] / 4,
            xlab = xlab,
            ylab = ylab,
            ylim = ylim,
            add = add,
            ...)
        add <- TRUE
    }

    ## Plot the probabilities
    if(length(p$values) > 1){
        flexiplot(x = x, y = p$values,
            col = col,
            lty = lty,
            lwd = lwd,
            xlab = xlab,
            ylab = ylab,
            ylim = ylim,
            add = add,
            ...)

        ## Plot legends
        if(!is.null(leg) && legend){
            graphics::legend(x = 'topright',
                legend = apply(leg, 1, paste0, collapse = ', '),
                bty = 'n',
                col = col,
                lty = lty,
                lwd = lwd,
                ...)
        }
    }
}


#### Old functions below, deleting them soon


## #' Plot probabilities
## #'
## #' Utility function to plot probabilities obtained with \code{\link{Pr}}. The advantage of this function over \code{\link[base]{plot}} is that its `x` argument can be a vector of character values, as typically the case for nominal or ordinal variates (see \code{\link{metadata}}). The plot will have these characters as labels of the x-axis.
## #'
## #' @param x Numeric or character: Set of values that were used as `Y` argument to \code{\link{Pr}} (typically they make sense if `Y` only was a single variate rather than aconjunction of several).
## #' @param y Numeric: a vector of same length as `x`, or a matrix having as many rows as `x`. Typically this is a vector of probabilities obtained with \code{\link{Pr}}, as its `$values` value, on as its `$samples` value.
## #'
## #' @export
## tplot <- function(x, y, xlim = c(NA, NA), ylim = c(NA, NA), asp = NA,
##     n = 10, family = '', xticks = NULL, xlabels = TRUE,
##     yticks = NULL, ylabels = TRUE, cex = 1.5, ly = NULL,
##     lx = NULL, mar = NULL, lty.axis = 1, lwd.axis = 0,
##     lwd.ticks = 1, col.ticks = '#bbbbbb80', col.lab = 'black',
##     cex.axis = 1.12, las.y = 1, xgrid = NULL, ygrid = NULL,
##     main = NULL, cex.main = 1.5, xlab = NULL, ylab = NULL,
##     cex.lab = 1.5, type = 'l', col = palette(),
##     pch = c(1, 0, 2, 5, 6, 3, 4), lty = 1:4, lwd = 2, alpha = NA,
##     border = palette(), border.alpha = NA, xtransf = NULL,
##     ytransf = NULL, add = FALSE) {
##     ## palette(khroma::colour('bright')())
##     ## scale_colour_discrete <- khroma::scale_colour_bright
##     ## if (missing(x)) {
##     ##     if (missing(y))
##     ##         stop("must specify at least one of 'x' and 'y'")
##     ##     else x <- seq_len(NROW(y))
##     ## }
##     ## else if (missing(y)) {
##     ##     if (missing(x))
##     ##         stop("must specify at least one of 'x' and 'y'")
##     ##     else y <- seq_len(NROW(x))
##     ## }
##     if (!missing(y) && !missing(x)) {
##         if (!is.list(x)) {
##             x <- apply(cbind(x), 2, identity, simplify = 'list')
##         }
##         if (!is.list(y)) {
##             y <- apply(cbind(y), 2, identity, simplify = 'list')
##         }
##     }
##     ##
##     else if (missing(x) && !missing(y)) {
##         if (!is.list(y)) {
##             y <- apply(cbind(y), 2, identity, simplify = 'list')
##         }
##         x <- lapply(y, seq_along)
##     } else if (missing(y) && !missing(x)) {
##         if (!is.list(x)) {
##             x <- apply(cbind(x), 2, identity, simplify = 'list')
##         }
##         y <- lapply(x, seq_along)
##     }
##     ##
##     xx <- unlist(x)
##     yy <- unlist(y)
##     if (!is.character(xx)) {
##         temp <- unique(xx[is.finite(xx)])
##         if(length(temp) > 1) {
##             xlim0 <- range(temp)
##         } else if(length(temp) == 1){
##             xlim0 <- range(temp) + c(-1, 1)
##         } else {
##             xlim0 <- c(0, 1)
##         }
##     } else {
##         uxx <- unique(xx)
##         if (is.character(xlabels) && all(uxx %in% xlabels)) {
##             uxx <- intersect(xlabels, uxx)
##         }
##         if (any(type == 'h')) {
##             xlim0 <- c(0.5, length(uxx) + 0.5)
##         } else {
##             xlim0 <- c(1, length(uxx))
##         }
##     }
##     if (!is.character(yy)) {
##         temp <- unique(yy[is.finite(yy)])
##         if(length(temp) > 1) {
##             ylim0 <- range(temp)
##         } else if(length(temp) == 1){
##             ylim0 <- range(temp) + c(-1, 1)
##         } else {
##             ylim0 <- c(0, 1)
##         }
##     } else {
##         uyy <- unique(yy)
##         if (is.character(ylabels) && all(uyy %in% ylabels)) {
##             uyy <- intersect(ylabels, uyy)
##         }
##         if (any(type == 'h')) {
##             ylim0 <- c(0.5, length(uyy) + 0.5)
##         } else {
##             ylim0 <- c(1, length(uyy))
##         }
##     }
##     if (is.na(ylim[1]) & any(type == 'h')) {
##         ylim[1] <- 0
##     }
##     xlim[is.na(xlim)] <- xlim0[is.na(xlim)]
##     ylim[is.na(ylim)] <- ylim0[is.na(ylim)]
##     if (length(n) < 2) {
##         n <- rep(n, 2)
##     }
##     if (is.null(xticks)) {
##         if (!is.character(xx)) {
##             xticks <- pretty(xlim, n = n[1])
##         } else {
##             xticks <- seq_along(uxx)
##         }
##     }
##     if (is.character(xx) && length(xlabels) == 1 && xlabels == TRUE) {
##         xlabels <- uxx
##     }
##     ## if(length(xlabels)==1 && xlabels){
##     ##     xlabels <- xticks
##     ##     xlabels[!(xticks %in% pretty(xlim, n=round(n[1]/2)))] <- NA
##     ## }
##     if (is.null(yticks)) {
##         if (!is.character(yy)) {
##             yticks <- pretty(ylim, n = n[2])
##         } else {
##             yticks <- seq_along(uyy)
##         }
##     }
##     if (is.character(yy) && length(ylabels) == 1 && ylabels == TRUE) {
##         ylabels <- uyy
##     }
##     ## if(length(ylabels)==1 && ylabels){
##     ##     ylabels <- yticks
##     ##     ylabels[!(yticks %in% pretty(ylim, n=round(n[2]/2)))] <- NA
##     ## }
##     if (is.null(lx)) {
##         lx <- 0
##     }
##     if (is.null(ly)) {
##         ly <- 1
##         if (!(length(yticks) == 1 && (any(is.na(yticks) | yticks == FALSE)))) {
##             ly <- 1 + max(nchar(sprintf('%.7g', yticks))) * 0.75
##         } else if (length(ylabels) > 1) {
##             ly <- 1 + max(nchar(ylabels)) * 0.75
##         }
##     }
##     if (is.null(xlab)) {
##         xlab <- names(x)[1]
##     }
##     if (is.null(ylab)) {
##         ylab <- names(y)[1]
##     }
##     ##
##     if (!add) {
##         plot.new()
##         ## par(mai=c(2, 3.5, 2, 0)/2.54, family='Palatino')#, mar=c(4,6,4,0)+0.1)
##         if (is.null(main)) {
##             marup <- 0
##         } else {
##             marup <- 3.5
##         }
##         if (is.null(mar)) {
##             mar <- c(3.25, ly, marup, 1) + c(1, 1.5, 1, 1)
##         }
##         mar[is.na(mar)] <- (c(3.25, ly, marup, 1) + c(1, 1.5, 1, 1))[is.na(mar)]
##         par(mar = mar, family = family) # , mar=c(4,6,4,0)+0.1)
##         ##
##         plot.window(xlim = xlim, ylim = ylim, xaxs = 'r', yaxs = 'r', asp = asp)
##         ##
##         if (!is.null(xtransf)) {
##             xlabels <- xtransf(xticks)
##         }
##         if (!is.null(ytransf)) {
##             ylabels <- ytransf(yticks)
##         }
##         if (!(length(xticks) == 1 & (any(is.na(xticks) | xticks == FALSE)))) {
##             axis(side = 1, at = xticks, labels = xlabels, tick = TRUE, lty = lty.axis,
##                 lwd = lwd.axis, lwd.ticks = lwd.ticks, col.ticks = col.ticks,
##                 gap.axis = NA, cex.axis = cex.axis, line = 0)
##         }
##         if (!(length(yticks) == 1 & (any(is.na(yticks) | yticks == FALSE)))) {
##             axis(side = 2, at = yticks, labels = ylabels, tick = TRUE, lty = lty.axis,
##                 lwd = lwd.axis, lwd.ticks = lwd.ticks, col.ticks = col.ticks,
##                 gap.axis = NA,las = las.y, cex.axis = cex.axis, line = 0)
##         }
##         ##
##         if (length(cex.lab) == 1) {
##             cex.lab <- rep(cex.lab, 2)
##         }
##         if (length(col.lab) == 1) {
##             col.lab <- rep(col.lab, 2)
##         }
##         if (!is.null(main)) {
##             title(main = main, cex.main = cex.main, line = 3)
##         }
##         if (!is.null(xlab)) {
##             title(xlab = xlab, cex.lab = cex.lab[1], line = 3 + lx,
##                 col.lab = col.lab[1])
##         }
##         if (!is.null(ylab)) {
##             title(ylab = ylab, cex.lab = cex.lab[2], line = ly, col.lab = col.lab[2])
##         }
##     }
##     if (is.null(xgrid)) {
##         xgrid <- !add
##     }
##     if (is.null(ygrid)) {
##         ygrid <- !add
##     }
##     if (xgrid) {
##         for (i in xticks) {
##             abline(v = i, lty = lty.axis, lwd = lwd.ticks, col = col.ticks)
##         }
##     }
##     if (ygrid) {
##         for (i in yticks) {
##             abline(h = i, lty = lty.axis, lwd = lwd.ticks, col = col.ticks)
##         }
##     }
##     ##
##     ## col[!grepl('^#',col)] <- palette()[as.numeric(col[!grepl('^#',col)])]
##     ## border[!grepl('^#',border)] <- palette()[as.numeric(border[!grepl('^#',border)])]
##     if (is.numeric(col)) {
##         col <- palette()[col]
##     }
##     if (is.numeric(border)) {
##         col <- palette()[border]
##     }
##     ##
##     nx <- length(x)
##     ny <- length(y)
##     if (all(!is.na(x) & !is.na(y))) {
##         for (j in 1:max(nx, ny)) {
##             xx <- x[[(j - 1) %% nx + 1]]
##             if (is.character(xx)) {
##                 xx <- match(xx, uxx)
##             }
##             yy <- y[[(j - 1) %% ny + 1]]
##             if (is.character(yy)) {
##                 yy <- match(yy, uyy)
##             }
##             if (length(xx) > length(yy) + 1 || length(yy) > length(xx) + 1) {
##                 stop(paste0('plot ', j, ': "x" and "y" must have same number ',
##                     'of rows or differ by 1'))
##             }
##             ialpha <- alpha[(j - 1) %% length(alpha) + 1]
##             icol <- col[(j - 1) %% length(col) + 1]
##             if (!(type[(j - 1) %% length(type) + 1] == 'h' ||
##                       length(xx) == length(yy) + 1 ||
##                       length(yy) == length(xx) + 1)) { # not a histogram
##                 if (is.na(ialpha)) {
##                     ialpha <- 1
##                 }
##                 icol <- adjustcolor(icol, alpha.f = ialpha)
##                 ## else if(!is.character(ialpha)){ialpha <- alpha2hex(ialpha)}
##                 ## if(!(is.na(icol) | nchar(icol)>7)){icol <- paste0(icol, ialpha)}
##                 ##
##                 plot.xy(xy.coords(x = xx, y = yy),
##                     type = type[(j - 1) %% length(type) + 1],
##                     col = icol,
##                     pch = pch[[(j - 1) %% length(pch) + 1]],
##                     lty = lty[(j - 1) %% length(lty) + 1],
##                     lwd = lwd[(j - 1) %% length(lwd) + 1],
##                     cex = cex[(j - 1) %% length(cex) + 1]
##                 )
##             } else { # histogram
##                 iborder <- border[(j - 1) %% length(border) + 1]
##                 iborder.alpha <- border.alpha[(j - 1) %% length(border.alpha) + 1]
##                 ##
##                 if (is.na(ialpha)) {
##                     ialpha <- 0.5
##                 }
##                 icol <- adjustcolor(icol, alpha.f = ialpha)
##                 ##
##                 if (is.na(iborder.alpha)) {
##                     iborder.alpha <- 0.5
##                 }
##                 iborder <- adjustcolor(iborder, alpha.f = iborder.alpha)
##                 ##
##                 if (length(yy) == length(xx)) {
##                     xx <- c(
##                     (3 * xx[1] - xx[2]) / 2,
##                     xx[-length(xx)] + diff(xx) / 2,
##                     (3 * xx[length(xx)] - xx[length(xx) - 1]) / 2
##                     )
##                 }
##                 if (length(xx) == length(yy) + 1) {
##                     ## if(is.na(ialpha)){ialpha <- '80'}
##                     ## else if(!is.character(ialpha)){ialpha <- alpha2hex(ialpha)}
##                     ## if(!(is.na(icol) | nchar(icol)>7)){icol <- paste0(icol, ialpha)}
##                     ## if(is.na(iborder.alpha)){iborder.alpha <- '80'}
##                     ## else if(!is.character(iborder.alpha)){iborder.alpha <- alpha2hex(iborder.alpha)}
##                     ## if(!(is.na(iborder) | nchar(iborder)>7)){iborder <- paste0(iborder, iborder.alpha)}
##                     for (i in 1:(length(xx) - 1)) {
##                         polygon(
##                             x = rbind(xx[i], xx[i], xx[i + 1], xx[i + 1]),
##                             y = rbind(ylim[1], yy[i], yy[i], ylim[1]),
##                             col = icol, border = iborder,
##                             lty = lty[(j - 1) %% length(lty) + 1],
##                             lwd = lwd[(j - 1) %% length(lwd) + 1]
##                         )
##                     }
##                 } else {
##                     ## iborder <- border[(j-1)%%length(border)+1]
##                     ## iborder.alpha <- border.alpha[(j-1)%%length(border.alpha)+1]
##                     ##
##                     ## if(is.na(ialpha)){ialpha <- 0.5}
##                     ## icol <- alpha2hex(icol, ialpha)
##                     ## if(is.na(ialpha)){ialpha <- '80'}
##                     ## else if(!is.character(ialpha)){ialpha <- alpha2hex(ialpha)}
##                     ## if(!(is.na(icol) | nchar(icol)>7)){icol <- paste0(icol, ialpha)}
##                     ##
##                     ## if(is.na(iborder.alpha)){iborder.alpha <- 0.5}
##                     ## iborder <- alpha2hex(iborder, iborder.alpha)
##                     ## if(is.na(iborder.alpha)){iborder.alpha <- '80'}
##                     ## else if(!is.character(iborder.alpha)){iborder.alpha <- alpha2hex(iborder.alpha)}
##                     ## if(!(is.na(iborder) | nchar(iborder)>7)){iborder <- paste0(iborder, iborder.alpha)}
##                     for (i in 1:(length(yy) - 1)) {
##                         polygon(
##                             y = rbind(yy[i], yy[i], yy[i + 1], yy[i + 1]),
##                             x = rbind(xlim[1], xx[i], xx[i], xlim[1]),
##                             col = icol, border = iborder,
##                             lty = lty[(j - 1) %% length(lty) + 1],
##                             lwd = lwd[(j - 1) %% length(lwd) + 1]
##                         )
##                     }
##                 }
##             }
##         }
##     }
## }

## #' @keywords internal
## tlegend <- function(x, y=NULL, legend, col=palette(), pch=c(1,0,2,5,6,3,4), lty=1:4, lwd=2, alpha=1, cex=1.5, ...){
##     suppressWarnings(col <- mapply(function(i,j)adjustcolor(i,j),col,alpha))
##     legend(x=x, y=y, legend=legend, col=col, pch=pch, lty=lty, lwd=lwd, bty='n', cex=cex, ...)
## }
##
## #' @keywords internal
## fivenumaxis <- function(side, x, col='#555555', type=6){
##     x <- x[!is.na(x) && is.finite(x)]
##     if(length(x)==0){x <- c(0,1)}
##     if(diff(range(x))==0){x <- range(x) + c(-1,1)}
##     ylim <- par('usr')
##     five <- c(min(x), quantile(x=x, probs=(1:3)/4, type=type), max(x))
##     ##
##     if(side==1){
##         xl <- rbind(five[c(1,4)], five[c(2,5)])
##         yl <- rep(ylim[3],2)
##         xp <- five[3]
##         yp <- ylim[3]
##     }else if(side==2){
##         yl <- rbind(five[c(1,4)], five[c(2,5)])
##         xl <- rep(ylim[1],2)
##         yp <- five[3]
##         xp <- ylim[1]
##     }else if(side==3){
##         xl <- rbind(five[c(1,4)], five[c(2,5)])
##         yl <- rep(ylim[2],2)
##         xp <- five[3]
##         yp <- ylim[2]
##     }else if(side==4){
##         yl <- rbind(five[c(1,4)], five[c(2,5)])
##         xl <- rep(ylim[4],2)
##         yp <- five[3]
##         xp <- ylim[4]
##     }
##     ##
##     matlines(x=xl, y=yl, lty=1, lwd=2, col=col)
##     matpoints(x=xp, y=yp, pch=18, cex=2, col=col)
## }

## #' Construct and plot histograms
## #'
## #' Utility function to construct and plot histograms, from samples obtained with \code{\link{Pr}}. The advantage of this function over \code{\link[base]{hist}} is that its `x` argument can be a vector of character values, as typically the case for nominal or ordinal variates (see \code{\link{metadata}}). The plot will have these characters as labels of the x-axis.
## #'
## #' @param x Numeric or character: a collection of data or simulated values.
## #' @param n integer or vector: number of bins to use, or boundaries of bins.
## #' @param plot logical: plot the results? Default `FALSE`.
## #'
## #' @return A list of several elements. Most important: `counts`: the value count in each bin (frequency). `density`: the value count dividen by the bin width (frequency density). `mids`: the values of the bins' centres. `breaks`: the boundaries of the bins.
## #'
## #' @export
## thist <- function(x, n=NULL, type=6, pretty=FALSE, plot=FALSE, extendbreaks=FALSE, ...){
##     if(!is.list(x)){x <- list(x)}
##     if(!is.list(n)){n <- list(n)}
##     out <- list()
##     for(i in 1:length(x)){
##         ax <- x[[i]]
##         an <- n[[(i-1)%%length(n)+1]]
##         if(is.character(ax)){
##             nextout <- c(table(ax[!is.na(ax)]))
##             nextout <- list(
##                 breaks=NA,
##                 counts=unname(nextout),
##                 density=unname(nextout)/sum(nextout),
##                 mids=names(nextout),
##                 xname=names(x)[i],
##                 equidist=NA
##             )
##         }else{
##             ax <- ax[!is.na(ax) & is.finite(ax)]
##             if(is.null(an)){an <- (round(sqrt(length(ax))/2))}
##             if(length(an)==1 && (is.na(an) || an=='i' || an=='integer')){breaks <- (round(min(ax))-0.5):(round(max(ax))+0.5)}
##             else if(length(an) > 1 || is.character(an)){breaks <- an}
##             else if(length(an) == 1 && an > 0){
##                 rg <- range(ax)
##                 if(diff(rg)==0){rg <- rg + c(-0.5,0.5)}
##                 breaks <- seq(rg[1], rg[2], length.out=an+1)}
##             else if(length(an) == 1 && an < 0){
##                     rg <- range(ax)
##                     if(diff(rg)==0){rg <- rg + c(-0.5,0.5)}
##                     breaks <- seq(rg[1], rg[2]-an, by=-an)
##                     breaks <- breaks - (breaks[length(breaks)]-rg[2])/2
##                 }
##             else {print('Error with n')}
##             if(!is.null(pretty) && pretty){
##                 breaks <- pretty(ax, n=length(breaks)-1)
##             }
##             if(extendbreaks){
##                 breaks <- c(-Inf,breaks,+Inf)
##             }
##             nextout <- hist(x=ax, breaks=breaks, plot=FALSE)
##         }
##         out <- c(out,list(nextout))
##     }
##     if(plot){
##         tplot(x=lapply(out,function(xx){
##             if(length(xx$breaks)==1){xx$mids}else{xx$breaks}
##         } ),
##         y=lapply(out,function(xx)xx$density),ylim=c(0,NA),type='h', ...)
##     }else{
##         if(length(out)==1){unlist(out,recursive=F)}else{out}
##     }
## }

## tquant <- function(x, probs=c(1:3)/4, na.rm=TRUE, names=TRUE, type=6, ...){
##     quantile(x=x, probs=probs, na.rm=na.rm, names=names, type=type, ...)
## }
##
## tmad <- function(x){mad(x, constant=1, na.rm=TRUE)}

## tsummary <- function(x){
##     x <- cbind(x)
##     apply(x, 2, function(xx){
##         c(tquant(xx, c(2.5/100,1/8,2/8,4/8,6/8,7/8,97.5/100)), MAD=mad(xx,constant=1,na.rm=T), IQR=IQR(xx,na.rm=T), mean=mean(xx,na.rm=T), sd=sd(xx,na.rm=T), hr=diff(range(xx,na.rm=T))/2, min=min(xx,na.rm=T), max=max(xx,na.rm=T), NAs=sum(is.na(xx)))
##     })
## }

## Function to build powerset
## powerset <<- function(set){
##     n <- length(set)
##     masks <- 2^(1:n-1)
##     lapply( 1:2^n-1, function(u) set[ bitwAnd(u, masks) != 0 ] )
## }

## Greatest common denominator
## gcd <<- function(...){Reduce(function(a, b){if (b == 0) a else Recall(b, a %% b)}, c(...))}

## Normalize according to row
## tnormalize <- function(x){
##     if(is.null(dim(x)) || is.table(x)){
##         x/sum(x,na.rm=T)
##     }else{
##         aperm(aperm(x)/c(aperm(cbind(colSums(x,na.rm=T)))))
##     }
## }
##
## ## Table with list of values
## tablev <- function(x, values=NULL, norm=FALSE){
##     if(norm){
##         tnormalize(table(c(x,values))-!(is.null(values)))
##     }else{
##         table(c(x,values))-!(is.null(values))
##     }
## }
