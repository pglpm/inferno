suppressPackageStartupMessages(require(grDevices))
options(bitmapType='cairo')
## Colour-blind friendly palettes, from https://personal.sron.nl/~pault/
library('khroma')
palette(colour('bright')())
bluepurple <- palette()[1]
red <- palette()[2]
green <- palette()[3]
yellow <- palette()[4]
blue <- palette()[5]
redpurple <- palette()[6]
grey <- palette()[7]
darkgrey <- '#555555'
black <- '#000000'
scale_colour_discrete <- scale_colour_bright

pdff <- function(file='Rplot', height=210/25.4, width=297/25.4, asp=NA, ...){
    if(!is.na(asp)){
        width <- height*asp
    }
    pdf(file=paste0(file,'.pdf'), paper='special', height=height, width=width, ...)
} # to output in pdf format
pngf <- function(filename='Rplot',res=300){png(file=paste0(filename,'.png'), height=11.7*1.2, width=16.5*1.2, units='in', res=res, pointsize=36)} # to output in png format

alpha2hex <- function(alpha){
    if(!is.character(alpha)){alpha <- sprintf('%02x', round((1-alpha)*255))}
    alpha
}

tticks <- pretty
## tticks <- function(x, n=10){
##     x <- x[!is.na(x) && is.finite(x)]
##     if(length(x)==0){x <- c(0,1)}
##     rg <- range(x, na.rm=T)
##     if(diff(rg)==0){rg <- rg + c(-1,1)}
##     ext <- diff(rg)
##     deltas <- sort(c(c(1,2,5) * 10^ceiling(log10(ext/(n*c(1,2,5)))),
##                 c(1,2,5) * 10^floor(log10(ext/(n*c(1,2,5))))))
##     ## print(deltas)
##     ## print(abs(sapply(deltas, function(d){
##     ##     length((rg[1]%/%d):((rg[2]%/%d) + (rg[2]%%d > 0)))}) -n
##     ## ))
##     ##deltas <- c(0.1, 0.2, 0.5, 1) * n^floor(log(ext,base=n))
##     delta <- deltas[which.min(abs(sapply(deltas, function(d){
##         diff(rg%/%d) + (rg[2]%%d > 0)}) + 1 - n
##     ))]
##     ## delta <- deltas[which.min(abs(sapply(deltas, function(d){
##     ##     length((rg[1]%/%d):((rg[2]%/%d) + (rg[2]%%d > 0)))}) -n
##     ## ))]
##     ((rg[1]%/%delta):((rg[2]%/%delta) + (rg[2]%%delta > 0)))*delta
## }

tplot <- function(x, y, xlim=c(NA,NA), ylim=c(NA,NA), asp=NA, n=10, family='', xticks=NULL, xlabels=TRUE, yticks=NULL, ylabels=TRUE, cex=1.5, ly=NULL, lx=NULL, mar=NULL, lty.axis=1, lwd.axis=0, lwd.ticks=1, col.ticks='#bbbbbb80', col.lab='black', cex.axis=1.5, las.y=1, xgrid=NULL, ygrid=NULL, main=NULL, cex.main=1.5, xlab=NULL, ylab=NULL, cex.lab=1.5, type='l', col=palette(), pch=c(1,0,2,5,6,3,4), lty=1:4, lwd=2, alpha=NA, border=palette(), border.alpha=NA, add=FALSE){
    ## if (missing(x)) {
    ##     if (missing(y)) 
    ##         stop("must specify at least one of 'x' and 'y'")
    ##     else x <- seq_len(NROW(y))
    ## }
    ## else if (missing(y)) {
    ##     if (missing(x)) 
    ##         stop("must specify at least one of 'x' and 'y'")
    ##     else y <- seq_len(NROW(x))
    ## }
    if(!missing(y) & !missing(x)){
        if(!is.list(x)){
            x <- apply(cbind(x), 2, identity, simplify='list')
        }
        if(!is.list(y)){
            y <- apply(cbind(y), 2, identity, simplify='list')
        }
    }
    ##
    else if(missing(x) & !missing(y)){
        if(!is.list(y)){
            y <- apply(cbind(y), 2, identity, simplify='list')
        }
        x <- lapply(y,seq_along)
    }
    else if(missing(y) & !missing(x)){
        if(!is.list(x)){
            x <- apply(cbind(x), 2, identity, simplify='list')
        }
        y <- lapply(x,seq_along)
    }
    ##
    xx <- unlist(x)
    yy <- unlist(y)
    xlim0 <- range(xx[is.finite(xx)], na.rm=TRUE)
    ylim0 <- range(yy[is.finite(yy)], na.rm=TRUE)
    xlim[is.na(xlim)] <- xlim0[is.na(xlim)]
    ylim[is.na(ylim)] <- ylim0[is.na(ylim)]
    if(length(n)<2){n <- rep(n,2)}
    if(is.null(xticks)){xticks <- pretty(xlim, n=n[1])}
    if(is.null(yticks)){yticks <- pretty(ylim, n=n[2])}
    if(is.null(lx)){lx <- 0}
    if(is.null(ly)){
        ly <- 1
        if(!(length(yticks)==1 && (any(is.na(yticks) || yticks==FALSE))) && ylabels==TRUE){
            ly <- 1+max(nchar(sprintf('%.7g',yticks)))*0.75
        } else if(length(ylabels)>1){
            ly <- 1+max(nchar(ylabels))*0.75
        }
    }
    if(is.null(xlab)){xlab <- names(x)[1]}
    if(is.null(ylab)){ylab <- names(y)[1]}
    ##
    if(!add){plot.new()
    ##par(mai=c(2, 3.5, 2, 0)/2.54, family='Palatino')#, mar=c(4,6,4,0)+0.1)
    if(is.null(mar)){mar <- c(3.25, ly, 3, 1)+c(1,1.1,1,1)}
    par(mar=mar, family=family)#, mar=c(4,6,4,0)+0.1)
    ##
    plot.window(xlim=xlim, ylim=ylim, xaxs='r', yaxs='r', asp=asp)
    ##
        if(!(length(xticks)==1 && (any(is.na(xticks) || xticks==FALSE)))){ axis(side=1, at=xticks, labels=xlabels, tick=TRUE, lty=lty.axis, lwd=lwd.axis, lwd.ticks=lwd.ticks, col.ticks=col.ticks, gap.axis=0.25, cex.axis=cex.axis, line=0)}
        if(!(length(yticks)==1 && (any(is.na(yticks) || yticks==FALSE)))){ axis(side=2, at=yticks, labels=ylabels, tick=TRUE, lty=lty.axis, lwd=lwd.axis, lwd.ticks=lwd.ticks, col.ticks=col.ticks, las=las.y, cex.axis=cex.axis, line=0)}
        ##
        if(length(cex.lab)==1){cex.lab <- rep(cex.lab, 2)}
        if(length(col.lab)==1){col.lab <- rep(col.lab, 2)}
    if(!is.null(main)){title(main=main, cex.main=cex.main, line=3)}
    if(!is.null(xlab)){title(xlab=xlab, cex.lab=cex.lab[1], line=3+lx, col.lab=col.lab[1])}
    if(!is.null(ylab)){title(ylab=ylab, cex.lab=cex.lab[2], line=ly, col.lab=col.lab[2])}
    }
    if(is.null(xgrid)){xgrid <- !add}
    if(is.null(ygrid)){ygrid <- !add}
    if(xgrid){for(i in xticks){abline(v=i, lty=lty.axis, lwd=lwd.ticks, col=col.ticks)}}
    if(ygrid){for(i in yticks){abline(h=i, lty=lty.axis, lwd=lwd.ticks, col=col.ticks)}}
    ##
    col[!grepl('^#',col)] <- palette()[as.numeric(col[!grepl('^#',col)])]
    border[!grepl('^#',border)] <- palette()[as.numeric(border[!grepl('^#',border)])]
    ##
    nx <- length(x)
    ny <- length(y)
    if(all(!is.na(x) & !is.na(y))){
        for(j in 1:max(nx,ny)){
            xx <- x[[(j-1)%%nx+1]]
            yy <- y[[(j-1)%%ny+1]]
            if(length(xx)>length(yy)+1 || length(yy)>length(xx)+1){stop(paste0("plot ",j,": 'x' and 'y' must have same number of rows or differ by 1"))}
            ialpha <- alpha[(j-1)%%length(alpha)+1]
            icol <- col[(j-1)%%length(col)+1]
            if(!(type[(j-1)%%length(type)+1]=='h' || length(xx)==length(yy)+1 || length(yy)==length(xx)+1)){# not a histogram
                if(is.na(ialpha)){ialpha <- ''}
                else if(!is.character(ialpha)){ialpha <- alpha2hex(ialpha)}
                if(!(is.na(icol) | nchar(icol)>7)){icol <- paste0(icol, ialpha)}
                ##
                plot.xy(xy.coords(x=xx, y=yy),
                        type=type[(j-1)%%length(type)+1],
                        col=icol,
                        pch=pch[(j-1)%%length(pch)+1], lty=lty[(j-1)%%length(lty)+1], lwd=lwd[(j-1)%%length(lwd)+1], cex=cex[(j-1)%%length(cex)+1])
            }else{# histogram
                if(length(xx)==length(yy)+1){
                iborder <- border[(j-1)%%length(border)+1]
                iborder.alpha <- border.alpha[(j-1)%%length(border.alpha)+1]
                ##
                if(is.na(ialpha)){ialpha <- '80'}
                else if(!is.character(ialpha)){ialpha <- alpha2hex(ialpha)}
                if(!(is.na(icol) | nchar(icol)>7)){icol <- paste0(icol, ialpha)}
                ##
                if(is.na(iborder.alpha)){iborder.alpha <- '80'}
                else if(!is.character(iborder.alpha)){iborder.alpha <- alpha2hex(iborder.alpha)}
                if(!(is.na(iborder) | nchar(iborder)>7)){iborder <- paste0(iborder, iborder.alpha)}
                for(i in 1:(length(xx)-1)){
                    polygon(x=rbind(xx[i], xx[i], xx[i+1], xx[i+1]),
                            y=rbind(ylim[1], yy[i], yy[i], ylim[1]),
                            col=icol, border=iborder,
                            lty=lty[(j-1)%%length(lty)+1],
                            lwd=lwd[(j-1)%%length(lwd)+1])
                }
                }else{
                iborder <- border[(j-1)%%length(border)+1]
                iborder.alpha <- border.alpha[(j-1)%%length(border.alpha)+1]
                ##
                if(is.na(ialpha)){ialpha <- '80'}
                else if(!is.character(ialpha)){ialpha <- alpha2hex(ialpha)}
                if(!(is.na(icol) | nchar(icol)>7)){icol <- paste0(icol, ialpha)}
                ##
                if(is.na(iborder.alpha)){iborder.alpha <- '80'}
                else if(!is.character(iborder.alpha)){iborder.alpha <- alpha2hex(iborder.alpha)}
                if(!(is.na(iborder) | nchar(iborder)>7)){iborder <- paste0(iborder, iborder.alpha)}
                for(i in 1:(length(yy)-1)){
                    polygon(y=rbind(yy[i], yy[i], yy[i+1], yy[i+1]),
                            x=rbind(xlim[1], xx[i], xx[i], xlim[1]),
                            col=icol, border=iborder,
                            lty=lty[(j-1)%%length(lty)+1],
                            lwd=lwd[(j-1)%%length(lwd)+1])
                }
                }
            }
        }
    }
}


fivenumaxis <- function(side, x, col='#555555', type=8){
    x <- x[!is.na(x) && is.finite(x)]
    if(length(x)==0){x <- c(0,1)}
    if(diff(range(x))==0){x <- range(x) + c(-1,1)}
    ylim <- par('usr')
    five <- c(min(x), quantile(x=x, probs=(1:3)/4, type=type), max(x))
    ##
    if(side==1){
        xl <- rbind(five[c(1,4)], five[c(2,5)])
        yl <- rep(ylim[3],2)
        xp <- five[3]
        yp <- ylim[3]
    }else if(side==2){
        yl <- rbind(five[c(1,4)], five[c(2,5)])
        xl <- rep(ylim[1],2)
        yp <- five[3]
        xp <- ylim[1]
    }else if(side==3){
        xl <- rbind(five[c(1,4)], five[c(2,5)])
        yl <- rep(ylim[2],2)
        xp <- five[3]
        yp <- ylim[2]
    }else if(side==4){
        yl <- rbind(five[c(1,4)], five[c(2,5)])
        xl <- rep(ylim[4],2)
        yp <- five[3]
        xp <- ylim[4]
    }
    ##
    matlines(x=xl, y=yl, lty=1, lwd=2, col=col)
    matpoints(x=xp, y=yp, pch=18, cex=2, col=col)
}

scatteraxis <- function(side, x, n=256, col='#555555', alpha='88', lwd=0.1, ...){
    x <- x[!is.na(x) && is.finite(x)]
    if(is.na(n)){n <- length(x)}
    x <- x[round(seq(1, length(x), length.out=n))]
    ylim <- par('usr')
    exts <- diff(ylim)[-2]
    ##
    if(side==1){
        xl <- rbind(x, x)
        yl <- rbind(rep(ylim[3],length(x))+2*exts[2]/100,rep(ylim[3],length(x))+3*exts[2]/100)
    }else if(side==2){
        yl <- rbind(x, x)
        xl <- rbind(rep(ylim[1],length(x))+2*exts[1]/100,rep(ylim[1],length(x))+3*exts[1]/100)
    }else if(side==3){
        xl <- rbind(x, x)
        yl <- rbind(rep(ylim[4],length(x))-2*exts[2]/100,rep(ylim[4],length(x))-3*exts[2]/100)
    }else if(side==4){
        yl <- rbind(x, x)
        xl <- rbind(rep(ylim[2],length(x))-2*exts[1]/100,rep(ylim[2],length(x))-3*exts[1]/100)
    }
    ##
    if(sum(!is.na(col))>0 && !grepl('^#', col[!is.na(col)])){col[!is.na(col)] <- palette()[col[!is.na(col)]]}
    if(!is.null(alpha) && !is.character(alpha)){alpha <- alpha2hex(alpha)}
    col[!is.na(col)] <- paste0(col[!is.na(col)], alpha)
    matlines(x=xl, y=yl, lty=1, lwd=lwd, col=col, ...)
}

thist <- function(x, n=NULL, type=8){
    x <- x[!is.na(x) & is.finite(x)]
    if(is.null(n)){n <- (round(sqrt(length(x))/2))}
    if(length(n) > 1 | is.character(n)){
        breaks <- n
        x <- x[x>=min(breaks)&x<=max(breaks)]
    }
    else if(is.na(n) || n=='i' || n=='integer'){breaks <- (round(min(x))-0.5):(round(max(x))+0.5)}
    else if(length(n) == 1 && n<0){
        rg <- range(x)
        if(diff(rg)==0){rg <- rg + c(-0.5,0.5)}
        breaks <- seq(rg[1], rg[2], length.out=n+1)}
    else {breaks <- pretty(x, n=n)}
    hist(x=x, breaks=breaks, plot=FALSE)
}

quant <- function(x, probs=c(1:3)/4, na.rm=TRUE, names=TRUE, type=8, ...){
    quantile(x=x, probs=probs, na.rm=na.rm, names=names, type=type, ...)
}
