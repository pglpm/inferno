suppressPackageStartupMessages(require(grDevices))
options(bitmapType='cairo')
## Colour-blind friendly palettes, from https://personal.sron.nl/~pault/
library('khroma')
## palette(colour('bright')())
cc <- colour('bright')()
cc[8] <- '#000000'
names(cc)[8] <- 'black'
palette(cc)
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

pdff <- function(file='Rplot', apaper=5, portrait=FALSE, height=148/25.4, width=210/25.4, asp=NA, ...){
    if(is.numeric(apaper)){
        if(portrait){
            height <- floor(841/sqrt(2)^(apaper-1))/25.4
            width <- floor(841/sqrt(2)^(apaper))/25.4
        }else{
            width <- floor(841/sqrt(2)^(apaper-1))/25.4
            height <- floor(841/sqrt(2)^(apaper))/25.4
        }
    }
    if(!is.na(asp)){
        width <- height*asp
    }
    pdf(file=paste0(file,'.pdf'), paper='special', height=height, width=width, ...)
} # to output in pdf format
pngf <- function(filename='Rplot',res=300){png(file=paste0(filename,'.png'), height=11.7*1.2, width=16.5*1.2, units='in', res=res, pointsize=36)} # to output in png format

alpha2hex2 <- function(alpha,col=NULL){
    if(!is.character(alpha)){alpha <- sprintf('%02x', round((1-alpha)*255))}
    if(is.numeric(col)){col <- palette()[col]}
    paste0(col,alpha)
}

alpha2hex <- function(col, alpha=NULL){
    if(is.null(alpha)){alpha <- 0}
    if(!is.character(col)){col <- palette()[col]}
    do.call(rgb, c(as.list(col2rgb(col)),list((1-alpha)*255, maxColorValue=255)))
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

tplot <- function(x, y, xlim=c(NA,NA), ylim=c(NA,NA), asp=NA, n=10, family='', xticks=NULL, xlabels=TRUE, yticks=NULL, ylabels=TRUE, cex=1.5, ly=NULL, lx=NULL, mar=NULL, lty.axis=1, lwd.axis=0, lwd.ticks=1, col.ticks='#bbbbbb80', col.lab='black', cex.axis=1.35, las.y=1, xgrid=NULL, ygrid=NULL, main=NULL, cex.main=1.5, xlab=NULL, ylab=NULL, cex.lab=1.5, type='l', col=palette(), pch=c(1,0,2,5,6,3,4), lty=1:4, lwd=2, alpha=NA, border=palette(), border.alpha=NA, xtransf=NULL, ytransf=NULL, add=FALSE){
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
    if(!is.character(xx)){
        xlim0 <- range(xx[is.finite(xx)], na.rm=TRUE)
    }else{
        uxx <- unique(xx)
        if(any(type=='h')){
            xlim0 <- c(0.5,length(uxx)+0.5)
        }else{
            xlim0 <- c(1,length(uxx))
        }
    }
    if(!is.character(yy)){
        ylim0 <- range(yy[is.finite(yy)], na.rm=TRUE)
    }else{
        uyy <- unique(yy)
        if(any(type=='h')){
            ylim0 <- c(0.5,length(uyy)+0.5)
        }else{
            ylim0 <- c(1,length(uyy))
        }
    }
    if(is.na(ylim[1]) & any(type=='h')){ylim[1] <- 0}
    xlim[is.na(xlim)] <- xlim0[is.na(xlim)]
    ylim[is.na(ylim)] <- ylim0[is.na(ylim)]
    if(length(n)<2){n <- rep(n,2)}
    if(is.null(xticks)){
        if(!is.character(xx)){
            xticks <- pretty(xlim, n=n[1])
        }else{
            xticks <- 1:length(uxx)
        }
    }
    if(is.character(xx) && length(xlabels)==1 && xlabels==TRUE){
        xlabels <- uxx
    }
    ## if(length(xlabels)==1 && xlabels){
    ##     xlabels <- xticks
    ##     xlabels[!(xticks %in% pretty(xlim, n=round(n[1]/2)))] <- NA
    ## }
    if(is.null(yticks)){
        if(!is.character(yy)){
            yticks <- pretty(ylim, n=n[2])
        }else{
            yticks <- 1:length(uyy)
        }
    }
    if(is.character(yy) && length(ylabels)==1 && ylabels==TRUE){
        ylabels <- uyy
    }
    ## if(length(ylabels)==1 && ylabels){
    ##     ylabels <- yticks
    ##     ylabels[!(yticks %in% pretty(ylim, n=round(n[2]/2)))] <- NA
    ## }
    if(is.null(lx)){lx <- 0}
    if(is.null(ly)){
        ly <- 1
        if(!(length(yticks)==1 & (any(is.na(yticks) | yticks==FALSE)))){
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
        if(is.null(main)){marup <- 0}else{marup <- 3.5}
        if(is.null(mar)){mar <- c(3.25, ly, marup, 1)+c(1,1.1,1,1)}
        mar[is.na(mar)] <- (c(3.25, ly, marup, 1)+c(1,1.1,1,1))[is.na(mar)]
        par(mar=mar, family=family)#, mar=c(4,6,4,0)+0.1)
    ##
    plot.window(xlim=xlim, ylim=ylim, xaxs='r', yaxs='r', asp=asp)
        ##
        if(!is.null(xtransf)){xlabels <- xtransf(xticks)}
        if(!is.null(ytransf)){ylabels <- ytransf(yticks)}
        if(!(length(xticks)==1 & (any(is.na(xticks) | xticks==FALSE)))){ axis(side=1, at=xticks, labels=xlabels, tick=TRUE, lty=lty.axis, lwd=lwd.axis, lwd.ticks=lwd.ticks, col.ticks=col.ticks, gap.axis=NA, cex.axis=cex.axis, line=0)}
        if(!(length(yticks)==1 & (any(is.na(yticks) | yticks==FALSE)))){ axis(side=2, at=yticks, labels=ylabels, tick=TRUE, lty=lty.axis, lwd=lwd.axis, lwd.ticks=lwd.ticks, col.ticks=col.ticks, gap.axis=NA, las=las.y, cex.axis=cex.axis, line=0)}
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
    ## col[!grepl('^#',col)] <- palette()[as.numeric(col[!grepl('^#',col)])]
    ## border[!grepl('^#',border)] <- palette()[as.numeric(border[!grepl('^#',border)])]
    if(is.numeric(col)){col <- palette()[col]}
    if(is.numeric(border)){col <- palette()[border]}
    ##
    nx <- length(x)
    ny <- length(y)
    if(all(!is.na(x) & !is.na(y))){
        for(j in 1:max(nx,ny)){
            xx <- x[[(j-1)%%nx+1]]
            if(is.character(xx)){
                xx <- match(xx,uxx)
            }
            yy <- y[[(j-1)%%ny+1]]
            if(is.character(yy)){
                yy <- match(yy,uyy)
            }
            if(length(xx)>length(yy)+1 || length(yy)>length(xx)+1){stop(paste0("plot ",j,": 'x' and 'y' must have same number of rows or differ by 1"))}
            ialpha <- alpha[(j-1)%%length(alpha)+1]
            icol <- col[(j-1)%%length(col)+1]
            if(!(type[(j-1)%%length(type)+1]=='h' || length(xx)==length(yy)+1 || length(yy)==length(xx)+1)){# not a histogram
                if(is.na(ialpha)){ialpha <- 0}
                icol <- alpha2hex(icol, ialpha)
                ## else if(!is.character(ialpha)){ialpha <- alpha2hex(ialpha)}
                ## if(!(is.na(icol) | nchar(icol)>7)){icol <- paste0(icol, ialpha)}
                ##
                plot.xy(xy.coords(x=xx, y=yy),
                        type=type[(j-1)%%length(type)+1],
                        col=icol,
                        pch=pch[[(j-1)%%length(pch)+1]], lty=lty[(j-1)%%length(lty)+1], lwd=lwd[(j-1)%%length(lwd)+1], cex=cex[(j-1)%%length(cex)+1])
            }else{# histogram
                iborder <- border[(j-1)%%length(border)+1]
                iborder.alpha <- border.alpha[(j-1)%%length(border.alpha)+1]
                ##
                if(is.na(ialpha)){ialpha <- 0.5}
                icol <- alpha2hex(icol, ialpha)
                ##
                if(is.na(iborder.alpha)){iborder.alpha <- 0.5}
                iborder <- alpha2hex(iborder, iborder.alpha)
                ##
                if(length(yy)==length(xx)){
                    xx <- c((3*xx[1]-xx[2])/2,
                            xx[-length(xx)] + diff(xx)/2,
                            (3*xx[length(xx)]-xx[length(xx)-1])/2
                            )
                    }
                if(length(xx)==length(yy)+1){
                ## if(is.na(ialpha)){ialpha <- '80'}
                ## else if(!is.character(ialpha)){ialpha <- alpha2hex(ialpha)}
                ## if(!(is.na(icol) | nchar(icol)>7)){icol <- paste0(icol, ialpha)}
                ## if(is.na(iborder.alpha)){iborder.alpha <- '80'}
                ## else if(!is.character(iborder.alpha)){iborder.alpha <- alpha2hex(iborder.alpha)}
                ## if(!(is.na(iborder) | nchar(iborder)>7)){iborder <- paste0(iborder, iborder.alpha)}
                for(i in 1:(length(xx)-1)){
                    polygon(x=rbind(xx[i], xx[i], xx[i+1], xx[i+1]),
                            y=rbind(ylim[1], yy[i], yy[i], ylim[1]),
                            col=icol, border=iborder,
                            lty=lty[(j-1)%%length(lty)+1],
                            lwd=lwd[(j-1)%%length(lwd)+1])
                }
                }else{
                ## iborder <- border[(j-1)%%length(border)+1]
                ## iborder.alpha <- border.alpha[(j-1)%%length(border.alpha)+1]
                ##
                ## if(is.na(ialpha)){ialpha <- 0.5}
                ## icol <- alpha2hex(icol, ialpha)
                ## if(is.na(ialpha)){ialpha <- '80'}
                ## else if(!is.character(ialpha)){ialpha <- alpha2hex(ialpha)}
                ## if(!(is.na(icol) | nchar(icol)>7)){icol <- paste0(icol, ialpha)}
                ##
                ## if(is.na(iborder.alpha)){iborder.alpha <- 0.5}
                ## iborder <- alpha2hex(iborder, iborder.alpha)
                ## if(is.na(iborder.alpha)){iborder.alpha <- '80'}
                ## else if(!is.character(iborder.alpha)){iborder.alpha <- alpha2hex(iborder.alpha)}
                ## if(!(is.na(iborder) | nchar(iborder)>7)){iborder <- paste0(iborder, iborder.alpha)}
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

tlegend <- function(x, y=NULL, legend, col=palette(), pch=c(1,0,2,5,6,3,4), lty=1:4, lwd=2, alpha=0, cex=1.5, ...){
    suppressWarnings(col <- mapply(function(i,j)alpha2hex(i,j),col,alpha))
    legend(x=x, y=y, legend=legend, col=col, pch=pch, lty=lty, lwd=lwd, bty='n', ...)
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

plotquantiles <- function(x, y, col=7, alpha=0.75, border=NA){
    if(dim(y)[2]==2){y <- t(y)}
    ##
    ## col[!grepl('^#',col)] <- palette()[as.numeric(col[!grepl('^#',col)])]
    if(is.na(alpha)){alpha <- 0}
    col <- alpha2hex(col, alpha)
    ## if(is.na(alpha)){alpha <- ''}
    ## else if(!is.character(alpha)){alpha <- alpha2hex(alpha)}
    ## if(!(is.na(col) | nchar(col)>7)){col <- paste0(col, alpha)}
    ##
    polygon(x=c(x,rev(x)), y=c(y[1,], rev(y[2,])),
            col=col, border=border)
}

scatteraxis <- function(x, side=1, n=128, col='#555555', alpha=0.5, ext=5, pos=NULL, exts=NULL, lwd=0.1, ...){
    x <- x[!is.na(x) & is.finite(x)]
    if(is.na(n)){n <- length(x)}
    x <- x[round(seq(1, length(x), length.out=n))]
    if(is.null(pos)){ pos <- par('usr') }
    if(is.null(exts)){ exts <- diff(pos)[-2]/100}
    ##
    if(side==1){
        xl <- rbind(x, x)
        yl <- rbind(rep(pos[3],length(x))+2*exts[2],rep(pos[3],length(x))+ext*exts[2])
    }else if(side==2){
        yl <- rbind(x, x)
        xl <- rbind(rep(pos[1],length(x))+2*exts[1],rep(pos[1],length(x))+ext*exts[1])
    }else if(side==3){
        xl <- rbind(x, x)
        yl <- rbind(rep(pos[4],length(x))-2*exts[2],rep(pos[4],length(x))-ext*exts[2])
    }else if(side==4){
        yl <- rbind(x, x)
        xl <- rbind(rep(pos[2],length(x))-2*exts[1],rep(pos[2],length(x))-ext*exts[1])
    }
    ##
    if(sum(!is.na(col))>0 && !grepl('^#', col[!is.na(col)])){col[!is.na(col)] <- palette()[col[!is.na(col)]]}
    col[!is.na(col)] <- alpha2hex(col[!is.na(col)], alpha)
    ## if(!is.null(alpha) && !is.character(alpha)){alpha <- alpha2hex2(alpha)}
    ## col[!is.na(col)] <- paste0(col[!is.na(col)], alpha)
    matlines(x=xl, y=yl, lty=1, lwd=lwd, col=col, ...)
}

thist <- function(x, n=NULL, type=8, pretty=FALSE, plot=FALSE, extendbreaks=FALSE){
    if(!is.list(x)){x <- list(x)}
    if(!is.list(n)){n <- list(n)}
    out <- list()
    for(i in 1:length(x)){
        ax <- x[[i]]
        an <- n[[(i-1)%%length(n)+1]]
    ax <- ax[!is.na(ax) & is.finite(ax)]
    if(is.null(an)){an <- (round(sqrt(length(ax))/2))}
    if(length(an)==1 && (is.na(an) || an=='i' || an=='integer')){breaks <- (round(min(ax))-0.5):(round(max(ax))+0.5)}
    else if(length(an) > 1 || is.character(an)){breaks <- an}
    else if(length(an) == 1 && an > 0){
        rg <- range(ax)
        if(diff(rg)==0){rg <- rg + c(-0.5,0.5)}
        breaks <- seq(rg[1], rg[2], length.out=an+1)}
    else if(length(an) == 1 && an < 0){
        rg <- range(ax)
        if(diff(rg)==0){rg <- rg + c(-0.5,0.5)}
        breaks <- seq(rg[1], rg[2]-an, by=-an)
        breaks <- breaks - (breaks[length(breaks)]-rg[2])/2
    }
    else {print('Error with n')}
        if(!is.null(pretty) && pretty){
            breaks <- pretty(ax, n=length(breaks)-1)
        }
        if(extendbreaks){
         breaks <- c(-Inf,breaks,+Inf)
        }
        out <- c(out,list(hist(x=ax, breaks=breaks, plot=FALSE)))
    }
    if(plot){
        tplot(x=lapply(out,function(xx)xx$breaks),
              y=lapply(out,function(xx)xx$density),ylim=c(0,NA))
    }else{
        if(length(out)==1){unlist(out,recursive=F)}else{out}
    }
}

tquant <- function(x, probs=c(1:3)/4, na.rm=TRUE, names=TRUE, type=6, ...){
    quantile(x=x, probs=probs, na.rm=na.rm, names=names, type=type, ...)
}

tmad <- function(x){mad(x, constant=1, na.rm=TRUE)}

tsummary <- function(x){
    x <- cbind(x)
    apply(x, 2, function(xx){
        c(tquant(xx, c(2.5/100,1/8,2/8,4/8,6/8,7/8,97.5/100)), MAD=mad(xx,constant=1,na.rm=T), IQR=IQR(xx,na.rm=T), mean=mean(xx,na.rm=T), sd=sd(xx,na.rm=T), min=min(xx,na.rm=T), max=max(xx,na.rm=T), NAs=sum(is.na(xx)))
    })
}




## > trythese <- foreach(i=1:1e6, .combine=cbind)%dorng%{tests <- rt(n=4096,df=4); cbind(1-tquant(tests,qlev)/qt(qlev,df=4))*100}
## > tsummary(t(trythese))
##             0.5%      99.5%       2.5%      97.5%       12.5%       87.5%         25%
## 12.5%  -8.367148  -8.360343  -4.183885  -4.187169  -3.0651119  -3.0556805  -3.9053598
## 25%    -4.968965  -4.963558  -2.469266  -2.474195  -1.8028822  -1.7991997  -2.2961003
## 50%    -0.480020  -0.473392  -0.116171  -0.111059  -0.0326341  -0.0324019  -0.0218354
## 75%     3.652310   3.660276   2.165386   2.169362   1.7176493   1.7144059   2.2437016
## 87.5%   6.379180   6.385466   3.724190   3.730104   2.9415736   2.9345811   3.8349984
## MAD     6.373633   6.374755   3.433849   3.439678   2.6096926   2.6042186   3.3651678
## IQR     8.621234   8.623820   4.634616   4.643554   3.5205249   3.5135955   4.5397950
## mean   -0.874277  -0.864254  -0.199958  -0.199712  -0.0542001  -0.0541560  -0.0351858
## sd      6.475902   6.481717   3.443751   3.447723   2.6056377   2.6032365   3.3624852
## min   -45.108882 -43.715590 -18.702222 -19.914109 -13.5589661 -12.7060161 -17.5715736
## max    22.906166  23.319217  14.659703  15.753545  12.0475251  11.3247566  16.5761713
## NAs     0.000000   0.000000   0.000000   0.000000   0.0000000   0.0000000   0.0000000
##               75%  50%
## 12.5%  -3.9120108 -Inf
## 25%    -2.2984156 -Inf
## 50%    -0.0249747 -Inf
## 75%     2.2335347  Inf
## 87.5%   3.8200256  Inf
## MAD     3.3594907   NA
## IQR     4.5319444  Inf
## mean   -0.0388496  NaN
## sd      3.3564239  NaN
## min   -17.2698215 -Inf
## max    16.4433934  Inf
## NAs     0.0000000    0
