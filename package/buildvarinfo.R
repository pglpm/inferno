buildvarinfo <- function(data, file=NULL){
    gcd2 <- function(a, b){ if (b == 0) a else Recall(b, a %% b) }
    gcd <- function(...) Reduce(gcd2, c(...))
    ##
    if(is.character(data) && file.exists(data)){data <- fread(data, na.strings='')}
    data <- as.data.table(data)
    varinfo <- data.table()
    for(x in data){
        x <- x[!is.na(x)]
        transf <- 'identity' # temporary
        if(is.numeric(x)){
            Q1 <- loval <- quantile(x, probs=0.25, type=6)
            Q2 <- meval <- quantile(x, probs=0.5, type=6)
            Q3 <- hival <- quantile(x, probs=0.75, type=6)
            if(loval == hival){
                loval <- (Q1 + min(x))/2
                hival <- (Q1 + max(x))/2
            }
        }else{
            loval <- meval <- hival <- NA
        }
        if(length(unique(x)) == 2){# seems binary variate
            vtype <- 'binary'
            vn <- 2
            vd <- NA
            vmin <- NA
            vmax <- NA
            tmin <- NA
            tmax <- NA
            vval <- sort(as.character(unique(x)))
            names(vval) <- paste0('V',1:2)
            plotmin <- NA
            plotmax <- NA
        }else if(!is.numeric(x)){# nominal variate
            vtype <- 'nominal'
            vn <- length(unique(x))
            vd <- NA
            vmin <- NA # Nimble index categorical from 1
            vmax <- NA
            tmin <- NA
            tmax <- NA
            vval <- sort(as.character(unique(x)))
            names(vval) <- paste0('V',1:vn)
            plotmin <- NA
            plotmax <- NA
        }else{# ordinal, continuous, censored, or discretized variate
            ud <- unique(signif(diff(sort(unique(x))),3)) # differences
            rx <- diff(range(x))
            multi <- 10^(-min(floor(log10(ud))))
            dd <- round(gcd(ud*multi))/multi # greatest common difference
            ##
            if(dd/rx < 1e-3){ # consider it as continuous
                ## temporary values
                vtype <- 'continuous'
                vn <- Inf
                vd <- 0
                vmin <- -Inf
                vmax <- +Inf
                tmin <- -Inf
                tmax <- +Inf
                plotmin <- min(x) - (Q3-Q1)/2
                plotmax <- max(x) + (Q3-Q1)/2
                ##
                ix <- x[!(x %in% range(x))] # exclude boundary values
                repindex <- mean(table(ix)) # average of repeated inner values
                ## contindex <- length(unique(diff(sort(unique(ix)))))/length(ix) # check for repeated values
                if(sum(x == min(x)) > repindex){ # seems to be left-singular
                    tmin <- min(x)
                    plotmin <- tmin
                }
                if(sum(x == max(x)) > repindex){ # seems to be right-singular
                    tmax <- max(x)
                    plotmax <- tmax
                }
                if(all(x > 0)){ # seems to be strictly positive
                    transf <- 'log'
                    vmin <- 0
                    tmin <- max(vmin, tmin)
                    ## location <- log(Q2)
                    ## scale <- (log(Q3) - log(Q1))/2
                    plotmin <- max(vmin, plotmin)
                }
            }else{# ordinal
                vtype <- 'ordinal'
                if(dd >= 1){ # seems originally integer
                    transf <- 'Q'
                    vmin <- min(1, x)
                    vmax <- max(x)
                    vn <- vmax - vmin + 1
                    vd <- 1
                    tmin <- NA
                    tmax <- NA
                    ## location <- NA # (vn*vmin-vmax)/(vn-1)
                    ## scale <- NA # (vmax-vmin)/(vn-1)
                    plotmin <- max(vmin, min(x) - IQR(x)/2)
                    plotmax <- max(x)
                }else{ # seems a rounded continuous variate
                    vtype <- 'continuous'
                    vn <- Inf
                    vd <- dd
                    vmin <- -Inf
                    vmax <- +Inf
                    tmin <- -Inf
                    tmax <- +Inf
                    ## location <- Q2
                    ## scale <- (Q3-Q1)/2
                    plotmin <- min(x) - (Q3-Q1)/2
                    plotmax <- max(x) + (Q3-Q1)/2
                    if(all(x > 0)){ # seems to be strictly positive
                        transf <- 'log'
                        vmin <- 0
                        tmin <- max(vmin, tmin)
                        ## location <- log(Q2)
                        ## scale <- (log(Q3) - log(Q1))/2
                        plotmin <- max(vmin, plotmin)
                    }
                }# end rounded
            }# end integer
            vval <- NULL
        }# end numeric
        ##
        varinfo <- rbind(varinfo,
                         c(list(type=vtype, Nvalues=vn, rounding=vd, domainmin=vmin, domainmax=vmax, censormin=tmin, censormax=tmax, centralvalue=meval, lowvalue=loval, highvalue=hival, plotmin=plotmin, plotmax=plotmax),
                           as.list(vval)
                         ), fill=TRUE)
    }
    varinfo <- cbind(name=names(data), varinfo)
    if(!is.null(file)){
        file <- paste0(sub('.csv$', '', file), '.csv')
        fwrite(varinfo, file)
        cat(paste0('Saved proposal variate-info file to ', file, '\n'))
    }else{
        varinfo
    }
}
