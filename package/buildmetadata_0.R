buildmetadata <- function(data, file=NULL){
    require('data.table')
    gcd2 <- function(a, b){ if (b == 0) a else Recall(b, a %% b) }
    gcd <- function(...) Reduce(gcd2, c(...))
    ##
    datafile <- NULL
    if(is.character(data) && file.exists(data)){
        datafile <- data
        data <- fread(datafile, na.strings='')
    }
    data <- as.data.table(data)
    metadata <- data.table()
    for(xn in names(data)){
        x <- data[[xn]]
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
            domainmin <- NA
            domainmax <- NA
            censormin <- NA
            censormax <- NA
            vval <- sort(as.character(unique(x)))
            names(vval) <- paste0('V',1:2)
            loval <- meval <- hival <- NA
            plotmin <- NA
            plotmax <- NA
        }else if(!is.numeric(x)){# nominal variate
            vtype <- 'nominal'
            vn <- length(unique(x))
            vd <- NA
            domainmin <- NA # Nimble index categorical from 1
            domainmax <- NA
            censormin <- NA
            censormax <- NA
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
                domainmin <- -Inf
                domainmax <- +Inf
                censormin <- NA
                censormax <- NA
                plotmin <- min(x) - (Q3-Q1)/2
                plotmax <- max(x) + (Q3-Q1)/2
                ##
                ix <- x[!(x %in% range(x))] # exclude boundary values
                repindex <- mean(table(ix)) # average of repeated inner values
                ## contindex <- length(unique(diff(sort(unique(ix)))))/length(ix) # check for repeated values
                if(sum(x == min(x)) > repindex){ # seems to be left-singular
                    censormin <- min(x)
                    plotmin <- censormin
                }
                if(sum(x == max(x)) > repindex){ # seems to be right-singular
                    censormax <- max(x)
                    plotmax <- censormax
                }
                if(all(x > 0)){ # seems to be strictly positive
                    transf <- 'log'
                    domainmin <- 0
                    censormin <- max(domainmin, censormin)
                    ## location <- log(Q2)
                    ## scale <- (log(Q3) - log(Q1))/2
                    plotmin <- max((domainmin+min(x))/2, plotmin)
                }
            }else{# ordinal
                vtype <- 'ordinal'
                if(dd >= 1){ # seems originally integer
                    transf <- 'Q'
                    domainmin <- min(1, x)
                    domainmax <- max(x)
                    vn <- domainmax - domainmin + 1
                    vd <- NA
                    censormin <- NA
                    censormax <- NA
                    ## location <- NA # (vn*domainmin-domainmax)/(vn-1)
                    ## scale <- NA # (domainmax-domainmin)/(vn-1)
                    plotmin <- max(domainmin, min(x) - (Q3-Q1)/2)
                    plotmax <- max(x)
                }else{ # seems a rounded continuous variate
                    vtype <- 'continuous'
                    vn <- Inf
                    vd <- dd
                    if(diff(range(x))/vd > 256){
                    message('\nNOTE: variate ',xn,' is reported as "rounded",\nbut consider the possibility of treating it as continuous,\nby setting its "rounding" to 0 in the metadata file.\n')
                    }
                    domainmin <- -Inf
                    domainmax <- +Inf
                    censormin <- NA
                    censormax <- NA
                    ## location <- Q2
                    ## scale <- (Q3-Q1)/2
                    plotmin <- min(x) - (Q3-Q1)/2
                    plotmax <- max(x) + (Q3-Q1)/2
                    ix <- x[!(x %in% range(x))] # exclude boundary values
                    repindex <- mean(table(ix)) # average of repeated inner values
                    ## contindex <- length(unique(diff(sort(unique(ix)))))/length(ix) # check for repeated values
                    if(sum(x == min(x)) > repindex){ # seems to be left-singular
                        censormin <- min(x)
                        plotmin <- max(plotmin, censormin)
                    }
                    if(sum(x == max(x)) > repindex){ # seems to be right-singular
                        censormax <- max(x)
                        plotmax <- min(plotmax, censormax)
                    }
                    if(all(x > 0)){ # seems to be strictly positive
                        transf <- 'log'
                        domainmin <- 0
                        censormin <- max(domainmin, censormin)
                        ## location <- log(Q2)
                        ## scale <- (log(Q3) - log(Q1))/2
                        plotmin <- max(domainmin, plotmin)
                    }
                }# end rounded
            }# end integer
            vval <- NULL
        }# end numeric
        ##
        metadata <- rbind(metadata,
                          c(list(type=vtype, Nvalues=vn, rounding=vd, domainmin=domainmin, domainmax=domainmax, censormin=censormin, censormax=censormax, centralvalue=meval, lowvalue=loval, highvalue=hival, plotmin=plotmin, plotmax=plotmax),
                            as.list(vval)
                            ), fill=TRUE)
    }
    metadata <- cbind(name=names(data), metadata)
    
    if(!missing(file) && file!=FALSE){# must save to file
        if(is.character(file)){
            file <- paste0(sub('.csv$', '', file), '.csv')
        }else{
            file <- paste0('metadata_', datafile)
            file <- paste0(sub('.csv$', '', file), '.csv')
        }
        if(file.exists(file)){
            file.rename(from=file, to=paste0(sub('.csv$', '', file), '_bak',format(Sys.time(), '%y%m%dT%H%M%S'),'.csv'))
        }
        fwrite(metadata, file)
        cat(paste0('Saved proposal metadata file as ', file, '\n'))
    }else{
        metadata
    }
}
