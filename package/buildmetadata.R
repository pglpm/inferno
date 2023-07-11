buildmetadata <- function(data, file=NULL){
    gcd2 <- function(a, b){suppressWarnings( if (b == 0) a else Recall(b, a %% b) )}
    gcd <- function(...){suppressWarnings(Reduce(gcd2, c(...)))}
    ##
    datafile <- NULL
    if(is.character(data) && file.exists(data)){
        datafile <- data
        data <- fread(datafile, na.strings='')
    }
    data <- as.data.table(data)
    metadata <- data.table()
    for(xn in names(data)){
        ## print(xn)
        x <- data[[xn]]
        x <- x[!is.na(x)]
        if(is.numeric(x)){
            ## print('numeric')
            Q1 <- loval <- quantile(x, probs=0.25, type=6)
            Q2 <- meval <- quantile(x, probs=0.5, type=6)
            Q3 <- hival <- quantile(x, probs=0.75, type=6)
            dmin <- min(x)
            dmax <- max(x)
            if(loval == hival){
                loval <- if(sum(x < Q1)>0){max(x[x < Q1])}else{max(x[x <= Q1])}
                hival <- if(sum(x > Q3)>0){min(x[x > Q3])}else{min(x[x >= Q3])}
            }
        }else{
            loval <- meval <- hival <- dmin <- dmax <- NA
        }
        if(length(unique(x)) == 2){# seems binary variate
            vtype <- 'binary'
            vn <- 2
            vd <- NA
            domainmin <- NA
            domainmax <- NA
            censormin <- TRUE
            censormax <- TRUE
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
            censormin <- TRUE
            censormax <- TRUE
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
                domainmin <- signif(min(x) - 3*diff(range(x)),1)
                domainmax <- signif(max(x) + 3*diff(range(x)),1)
                censormin <- FALSE
                censormax <- FALSE
                plotmin <- min(x) - (Q3-Q1)/2
                plotmax <- max(x) + (Q3-Q1)/2
                ##
                if(all(x >= 0)){ # seems to be strictly positive
                    domainmin <- 0
                    plotmin <- max((domainmin+min(x))/2, plotmin)
                }
                ix <- x[!(x %in% range(x))] # exclude boundary values
                repindex <- mean(table(ix)) # average of repeated inner values
                ## contindex <- length(unique(diff(sort(unique(ix)))))/length(ix) # check for repeated values
                if(sum(x == min(x)) > repindex){ # seems to be left-singular
                    domainmin <- plotmin <- min(x)
                    censormin <- TRUE
                }
                if(sum(x == max(x)) > repindex){ # seems to be right-singular
                    domainmax <- plotmax <- max(x)
                    censormax <- TRUE
                }
            }else{# ordinal
                vtype <- 'ordinal'
                if(dd >= 1){ # seems originally integer
                    domainmin <- min(1, x)
                    domainmax <- max(x)
                    vn <- domainmax - domainmin + 1
                    vd <- NA
                    censormin <- TRUE
                    censormax <- TRUE
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
                    domainmin <- signif(min(x) - 4*diff(range(x)),1)
                    domainmax <- signif(max(x) + 4*diff(range(x)),1)
                    censormin <- FALSE
                    censormax <- FALSE
                    ## location <- Q2
                    ## scale <- (Q3-Q1)/2
                    plotmin <- min(x) - (Q3-Q1)/2
                    plotmax <- max(x) + (Q3-Q1)/2
                    ix <- x[!(x %in% range(x))] # exclude boundary values
                    repindex <- mean(table(ix)) # average of repeated inner values
                    ## contindex <- length(unique(diff(sort(unique(ix)))))/length(ix) # check for repeated values
                    if(all(x >= 0)){ # seems to be strictly positive
                        domainmin <- 0
                        plotmin <- max((domainmin+min(x))/2, plotmin)
                    }
                    ix <- x[!(x %in% range(x))] # exclude boundary values
                    repindex <- mean(table(ix)) # average of repeated inner values
                    ## contindex <- length(unique(diff(sort(unique(ix)))))/length(ix) # check for repeated values
                    if(sum(x == min(x)) > repindex){ # seems to be left-singular
                        domainmin <- plotmin <- min(x)
                        censormin <- TRUE
                    }
                    if(sum(x == max(x)) > repindex){ # seems to be right-singular
                        domainmax <- plotmax <- max(x)
                        censormax <- TRUE
                    }
                }# end rounded
            }# end integer
            vval <- NULL
        }# end numeric
        ##
        metadata <- rbind(metadata,
                          c(list(name=xn, type=vtype, datamin=dmin, datamax=dmax, datamaxrep=max(table(x)), dataNvalues=length(unique(x)), Nvalues=vn, rounding=vd, domainmin=domainmin, domainmax=domainmax, minincluded=censormin, maxincluded=censormax, centralvalue=meval, lowvalue=loval, highvalue=hival, plotmin=plotmin, plotmax=plotmax),
                            as.list(vval)
                            ), fill=TRUE)
    }
    ## metadata <- cbind(name=names(data), metadata)
    
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
