buildvarinfoaux <- function(data, varinfo, file=TRUE){
    if(is.character(data) && file.exists(data)){data <- fread(data, na.strings='')}
    data <- as.data.table(data)
    if(is.character(varinfo) && file.exists(varinfo)){
        varinfoname <- varinfo
        varinfo <- fread(varinfo)
    }
    varinfo <- as.data.table(varinfo)
    ## consistency checks
    if(!identical(varinfo$name, colnames(data))){
        stop('ERROR: mismatch in variate names or order')
    }
    ##
    ##Q <- readRDS('Qfunction512.rds')
    ##
    varinfoaux <- data.table()
    for(xn in colnames(data)){
        x <- data[[xn]]
        x <- x[!is.na(x)]
        xinfo <- as.list(varinfo[name == xn])
        xinfo$type <- tolower(xinfo$type)
        ordinal <- NA
        cens <- NA
        rounded <- NA
        transf <- 'identity' # temporary
        vval <- xinfo[grep('^V[0-9]+$', names(xinfo))]
        ## print(xn)
        ## str(vval)
        Q1 <- NA
        Q2 <- NA
        Q3 <- NA
        if(xinfo$type == 'binary'){# seems binary variate
            if(length(unique(x)) != 2){
                cat('Warning: inconsistencies with variate ', xn, '\n')
            }
            vtype <- 'B'
            vn <- xinfo$Nvalues
            vd <- xinfo$rounding/2
            vmin <- 0
            vmax <- 1
            tmin <- -Inf
            tmax <- +Inf
            location <- 0
            scale <- 1
            plotmin <- 0
            plotmax <- 1
        }else if(xinfo$type == 'nominal'){# nominal variate
            vtype <- 'N'
            vn <- xinfo$Nvalues
            vd <- 0.5
            vmin <- 1 # Nimble index categorical from 1
            vmax <- vn
            tmin <- -Inf
            tmax <- +Inf
            location <- 0
            scale <- 1
            plotmin <- 1
            plotmax <- vn
        }else if(xinfo$type == 'ordinal'){
            vtype <- 'O'
            transf <- 'Q'
            ordinal <- TRUE
            vn <- xinfo$Nvalues
            vd <- 0.5
            vmin <- xinfo$domainmin
            vmax <- xinfo$domainmax
            tmin <- -Inf
            tmax <- +Inf
            ##vval <- as.vector(xinfo[paste0('V',1:vn)], mode='character')
            location <- (vn*vmin - vmax)/(vn - 1)
            scale <- (vmax - vmin)/(vn - 1)
            plotmin <- 1
            plotmax <- vn
        }else if(xinfo$type == 'continuous'){
            vn <- +Inf
            vd <- xinfo$rounding/2
            rounded <- (vd > 0)
            vmin <- xinfo$domainmin
            vmax <- xinfo$domainmax
            tmin <- max(xinfo$censormin, -Inf, na.rm=TRUE)
            tmax <- min(xinfo$censormax, +Inf, na.rm=TRUE)
            cens <- any(is.finite(c(tmin,tmax)))
            location <- xinfo$location
            scale <- xinfo$scale
            Q1 <- quantile(x, probs=0.25, type=6)
            Q2 <- quantile(x, probs=0.5, type=6)
            Q3 <- quantile(x, probs=0.75, type=6)
            plotmin <- xinfo$plotmin
            plotmax <- xinfo$plotmax
            if(is.finite(xinfo$domainmin) && is.finite(xinfo$domainmax)){ # needs transformation
                transf <- 'probit'
                location <- qnorm((location-vmin)/(vmax-vmin))
                scale <- scale/((vmax-vmin)*dnorm(location))
            }else if(is.finite(xinfo$domainmin)){
                transf <- 'log'
                scale <- scale/(location-vmin)
                location <- log(location-vmin)
            }else if(is.finite(xinfo$domainmax)){
                transf <- 'logminus'
                scale <- scale/(vmax-location)
                location <- log(vmax-location)
            }
            if(xinfo$rounding > 0){ # continuous discretized
                vtype <- 'D'
            }
            else if(is.finite(xinfo$censormin) || is.finite(xinfo$censormax)){ # censored
                vtype <- 'C'
            }else{ # continuous
                vtype <- 'R'
            }
        }else{
            stop(paste0('ERROR: unknown variate type for ', xn))
        }
        ##
        ## print(varinfoaux[nrow(varinfoaux)])
        ## print(                         as.data.table(c(list(name=xn, type=vtype, transform=transf, Nvalues=vn, step=vd, domainmin=vmin, domainmax=vmax, censormin=tmin, censormax=tmax, tlocation=location, tscale=scale, plotmin=plotmin, plotmax=plotmax, Q1=Q1, Q2=Q2, Q3=Q3),
        ##                    vval
        ##                    )))
        varinfoaux <- rbind(varinfoaux,
                         c(list(name=xn, mcmctype=vtype, censored=cens, rounded=rounded, transform=transf, Nvalues=vn, step=vd, domainmin=vmin, domainmax=vmax, censormin=tmin, censormax=tmax, tlocation=location, tscale=scale, plotmin=plotmin, plotmax=plotmax, Q1=Q1, Q2=Q2, Q3=Q3),
                           vval
                           ), fill=FALSE)
    }
    if(is.character(file) || (is.logical(file) && file)){ # must save to file
        if(is.character(file)){
            file <- paste0(sub('.rds$', '', file), '.rds')
        }else{
            if(file.exists('varinfoaux.rds')){
                file <- paste0('varinfoaux_',format(Sys.time(), '%y%m%dT%H%M%S'),'.rds')
            }else{
                file <- 'varinfoaux.rds'
            }
        }
        fwrite(varinfoaux, file)
        cat(paste0('Saved auxiliary variate-info file to ', file, '\n'))
    }else{
        varinfoaux
    }
}
