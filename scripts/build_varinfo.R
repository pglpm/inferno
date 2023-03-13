library('data.table')
library('png')
library('foreach')
library('doFuture')
library('doRNG')
registerDoFuture()
cat('\navailableCores: ')
cat(availableCores())
cat('\navailableCores-multicore: ')
cat(availableCores('multicore'))
if(Sys.info()['nodename']=='luca-HP-Z2-G9'){
    ncores <- 20}else{
    ncores <- 6}
cat(paste0('\nusing ',ncores,' cores\n'))
if(ncores>1){
    if(.Platform$OS.type=='unix'){
        plan(multicore, workers=ncores)
    }else{
        plan(multisession, workers=ncores)
    }
}else{
    plan(sequential)
}

dt <- fread('ingrid_data_nogds6.csv')
dt2 <- fread('testdata4.csv')
data(iris)
iris <- as.data.table(iris)
iris2 <- iris
iris2$Species <- as.integer(iris2$Species)


buildvarinfo <- function(data, file=NULL){
    gcd2 <- function(a, b){ if (b == 0) a else Recall(b, a %% b) }
    gcd <- function(...) Reduce(gcd2, c(...))
    ##
    if(is.character(data) && file.exists(data)){data <- fread(data)}
    data <- as.data.table(data)
    varinfo <- data.table()
    for(x in data){
        x <- x[!is.na(x)]
        transf <- 'identity' # temporary
        if(is.numeric(x)){
            Q1 <- quantile(x, probs=0.25, type=6)
            Q2 <- quantile(x, probs=0.5, type=6)
            Q3 <- quantile(x, probs=0.75, type=6)
            location <- Q2
            scale <- (Q3-Q1)/2
            if(scale == 0){scale <- diff(range(x))/2}
        }else{
            location <- NA
            scale <- NA
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
        }else{# discrete, continuous, or boundary-singular variate
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
                tmin <- NA
                tmax <- NA
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
                    tmin <- NA
                    tmax <- NA
                    ## location <- Q2
                    ## scale <- (Q3-Q1)/2
                    plotmin <- min(x) - (Q3-Q1)/2
                    plotmax <- max(x) + (Q3-Q1)/2
                    if(all(x > 0)){ # seems to be strictly positive
                        transf <- 'log'
                        vmin <- 0
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
                         c(list(type=vtype, Nvalues=vn, rounding=vd, domainmin=vmin, domainmax=vmax, censormin=tmin, censormax=tmax, location=location, scale=scale, plotmin=plotmin, plotmax=plotmax),
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

testf <- function(a){
    if(a==1){
        stop('abort')
    }else{
        print('continuing')
    }
    34
}


buildvarinfoaux <- function(data, varinfo, file=TRUE){
    if(is.character(data) && file.exists(data)){data <- fread(data)}
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
            vtype <- 'C'
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
            vtype <- 'L'
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
            if(!is.finite(xinfo$censormin) && !is.finite(xinfo$censormax) && (!is.finite(xinfo$rounding) || xinfo$rounding == 0)){ # no need for latent variate
                vtype <- 'R'
            }else{ # need latent variate
                vtype <- 'L'
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


## gcd <- function(vect){Reduce(function(x,y) ifelse(y, Recall(y, x %% y), x), as.list(vect))}
## gcdm <- function(...){Reduce(function(x,y) ifelse(y, Recall(y, x %% y), x), list(...))}


gcd2 <- function(a, b) {
  if (b == 0) a else Recall(b, a %% b)
}
gcd <- function(...) Reduce(gcd2, c(...))

dtx <- dt
dtx$extra <- rnorm(nrow(dtx))
dtx <- dtx[sample(1:nrow(dtx), min(10,nrow(dtx)))]
t(sapply(dtx, function(xx){
    xx <- xx[!is.na(xx)]
    testd <- unique(signif(diff(sort(unique(xx))),6))
    multi <- 10^(-min(floor(log10(testd))))
    round(gcd(testd*multi))/multi
}))





dtx <- iris
dtx$extra <- rnorm(nrow(dtx))
dtx <- dtx[sample(1:nrow(dtx), min(150,nrow(dtx)))]
t(sapply(dtx, function(xx){
    xx <- xx[!is.na(xx)]
    if(length(unique(xx)) > 2){ix <- xx[!(xx %in% range(xx))]}else{ix <- xx}
    dd0 <- diff(sort(xx))
    dd <- diff(sort(unique(xx)))
    q1 <- tquant(ix,0.25)
    q3 <- tquant(ix,0.75)
    c(
      left=sum(xx==min(xx)),
      right=sum(xx==max(xx)),
      diffratio=length(unique(dd))/length(dd),
      diffratio0=length(unique(dd0))/length(dd0),
      diffindex=length(unique(dd))/length(xx),
      diffindex0=length(unique(dd0))/length(xx),
      meanrep=mean(table(ix)),
      meanrepiqr=mean(table(ix[ix >= q1 & ix <=q3])),
      iqrrange=IQR(ix)/diff(range(xx)),
      min=min(xx),
      max=max(xx),
      int=is.integer(xx),
      unique=length(unique(xx)),
      uniqueratio=length(unique(xx))/length(xx),
      rg=diff(range(xx)),
      NULL
      )
}))
rm(dtx)

dtx <- dt
dtx$extra <- rnorm(nrow(dtx))
summary(t(sapply(1:100, function(xxx){set.seed(xxx)
    dtx <- dtx[sample(1:nrow(dtx), 10)]
    test <- t(sapply(dtx, function(xx){
        xx <- xx[!is.na(xx)]
        if(length(unique(xx)) > 2){ix <- xx[!(xx %in% range(xx))]}else{ix <- xx}
        dd <- diff(sort(unique(xx)))
        c(
            sum(xx==min(xx)),
            sum(xx==max(xx)),
            length(unique(dd))/length(dd),
            length(unique(dd))/length(xx),
            mean(table(ix))
        )
    }))
    c(sum(test['extra',5] < test[-c(13:14),5]),
      sum(test['extra',4] > test[-c(13:14),4]),
      sum(test['extra',3] > test[-c(13:14),3])
      )
})))

dtx <- dt
dtx$extra <- rnorm(nrow(dtx))
resu <- (t(sapply(1:100, function(xxx){set.seed(xxx)
    dtx <- dtx[sample(1:nrow(dtx), 150)]
    test <- t(sapply(dtx, function(xx){
        xx <- xx[!is.na(xx)]
        if(length(unique(xx)) > 2){ix <- xx[!(xx %in% range(xx))]}else{ix <- xx}
        dd <- diff(sort(unique(xx)))
        c(
            sum(xx==min(xx)),
            sum(xx==max(xx)),
            length(unique(dd))/length(dd),
            length(unique(dd))/length(xx),
            mean(table(ix))
        )
    }))
    c(
      min(test[13:14,4]),
      max(test[-c(13:14),4])
      )
})))
summary(resu)
