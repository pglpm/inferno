library('data.table')
library('nimble')

varnames <- c('Apoe4_', 'Subgroup_num_',  'Gender_num_', 'GDTOTAL_gds',  'AVDEL30MIN_neuro', 'AVDELTOT_neuro',  'CATANIMSC_neuro', 'ANARTERR_neuro',  'RAVLT_immediate', 'TRAASCOR_neuro',  'TRABSCOR_neuro', 'AGE',  'LRHHC_n_long')
##
dt0 <- fread('~/repositories/ADBayes/scripts/bl.csv')
all(varnames %in% colnames(dt0))
dt <- dt0[,..varnames]
##
nas <- (1:nrow(dt))[apply(dt,1,function(x)any(is.na(x)))]
gds6 <- (1:nrow(dt))[dt[['GDTOTAL_gds']] >= 6]
## ##
## dt3 <- fread('~/repositories/ADBayes/scripts/3.13_pycharet_K50.csv')
## diffs <- sapply(varnames,function(x){
##     x1 <- sort(dt[[x]][-nas])
##     x2 <- sort(dt3[[x]])
##     abs(x1-x2)/(x1+x2+(x1+x2==0))/2
## })
## range(diffs)
## ## [1] 0.00000e+00 8.94359e-15
## ##
## dtred <- dt[-nas]
## ##
## fwrite(dt, 'ingrid_data_all.csv')
## fwrite(dt[-gds6], 'ingrid_data_nogds6.csv')
## fwrite(dt[-nas], 'ingrid_data_nonan.csv')
## fwrite(dt[-c(nas,gds6)], 'ingrid_data_nonangds6.csv')
## rm(dt0,dt3)

varinfo <- cbind(
    min=c(0,0,0,0,0,0,0,0,0,0,0,0,0),
    max=c(1,1,1,6,15,15,63,50,75,150,300,Inf,Inf),
    n=c(2,2,2,7,16,16,64,51,76,Inf,Inf,Inf,Inf)
)
rownames(varinfo) <- varnames

subsel <- 4:nrow(varinfo)
nclusters <- 64
##
set.seed(229)
idata <- cbind(
    dt[[4]],
    dt[[5]],
    dt[[6]],
    dt[[7]],
    dt[[8]],
    dt[[9]]
)
rdata <- cbind(
    log(dt[[12]]),
    log(dt[[13]])
)
nivars <- ncol(idata)
nrvars <- ncol(rdata)
ndata <- nrow(rdata)
nint <- ncol(idata)
iint <- 4:9
##
initint <-  t((t(idata)-varinfo[iint,'min'])/(varinfo[iint,'max']-varinfo[iint,'min']))*2-1
initint[is.na(initint)] <- 0

##
iseqs <- apply(varinfo[iint,], 1, function(x) seq(-1,1,length.out=x['n']))
breaks <- sapply(iseqs, function(x){(x[-1]+x[-length(x)])/2})
nbreaks <- sapply(breaks,length)
maxnbreaks <- max(nbreaks)
mbreaks <- t(sapply(breaks,function(x){c(x,rep(+Inf,maxnbreaks-length(x)))}))
##
## censored <- foreach(datum=rdata[,1:ncens,drop=F], breaks=bounds, .combine=cbind)%do%{
##     sapply(datum, function(x){
##         which(x <= c(breaks,+Inf) & x > c(-Inf,breaks))-1
##     })
## }
## rdata2 <- rdata
## rdata[,1:ncens] <- NA
##

datapoints <- c(
    list(Rdata = rdata, Idata = idata)
)
##
constants <- list(ndata=ndata, nclusters=nclusters, nrvars=nrvars, nivars=nivars, mbreaks=mbreaks, maxnbreaks=maxnbreaks)
##
initsFunction <- function(){
    c(
        list(
            palphas = c(1/3,1/3,1/3),
            qalpha0 = rep(0.25/nclusters, nclusters), # cluster probabilities
            alpha0 = 3,
            q=rep(1/nclusters,nclusters),
            Labels = rep(1, ndata), # all data to first cluster
            Rrates = rinvgamma(nrvars, shape=0.5, scale=1),
            Irates = rinvgamma(nivars, shape=0.5, scale=1),
            Rmeans = matrix(rnorm(nrvars*nclusters, mean=0, sd=3), nrow=nrvars, ncol=nclusters),
            Imeans = matrix(rnorm(nivars*nclusters, mean=0, sd=3), nrow=nivars, ncol=nclusters),
            Rvariances = matrix(rinvgamma(nrvars*nclusters, shape=0.5, rate=1), nrow=nrvars, ncol=nclusters),
            Ivariances = matrix(rinvgamma(nivars*nclusters, shape=0.5, rate=1), nrow=nivars, ncol=nclusters),
            Ihidden=initint
        )
    )
}

infmixture <- nimbleCode({
    alpha0 ~ dcat(prob=palphas[1:3])
    qalpha[1:nclusters] <- 2^alpha0 * qalpha0[1:nclusters]
    q[1:nclusters] ~ ddirch(alpha=qalpha[1:nclusters])
    ##
    for(variate in 1:nrvars){
        Rrates[variate] ~ dinvgamma(shape=0.5, scale=1)
    }
    for(variate in 1:nivars){
        Irates[variate] ~ dinvgamma(shape=0.5, scale=1)
    }
    for(cluster in 1:nclusters){
        for(variate in 1:nrvars){# real variates
            Rmeans[variate, cluster] ~ dnorm(mean=0, var=3*3)
            Rvariances[variate, cluster] ~ dinvgamma(shape=0.5, rate=Rrates[variate])
        }
        for(variate in 1:nivars){# real variates
            Imeans[variate, cluster] ~ dnorm(mean=0, var=3*3)
            Ivariances[variate, cluster] ~ dinvgamma(shape=0.5, rate=Irates[variate])
        }
    }
    ##
    for(datum in 1:ndata){
        Labels[datum] ~ dcat(prob=q[1:nclusters])
        ##
        for(variate in 1:nrvars){
            Rdata[datum, variate] ~ dnorm(mean=Rmeans[variate, Labels[datum]], var=Rvariances[variate, Labels[datum]])
        }
        for(variate in 1:nivars){# real variates
            Idata[datum, variate] ~ dinterval(Ihidden[datum, variate], mbreaks[variate,1:maxnbreaks])
            Ihidden[datum, variate] ~ dnorm(mean=Imeans[variate, Labels[datum]], var=Ivariances[variate, Labels[datum]])
        }
    }
})

model <- nimbleModel(code=infmixture, name='model',
                     constants=constants,
                     inits=initsFunction(),
                     data=datapoints,
                     dimensions=c(
                         list(
                              Labels=ndata,
                              q=nclusters,
                              Rrates=nrvars,
                              Irates=nivars,
                              Rmeans=c(nrvars, nclusters),
                              Imeans=c(nivars, nclusters),
                              Rvariances=c(nrvars, nclusters),
                             Ivariances=c(nivars, nclusters),
                             Ihidden=c(ndata,nivars)
                         )
                     )
                     )

Cmodel <- compileNimble(model, showCompilerOutput=FALSE)


confmodel <- configureMCMC(Cmodel,
                           monitors=c('q','alpha0',
                                      ## 'Labels',
                                      ## 'Rdata[2, 1]',
                                      ## 'Rdata[4, 2]',
                                      'Rrates', 'Rmeans', 'Rvariances',
                                      'Irates', 'Imeans', 'Ivariances'
                                      ),
                           monitors2='Idata')

mcsampler <- buildMCMC(confmodel)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

##
timestart <- Sys.time()
Cmodel$setInits(initsFunction())
output <- Cmcsampler$run(niter=100, thin=1, reset=TRUE, resetMV=TRUE, nburnin=1, time=T)

output <- Cmcsampler$run(niter=1024, thin=1, reset=FALSE, resetMV=TRUE, nburnin=0, time=T)
timestop <- Sys.time() - timestart
timestop
##
dpsamples <- as.matrix(Cmcsampler$mvSamples)
dpsamples2 <- as.matrix(Cmcsampler$mvSamples2)
times <- Cmcsampler$getTimes()
names(times) <- sapply(confmodel$getSamplers(),function(x)x$target)
## sum(times[grepl('^Rdata',names(times))])
##
prefs <- unique(sub('^([^[]+)(\\[.*\\])', '\\1', names(times)))
sort(sapply(prefs, function(x)sum(times[grepl(x,names(times))])),decreasing=T)
## with hidden:
## Time difference of 1.51792 mins
## Ivariances     Labels     Imeans Rvariances     Rmeans    Ihidden          q 
##  26.372392  21.895887  21.326742   8.408666   7.174604   4.418221   0.340741 
##     Irates     Rrates     alpha0      Rdata 
##   0.143012   0.053464   0.015036   0.009979 



## Cmodel$plotGraph()
