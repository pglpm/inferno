library('nimble')

## Dirichlet-process mixture of product kernels
## Two kinds of variates, two variates per kind:
## real (kernel: gaussian, prior: normal + compound-gamma)
## categorical (kernel: categorical distr, prior: Dirichlet)

nvars <- 10 # switch to include/exclude categorical components
ndata <- 500
ncens <- 8
nclusters <- 64
##
set.seed(229)
rdata <- matrix(rnorm(ndata*nvars, mean=0, sd=2), nrow=ndata)
## 1: 2.9, 2: 5
bounds <- matrix(seq(2.9,5,length.out=ncens), nrow=ndata,ncol=ncens,byrow=T)
xbounds <- matrix(+Inf,nrow=ndata,ncol=nvars-ncens)
## censored <- cbind(rep(1,ndata), rep(0,ndata))
censored <- rdata[,1:ncens] >= bounds
allcensored <- rdata >= cbind(bounds,xbounds)
##
rdata2 <- rdata
rdata[allcensored] <- NA

## cbind(rdata2,rdata,censored)
## bounds[!is.na(rdata)] <- signs[!is.na(rdata)]*Inf


datapoints <- c(
    list(Rdata = rdata, Censored=censored*1L)
)

constants <- list(ndata=ndata, nclusters=nclusters, bounds=bounds, nvars=nvars, ncens=ncens, ncens1=ncens+1)

initsFunction <- function(){
    c(
        list(
            palphas = c(1/3,1/3,1/3),
            qalpha0 = rep(0.25/nclusters, nclusters), # cluster probabilities
            alpha0 = 2,
            q=rep(1/nclusters,nclusters),
            Labels = rep(1, ndata), # all data to first cluster
            Rrates = rinvgamma(nvars, shape=0.5, scale=1),
            Rmeans = matrix(rnorm(nvars*nclusters, mean=0, sd=3), nrow=nvars, ncol=nclusters),
            Rvariances = matrix(rinvgamma(nvars*nclusters, shape=0.5, rate=1), nrow=nvars, ncol=nclusters)
        )
    )
}

infmixture <- nimbleCode({
    alpha0 ~ dcat(prob=palphas[1:3])
    qalpha[1:nclusters] <- 2^alpha0 * qalpha0[1:nclusters]
    q[1:nclusters] ~ ddirch(alpha=qalpha[1:nclusters])
    ##
    for(variate in 1:nvars){
        Rrates[variate] ~ dinvgamma(shape=0.5, scale=1)
    }
    for(cluster in 1:nclusters){
        for(variate in 1:nvars){# real variates
            Rmeans[variate, cluster] ~ dnorm(mean=0, var=3*3)
            Rvariances[variate, cluster] ~ dinvgamma(shape=0.5, rate=Rrates[variate])
        }
    }
    ##
    for(datum in 1:ndata){
        Labels[datum] ~ dcat(prob=q[1:nclusters])
        ##
        for(variate in 1:ncens){# real variates
            Censored[datum, variate] ~ dinterval(Rdata[datum, variate], bounds[datum, variate])
        }
        for(variate in 1:nvars){
            Rdata[datum, variate] ~ dnorm(mean=Rmeans[variate, Labels[datum]], var=Rvariances[variate, Labels[datum]])
        }
    }
})

model <- nimbleModel(code=infmixture, name='model',
                     constants=constants,
                     inits=initsFunction(),
                     data=datapoints,
                     dimensions=c(
                         list(Censored=c(ndata,ncens),
                              Labels=ndata,
                              q=nclusters,
                              Rrates=nvars,
                              Rmeans=c(nvars, nclusters),
                              Rvariances=c(nvars, nclusters))
                     )
                     )

Cmodel <- compileNimble(model, showCompilerOutput=FALSE)


confmodel <- configureMCMC(Cmodel,
                           monitors=c('q','alpha0',
                                      ## 'Labels',
                                      ## 'Rdata[2, 1]',
                                      ## 'Rdata[4, 2]',
                                      'Rrates', 'Rmeans', 'Rvariances'
                           ))

mcsampler <- buildMCMC(confmodel)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

##
timestart <- Sys.time()
Cmodel$setInits(initsFunction())
output <- Cmcsampler$run(niter=1100, thin=1, reset=T, resetMV=TRUE, nburnin=100, time=T)
timestop <- Sys.time() - timestart
timestop
## ncens=2:
## Time difference of 1.28007 mins
## Rvariances     Rmeans     Labels          q     Rrates      Rdata     alpha0 
##  31.709330  27.369537  16.857415   0.246917   0.236185   0.048782   0.016083 
##
## ncens=8:
## Time difference of 1.20523 mins
## Rvariances     Rmeans     Labels          q     Rrates      Rdata     alpha0 
##  29.695191  25.041273  16.612910   0.237851   0.223879   0.163691   0.014620 

dpsamples <- as.matrix(Cmcsampler$mvSamples)
times <- Cmcsampler$getTimes()
names(times) <- sapply(confmodel$getSamplers(),function(x)x$target)

## sum(times[grepl('^Rdata',names(times))])
##
prefs <- unique(sub('^([^[]+)(\\[.*\\])', '\\1', names(times)))
sort(sapply(prefs, function(x)sum(times[grepl(x,names(times))])),decreasing=T)



## Cmodel$plotGraph()
