library('nimble')

## Dirichlet-process mixture of product kernels
## Two kinds of variates, two variates per kind:
## real (kernel: gaussian, prior: normal + compound-gamma)
## categorical (kernel: categorical distr, prior: Dirichlet)

nvars <- 10 # switch to include/exclude categorical components
ndata <- 500

nclusters <- 64

set.seed(229)
rdata <- matrix(rnorm(ndata*nvars, mean=0, sd=2), nrow=ndata)
## 1: 2.9, 2: 5
bounds <- matrix(c(2.9, 5), nrow=ndata,ncol=2,byrow=T)
xbounds <- matrix(+Inf,nrow=ndata,ncol=nvars-2)
## censored <- cbind(rep(1,ndata), rep(0,ndata))
censored <- rdata[,1:2] >= bounds
xcensored <- rdata >= cbind(bounds,xbounds)
##
rdata2 <- rdata
rdata[xcensored] <- NA

cbind(rdata2,rdata,censored)
bounds[!is.na(rdata)] <- signs[!is.na(rdata)]*Inf


datapoints <- c(
    list(Rdata = rdata, Censored=censored*1L)
)

constants <- list(ndata=ndata, nclusters=nclusters, bounds=bounds, nvars=nvars)

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
        for(variate in 1:2){# real variates
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
                         list(Censored=c(ndata,2),
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
Cmodel$setInits(initsFunction())
output <- Cmcsampler$run(niter=1100, thin=1, reset=T, resetMV=TRUE, nburnin=100, time=T)
##
dpsamples <- as.matrix(Cmcsampler$mvSamples)
times <- Cmcsampler$getTimes()
names(times) <- sapply(confmodel$getSamplers(),function(x)x$target)

Cmodel$plotGraph()
