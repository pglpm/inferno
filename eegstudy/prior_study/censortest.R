library('nimble')

## Dirichlet-process mixture of product kernels
## Two kinds of variates, two variates per kind:
## real (kernel: gaussian, prior: normal + compound-gamma)
## categorical (kernel: categorical distr, prior: Dirichlet)

includecat <- TRUE # switch to include/exclude categorical components
ndata <- 7
nclusters <- 4

set.seed(229)
##
rdata <- matrix(rnorm(ndata*2, mean=0, sd=2), nrow=ndata)
bounds <- cbind(rep(-0.5,ndata), rep(1,ndata))
signs <- cbind(rep(-1,ndata),rep(1,ndata))
## censored <- cbind(rep(1,ndata), rep(0,ndata))
censored <- rdata>=bounds
##
rdata2 <- rdata
rdata[signs*rdata>=signs*bounds] <- NA
cbind(rdata2,rdata,censored)
bounds[!is.na(rdata)] <- signs[!is.na(rdata)]*Inf


datapoints <- c(
    list(Rdata = rdata, Censored=censored*1L)
)

constants <- list(ndata=ndata, nclusters=nclusters, bounds=bounds)

initsFunction <- function(){
    c(
        list(
            palphas = c(1/3,1/3,1/3),
            qalpha0 = rep(0.25/nclusters, nclusters), # cluster probabilities
            alpha0 = 2,
            q=rep(1/nclusters,nclusters),
            Labels = rep(1, ndata), # all data to first cluster
            Rrates = rinvgamma(2, shape=0.5, scale=1),
            Rmeans = matrix(rnorm(2*nclusters, mean=0, sd=3), nrow=2),
            Rvariances = matrix(rinvgamma(2*nclusters, shape=0.5, rate=1), nrow=2),
            'Rdata[4, 1]'=-2,
            'Rdata[3, 2]'=2
        )
    )
}

infmixture <- nimbleCode({
    alpha0 ~ dcat(prob=palphas[1:3])
    qalpha[1:nclusters] <- 2^alpha0 * qalpha0[1:nclusters]
    q[1:nclusters] ~ ddirch(alpha=qalpha[1:nclusters])
    ##
    for(variate in 1:2){
        Rrates[variate] ~ dinvgamma(shape=0.5, scale=1)
    }
    for(cluster in 1:nclusters){
        for(variate in 1:2){# real variates
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
            Rdata[datum, variate] ~ dnorm(mean=Rmeans[variate, Labels[datum]], var=Rvariances[variate, Labels[datum]])
        }
    }
})

model <- nimbleModel(code=infmixture, name='model',
                     constants=constants,
                     inits=initsFunction(),
                     data=datapoints,
                     dimensions=c(
                         list(Labels=ndata,
                              q=nclusters,
                              Rrates=2,
                              Rmeans=c(2, nclusters),
                              Rvariances=c(2, nclusters))
                     )
                     )

Cmodel <- compileNimble(model, showCompilerOutput=FALSE)


confmodel <- configureMCMC(Cmodel,
                           monitors=c('Labels','q','alpha0',
                                      'Rdata[2, 1]',
                                      'Rdata[4, 2]',
                                      'Rrates', 'Rmeans', 'Rvariances'
                           ))

mcsampler <- buildMCMC(confmodel)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)
##
Cmodel$setInits(initsFunction())
output <- Cmcsampler$run(niter=10000, thin=1, reset=T, resetMV=TRUE, nburnin=0, time=T)
##
dpsamples <- as.matrix(Cmcsampler$mvSamples)
times <- Cmcsampler$getTimes()
names(times) <- sapply(confmodel$getSamplers(),function(x)x$target)

Cmodel$plotGraph()
