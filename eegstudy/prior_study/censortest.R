library('nimble')

## Dirichlet-process mixture of product kernels
## Two kinds of variates, two variates per kind:
## real (kernel: gaussian, prior: normal + compound-gamma)
## categorical (kernel: categorical distr, prior: Dirichlet)

includecat <- TRUE # switch to include/exclude categorical components
ndata <- 7
nclusters <- 4

set.seed(223)
##
rdata <- matrix(rnorm(ndata*2, mean=0, sd=2), nrow=ndata)
bounds <- cbind(rep(0.5,ndata), rep(1,ndata))
censored <- rdata>=bounds
##
rdata[censored] <- NA

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
            Rvariances = matrix(rinvgamma(2*nclusters, shape=0.5, rate=1), nrow=2)
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
                                      'Rrates', 'Rmeans', 'Rvariances'
                           ))

mcsampler <- buildMCMC(confmodel)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

Cmodel$setInits(initsFunction())
output <- Cmcsampler$run(niter=100, thin=1, reset=T, resetMV=TRUE, nburnin=0)

dpsamples <- as.matrix(Cmcsampler$mvSamples)

