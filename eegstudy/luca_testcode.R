library('nimble')

## Dirichlet-process mixture of product kernels
## Two kinds of variates, two variates per kind:
## real (kernel: gaussian, prior: normal + compound-gamma)
## categorical (kernel: categorical distr, prior: Dirichlet)

includecat <- TRUE # switch to include/exclude categorical components
ndata <- 5
nclusters <- 10
## NOTE:
## setting ndata=3 and nclusters=8 leads to an error
## even without categorical variates:
##
## free(): invalid next size (fast)
## aborted (core dumped)

set.seed(10)


## Data (rows: datapoints, cols: variates)
datapoints <- c(
    list(Rdata = matrix(rnorm(ndata*2, mean=0, sd=1), nrow=ndata)),
    if(includecat){
        list(Cdata = matrix(sample(1:4, ndata*2, prob=rep(1, 4), replace=T), nrow=ndata))
    }
)

## Parameters of Dirichlet prior for Cprobs
Calphas = matrix(1, nrow=2, ncol=4)

constants <- list(ndata=ndata, nclusters=nclusters)

initsFunction <- function(){
    c(
        list(
            alpha = 1,
            Labels = rep(1, ndata), # all data to first cluster
            Rrates = rinvgamma(2, shape=0.5, scale=1),
            Rmeans = matrix(rnorm(2*nclusters, mean=0, sd=3), nrow=2),
            Rvariances = matrix(rinvgamma(2*nclusters, shape=0.5, rate=1), nrow=2)
            ),
        if(includecat){
            list( Calphas = Calphas,
                 Cprobs = aperm(sapply(1:2, # loop over variates
                                       function(i){
                                           sapply(1:nclusters, # loop over clusters
                                                  function(cluster){
                                                      rdirch(1, alpha=Calphas[i,])
                                                  })
                                       }, simplify='array'
                                       ))
                 )
        }
    )
}


infmixture <- nimbleCode({
    Labels[1:ndata] ~ dCRP(alpha, size=ndata)
    ##
    ## the nested inverse-gammas are equivalent to a compound-gamma for the variance
    for(variate in 1:2){
        Rrates[variate] ~ dinvgamma(shape=0.5, scale=1)
    }
    ##
    for(cluster in 1:nclusters){
        for(variate in 1:2){# real variates
            Rmeans[variate, cluster] ~ dnorm(mean=0, var=3*3)
            Rvariances[variate, cluster] ~ dinvgamma(shape=0.5, rate=Rrates[variate])
        }
        if(includecat){# categorical variates
            for(variate in 1:2){
                Cprobs[variate, cluster, 1:4] ~ ddirch(alpha=Calphas[variate, 1:4])
            }
        }
    }
    ##
    for(datum in 1:ndata){
        for(variate in 1:2){# real variates
            Rdata[datum, variate] ~ dnorm(mean=Rmeans[variate, Labels[datum]], var=Rvariances[variate, Labels[datum]])
        }
        if(includecat){# categorical variates
            for(variate in 1:2){
                    Cdata[datum, variate] ~ dcat(prob=Cprobs[variate, Labels[datum], 1:4])
            }
        }
    }
})

model <- nimbleModel(code=infmixture, name='model',
                     constants=constants,
                     inits=initsFunction(),
                     data=datapoints,
                     dimensions=c(
                         list(Labels=ndata,
                              Rrates=2,
                              Rmeans=c(2, nclusters),
                              Rvariances=c(2, nclusters)),
                         if(includecat){
                             list(Cprobs=c(2, nclusters, 4))
                         }
                     )
                     )
Cmodel <- compileNimble(model, showCompilerOutput=FALSE)


confmodel <- configureMCMC(Cmodel,
                           monitors=c('Labels',
                                      'Rrates', 'Rmeans', 'Rvariances',
                                      if(includecat){ 'Cprobs' } )
                           )

mcsampler <- buildMCMC(confmodel)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

Cmodel$setInits(initsFunction())
output <- Cmcsampler$run(niter=100, thin=1, reset=T, resetMV=TRUE, nburnin=0)

dpsamples <- getSamplesDPmeasure(Cmcsampler)
