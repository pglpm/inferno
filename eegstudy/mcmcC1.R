library('nimble')
##
## Dirichlet-process mixture of product kernels
## Three kinds of variates, two variates per kind:
## real (kernel: gaussian, prior: normal + compound-gamma)
## binary (kernel: bernoulli, prior: beta)
## categorical (kernel: categorical distr, prior: Dirichlet)

includecat <- TRUE # switch to include/exclude categorical components
ndata <- 3
nclusters <- 8

## Data (rows: datapoints, cols: variates)
datapoints <- c(
    list(Rdata = matrix(rnorm(ndata*2, mean=0, sd=1), nrow=ndata)),
    list(Bdata = matrix(sample(0:1, ndata*2, replace=T), nrow=ndata)),
    if(includecat){
        list(Cdata = cbind( sample(1:3, ndata, prob=rep(1, 3), replace=T),
                           sample(1:4, ndata, prob=rep(1, 4), replace=T) ))
    }
)

## Alphas for the Dirichlet prior of the categorical components.
## The first categorical variate has 3 possible categories,
## the second has 4 possible categories.
## To be able to use both in Nimble
## I set the number of categories to 4 for each,
## but assign an almost-zero prior alpha to the 4th category of the first:
Calphas = matrix(c(
    1, 1, 1, 2^(-40),
    1, 1, 1, 1
), byrow=T, nrow=2)

constants <- list(ndata=ndata, nclusters=nclusters)

initsFunction <- function(){
    c(
        list(
            alpha = 1,
            Labels = rep(1, ndata), # all data to first cluster
            Rrates = rinvgamma(2, shape=0.5, scale=1),
            Rmeans = matrix(rnorm(2*nclusters, mean=0, sd=3), nrow=2),
            Rvariances = matrix(rinvgamma(2*nclusters, shape=0.5, rate=1), nrow=2),
            Bprobs = matrix(rbeta(2*nclusters, shape1=1, shape2=1), nrow=2)
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
        for(variate in 1:2){# binary variates
            Bprobs[variate, cluster] ~ dbeta(shape1=1, shape2=1)
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
        for(variate in 1:2){# binary variates
            Bdata[datum, variate] ~ dbern(prob=Bprobs[variate, Labels[datum]])
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
                              Rvariances=c(2, nclusters),
                              Bprobs=c(2, nclusters) ),
                         if(includecat){
                             list(Cprobs=c(2, nclusters, 4))
                         }
                     )
                     )
Cmodel <- compileNimble(model, showCompilerOutput=FALSE)


confmodel <- configureMCMC(Cmodel,
                           monitors=c('Labels',
                                      'Rrates', 'Rmeans', 'Rvariances',
                                      'Bprobs',
                                      if(includecat){ 'Cprobs' } )
                           )

mcsampler <- buildMCMC(confmodel)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

set.seed(22)
Cmodel$setInits(initsFunction())
Cmcsampler$run(niter=100, thin=1, reset=T, resetMV=TRUE, nburnin=0)

dpsamples <- getSamplesDPmeasure(Cmcsampler)
