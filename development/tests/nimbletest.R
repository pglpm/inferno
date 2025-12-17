library('nimble')

constants <- list(
    ncomponents = 32,
    npoints = 5,
    nalpha = 8,
    alphabase = sqrt(2),
    probalpha0 = (1:8)^2.25 / sum((1:8)^2.25),
    dirchalphas = rep((2^(-4 - 0.5)) / 32, 32),
    Nn = 1,
    Ncards = 5,
    Ni = c(1, NA_integer_),
    Nf = c(5, NA_integer_),
    Nalpha0 = rep(0.2, 5)
)

datapoints <- list(
    Ndata = matrix(1:5, nrow = 5, ncol = 1)
)


initsfn <- function() { list(
    Alpha = 4,
    W = rep(1/32, 32),
    K = 1:5,
    Nprob = matrix(0.2, nrow = 5, ncol = 32)
    ) }

finitemix <- nimbleCode({
    ## Component weights
    Alpha ~ dcat(prob = probalpha0[1:nalpha])
    alphas[1:ncomponents] <- dirchalphas[1:ncomponents] * alphabase^Alpha
    W[1:ncomponents] ~ ddirch(alpha = alphas[1:ncomponents])
    ## Probability density for the parameters of the components
    for (k in 1:ncomponents) {
        for (v in 1:Nn) {
            ## NB: in this case `v` has only one value,
            ## but in the general case it is used to select
            ## different chunks of `Nprob`, via `Ni` and `Nf`
            Nprob[Ni[v]:Nf[v], k] ~ ddirch(alpha = Nalpha0[Ni[v]:Nf[v]])
        }
    }
    ## Probability of data
    for (d in 1:npoints) {
        K[d] ~ dcat(prob = W[1:ncomponents])
        for (v in 1:Nn) {
            Ndata[d, v] ~ dcat(prob = Nprob[Ni[v]:Nf[v], K[d]])
        }
    }
})

finitemixnimble <- nimbleModel(
    code = finitemix,
    name = 'finitemixnimble1',
    constants = constants,
    data = datapoints,
    inits = initsfn()
)

Cfinitemixnimble <- compileNimble(finitemixnimble,
    showCompilerOutput = FALSE)

confnimble <- configureMCMC(
    Cfinitemixnimble, # nodes = NULL
    monitors = c( 'W', 'Nprob' ),
    monitors2 = 'K'
)

mcsampler <- buildMCMC(confnimble)
## print(confnimble$getUnsampledNodes())
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)

Cfinitemixnimble$setInits(initsfn())

Cmcsampler$run(
    niter = 3600,
    thin = 1,
    thin2 = 1800,
    nburnin = 0,
    time = FALSE,
    reset = TRUE,
    resetMV = TRUE
)

mcsamples <- as.list(Cmcsampler$mvSamples, iterationAsLastIndex = TRUE)
mcsamplesKA <- as.list(Cmcsampler$mvSamples2, iterationAsLastIndex = FALSE)
