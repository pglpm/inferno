nsamples <- 4
##
W <- c(0.5, 0.5)
dim(W) <- length(W)
Bprob <- array(c(0, 1), dim = c(1, 2))
Nprob <- array(rbind(c(0.5, 0.5, 0), c(0,0,1)), dim = c(1, 2, 3))
##
od <- dim(W)
W <- rep(W, nsamples)
dim(W) <- c(od, nsamples)
##
od <- dim(Bprob)
Bprob <- rep(Bprob, nsamples)
dim(Bprob) <- c(od, nsamples)
##
od <- dim(Nprob)
Nprob <- rep(Nprob, nsamples)
dim(Nprob) <- c(od, nsamples)
##
auxmetadata <- data.frame(
    name = c('B1', 'N1'),
    mcmctype = c('B', 'N'),
    id = c(1, 1),
    transform = c('identity', 'identity'),
    Nvalues = c(2, 3),
    halfstep = c(NA, NA),
    domainmin = c(NA, NA),
    domainmax = c(NA, NA),
    tdomainmin = c(0, 1),
    tdomainmax = c(1, 3),
    leftbound = c(NA, NA),
    rightbound = c(NA, NA),
    tleftbound = c(NA, NA),
    trightbound = c(NA, NA),
    tlocation = c(0, 0),
    tscale = c(1, 1),
    plotmin = c(NA, NA),
    plotmax = c(NA, NA),
    V1 = c('yes', 'A'),
    V2 = c('no', 'B'),
    V3 = c(NA, 'C')
)
##
testlearned <- list(
    Bprob = Bprob,
    Nprob = Nprob,
    W = W,
    auxmetadata = auxmetadata)
##
saveRDS(testlearned, 'learned.rds')

