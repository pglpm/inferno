library('inferno')

seed <- 16
parallel <- 4

outputdir <- '__tsf2_testprior_allvrt'
learntdir <- learn(
    data = NULL,
    metadata = 'metadata_allvrt_prior.csv',
    prior = FALSE,
    outputdir = outputdir,
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    output = 'directory',
    parallel = parallel,
    maxhours = 0,
    ## parameters for short test run:
    ## subsampledata = 10,
    ## maxhours = 0,
    ## nsamplesperchain = 60,
    ## nchains = parallel + 1,
    ##
    hyperparams = list(
        ncomponents = 64,
        minalpha = -4,
        maxalpha = 4,
        byalpha = 1,
        Rshapelo = 0.5,
        Rshapehi = 0.5,
        Rvarm1 = 3^2,
        Cshapelo = 0.5,
        Cshapehi = 0.5,
        Cvarm1 = 3^2,
        Dshapelo = 0.5,
        Dshapehi = 0.5,
        Dvarm1 = 3^2,
        Lshapelo = 0.5,
        Lshapehi = 0.5,
        Lvarm1 = 3^2,
        Bshapelo = 1,
        Bshapehi = 1,
        Dthreshold = 1,
        tscalefactor = 2,
        avoidzeroW = FALSE,
        initmethod = 'allinone'
        ## precluster, prior
    ),
    seed = seed
)



pRvrt <- Pr(Y = data.frame(Rvrt = seq(-3, 3, length.out = 129)),
    learnt = '/home/pglpm/repos/inferno/development/tests/priors/__tsf2_testprior_allvrt-250627T090902-vrt9_dat0_smp3600/learnt.rds',
    parallel = TRUE)


mypdf('prior_Rvrt', apaper = 1)
par(mfcol = c(20, 20))
for(asample in round(seq(1, dim(pRvrt$samples)[3], length.out = prod(par()$mfcol)))) {
    par(mar = rep(1, 4))
    flexiplot(
        x = seq(-3, 3, length.out = 129),
        y = pRvrt$samples[,1,asample],
        xlab = NA, ylab = NA,
        xdomain = FALSE, ydomain = FALSE
    )
}
dev.off()
