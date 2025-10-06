library('inferno')

seed <- 16
parallel <- 4

outputdir <- '__test_custom500'
learntdir <- learn(
    data = 'dataset_custom500.csv',
    metadata = 'metadata_custom.csv',
    prior = FALSE,
    outputdir = outputdir,
    appendtimestamp = TRUE,
    appendinfo = TRUE,
    outputvalue = 'directory',
    parallel = parallel,
    ## parameters for short test run:
    ## subsampledata = 10,
    ## maxhours = 0,
    ## nsamplesperchain = 60,
    ## nchains = parallel + 1,
    ##
    ## hyperparams = list(
    ##     ncomponents = 64,
    ##     minalpha = -4,
    ##     maxalpha = 4,
    ##     byalpha = 1,
    ##     Rshapelo = 0.5,
    ##     Rshapehi = 0.5,
    ##     Rvarm1 = 3^2,
    ##     Cshapelo = 0.5,
    ##     Cshapehi = 0.5,
    ##     Cvarm1 = 3^2,
    ##     Dshapelo = 0.5,
    ##     Dshapehi = 0.5,
    ##     Dvarm1 = 3^2,
    ##     Lshapelo = 0.5,
    ##     Lshapehi = 0.5,
    ##     Lvarm1 = 3^2,
    ##     Bshapelo = 1,
    ##     Bshapehi = 1,
    ##     Dthreshold = 1,
    ##     tscalefactor = 2,
    ##     initmethod = 'prior'
    ##     ## precluster, prior
    ## ),
    seed = seed
)

## mi <- mutualinfo(
##     Y1names = c('N2vrt'),
##     Y2names = c('Rvrt'),
##     X = cbind(Bvrt = 'no'),
##     learnt = currenttestdir,
##     nsamples = 3600,
##     parallel = 4
## )
## 
## print(mi)
## 
## warnings()
## 
## dataset <- read.csv('dataset_custom500.csv', na.strings='')
## nv <- round(ncol(dataset)/2)
## probs <- Pr(
##     Y = dataset[1:20,1:nv,drop=F],
##     X = dataset[1:20,(nv+1):ncol(dataset), drop=F],
##     learnt = currenttestdir,
##     parallel = 4
## )
## 
## print(probs$values)
## 
## print(probs$quantiles)
## 
## warnings()
## 
## 
## cat('\nEnd\n')
