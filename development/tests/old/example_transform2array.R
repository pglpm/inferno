#### Example procedure for using the OPM
startdir <- getwd()

#### Check and change working directory if necessary
if (basename(startdir) == 'tests') {
  setwd(file.path('..', 'R'))
  cat('\nSwitching to "R" directory\n')
} else if (basename(startdir) == 'bayes_nonparametric_inference') {
  setwd('R/')
  cat('\nSwitching to "R" directory\n')
} else if (basename(startdir) != 'R') {
  stop('Please run this script from directory "R" or "tests" or base')
}

set.seed(201) # random seed

### First of all, clone the repository
### https://github.com/pglpm/bayes_nonparametric_inference
### or at least the folder 'R' therein
###
### Then make sure that you're in the folder 'R'
###
### The example datafile is:
datafile <- file.path('..', 'tests', 'exampledata.csv')

## Load main functions and packages,
## including those for parallel processing
source('bnpi.R')

## For parallel processing (on Linux)
##
## Skip these lines if you don't want to use parallel computation
ncores <- 2 # number of CPUs to use
mycluster <- makeCluster(ncores, outfile = '')
registerDoParallel(mycluster)

## Read *all* available datapoints
alldata <- fread(datafile)

## Build a preliminary metadata file
##
## We use *all* datapoints,
## even those that can't be used for training because of computational reasons.
## This way our prior information contains at least some extra info
## that is available in the full dataset.
##
## This is especially important for the 'centralvalue'-'highvalue' metadata
buildmetadata(data = alldata, file = file.path('..', 'tests',
                                               'meta-exampledata'))

### Open the metadata .csv file and modify the metadata as appropriate
### In this case we do the following changes:
###
### 1. 'Truth' & 'Predictions' to type=nominal because they have no intrinsic ordering
### and create extra columns 'V1'-'V6' containing the variate values
###
### 2. 'domainmin/max' to +-Infinity for the continuous variates
### 'min/maxincluded' are already set to FALSE
###
### (See 'meta-exampledata-modified.csv' for the resulting modifications)


## Select a subset of data for training
ntrain <- 50
trainpoints <- sort(sample(1:nrow(alldata), ntrain))

## Call the main function that does the Monte Carlo sampling / 'training'
## note that we are feeding fewer datapoints
## Let's ask for 240 Monte Carlo samples
## from two Monte Carlo chains
## We must also specify an output directory

outputdir <- file.path('..', 'tests', '_testexampledata2')
learnt <- inferpopulation(data = alldata[trainpoints],
                           metadata = file.path('..', 'tests',
                                                'meta-exampledata-modified.csv'),
                           outputdir = outputdir,
                           nsamples = 240, nchains = 2,
                           appendtimestamp = FALSE,
                           appendinfo = FALSE)


## Now let's choose 100 datapoints from the full dataset
## excluding those used for training
ntest <- 100
testpoints <- sort(sample(setdiff(1:nrow(alldata), trainpoints), ntest))

## For each of these points
## we want to calculate and save the sample frequencies
## F(Truth=y | Output_class_0,..., ..., Output_class_5=...)
## for y=0,...,5
## and save them in a 3D array

## X variates
Xcols <- paste0('Output_class_', 0:5)
## Y values: each is repeated as many times as the number of testpoints
repY <- cbind(Truth = rep(0:5, each=ntest))
## repY and alldata[testpoints, ..Xcols] look like this, side-by-side:
## > cbind(repY, alldata[testpoints, ..Xcols])
##      Truth Output_class_0 Output_class_1 Output_class_2 Output_class_3 Output_class_4 Output_class_5
##      <int>          <num>          <num>          <num>          <num>          <num>          <num>
##   1:     0      -8.386725       1.431273      -1.481002       1.826596       -7.28851       1.820939
##   2:     0     -10.049677       1.010583      -3.203500       1.474996       -3.64449       1.766000
##   3:     0      -8.288185       0.449368      -5.052208       0.243640        2.70021      -0.268845
##   4:     0      -5.115926       1.345078      -2.022172       1.206520       -7.06938       1.321461
##   5:     0       2.027061       0.301541      -3.107857      -0.993193       -1.60768      -1.465528
##  ---                                                                                                
## 596:     5      -5.144355       0.200953       1.617125       0.813927       -9.50161       0.807825
## 597:     5     -16.285881      -0.268113      -5.206679       0.719011        3.20312       0.731005
## 598:     5       1.843565       0.875414      -1.588141      -0.113086       -7.58015      -0.448277
## 599:     5     -19.162981      -0.447744      -5.888669       0.134154        4.70466      -0.540910
## 600:     5       0.477277       0.876158      -0.348632       0.125635       -8.15816       0.033434



condfreqs <- samplesFDistribution(Y = repY,
                                  X = alldata[testpoints, ..Xcols],
                                  learnt = outputdir,
                                  silent = TRUE, parallel = TRUE)

## condfreqs is a matrix having
## one row for each of the combination of Y & X values shown above (600 in total)
## and the frequency samples as 2nd dim
##
## We must reshape it into a 3D array having:
## 1st dim: id of test datapoint
## 2nd dim: Truth=0:5
## 3rd dim: sample frequencies (as many as MC samples)
dim(condfreqs) <- c(ntest, 6, 240)
dimnames(condfreqs) <- list(datumid=1:ntest, Truth=0:5, sampleF=1:240)

pdff('_testplot')
for(i in 1:20){
  tplot(x=0:5, y=condfreqs[i,,], xlab='Truth', ylab='freq',
        ylim=c(0,1),
        col=7, lwd=1, lty=1, alpha=0.75)
  tplot(x=0:5, y=rowMeans(condfreqs[i,,]), col=1, add=T)
}
dev.off()


## Here is how this 3D array looks like:
##
## > condfreqs[1:2,,1:3]
## , , 1
## 
##         Y
## testdata         0        1         2        3        4        5
##     [1,] 0.0865464 0.281950 0.0566617 0.103286 0.216499 0.255057
##     [2,] 0.0325579 0.312384 0.0196845 0.151457 0.282931 0.200985
## 
## , , 2
## 
##         Y
## testdata         0        1         2        3         4        5
##     [1,] 0.0474548 0.105506 0.0150629 0.483467 0.0612551 0.287255
##     [2,] 0.0474756 0.104889 0.0144944 0.484677 0.0614629 0.287001
## 
## , , 3
## 
##         Y
## testdata          0        1         2        3        4        5
##     [1,] 0.00119039 0.130986 0.0245917 0.449805 0.151961 0.241466
##     [2,] 0.00120722 0.122206 0.0201994 0.459026 0.155482 0.241880
## 
## > dimnames(condfreqs) <- list(datumid=1:ntest, Truth=0:5, sampleF=1:240)
## > condfreqs[1:2,,1:3]
## , , sampleF = 1
## 
##        Truth
## datumid         0        1         2        3        4        5
##       1 0.0865464 0.281950 0.0566617 0.103286 0.216499 0.255057
##       2 0.0325579 0.312384 0.0196845 0.151457 0.282931 0.200985
## 
## , , sampleF = 2
## 
##        Truth
## datumid         0        1         2        3         4        5
##       1 0.0474548 0.105506 0.0150629 0.483467 0.0612551 0.287255
##       2 0.0474756 0.104889 0.0144944 0.484677 0.0614629 0.287001
## 
## , , sampleF = 3
## 
##        Truth
## datumid          0        1         2        3        4        5
##       1 0.00119039 0.130986 0.0245917 0.449805 0.151961 0.241466
##       2 0.00120722 0.122206 0.0201994 0.459026 0.155482 0.241880

