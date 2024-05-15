#### Example procedure for using the OPM

set.seed(201) # random seed

### First of all, clone the repository
### https://github.com/pglpm/bayes_nonparametric_inference
### or at least the folder 'package' therein
###
### Then make sure that you're in the folder 'package'
###
### The example datafile is:
datafile <- 'exampledata.csv'

## Load main functions and packages,
## including those for parallel processing
source('bnpi.R')

## For parallel processing (on Linux)
##
## Skip these lines if you don't want to use parallel computation
ncores <- 2 # number of CPUs to use
mycluster <- makeCluster(ncores, outfile="")
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
buildmetadata(data=alldata, file='meta-exampledata')

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
ntrain <- 100
trainpoints <- sort(sample(1:nrow(alldata), ntrain))

## Call the main function that does the Monte Carlo sampling / "training"
## note that we are feeding fewer datapoints
## Let's ask for 240 Monte Carlo samples
## from two Monte Carlo chains
## We must also specify an output directory

outputdir0 <- '_testexampledata2_Jef'
inferpopulation(data=alldata[trainpoints], metadata='meta-exampledata-modified_nopred.csv', outputdir=outputdir0)

## !!!! Check the directory name given in the messages of 'inferpopulation()'
outputdir <- paste0(outputdir0,'-V7-D',ntrain,'-K64-S',if(exists('nsamples')){nsamples}else{1200})



## The detailed Monte Carlo sampling can be monitored, for each core,
## by reading the files '_log-1.log', '_log-2.log', etc in the output dir
##
## The final output to be used in subsequent analysis is stored in
## [outputdir]/Fdistribution.rds
##
## Monte Carlo traces, to check convergence, are in
## [outputdir]/MCtraces.pdf
##
## Marginal distributions for each individual variate are in
## [outputdir]/plotsamples_Fdistribution.pdf
## [outputdir]/plotquantiles_Fdistribution.pdf

#### Example analysis using the updated probability distribution

## Let's say we want the posterior probability
## P(Truth=y | Output_class_0,=+1, Output_class_1=-1, ..., Output_class_5=-1)
## for y=0,...,5
##
## We call the predictand 'Y' and the predictor 'X'

## Predictand:
## we must define a matrix having 'Truth' as column
## and the specified x values as rows:
Y <- cbind(Truth=0:5)

## Predictor:
## we must define a matrix having 'Output_class_0,',...,'Output_class_5' as columns
## and the desired values as
X <- cbind(Output_class_0=+1,
           Output_class_1=-1,
           Output_class_2=-1,
           Output_class_3=-1,
           Output_class_4=-1,
           Output_class_5=-1)

## We use the function 'samplesFDistribution()'
## to calculate the posterior probability and its uncertainty.
## It requires the specification of the Monte Carlo output directory
## (don't worry about 'recycling' warnings)

testposterior <- samplesFDistribution(Y=Y, X=X, mcoutput=outputdir, parallel=TRUE)

## 'testposterior' has:
## - one row for each Y value
## - one column for each *population-frequency* f(Y|X)
##
## These limit frequencies are samples from
## the infinite-dimensional probability distribution over
## all possible population-frequencies, calculated from the data.
##
## So from these samples we can obtain:
## - the probability P(Y|X), which is just their mean
## - the variability of P(Y|X) if we had more data
## which can be calculated as quantiles or as a standard deviation
##
## In this case let's plot the plot the probability P(Y|X)
## and a sample of 100 population-frequencies,
## which give a visual indication of the variability

probdistr <- rowMeans(testposterior)
samplefreqs <- testposterior[,sample(1:ncol(testposterior), 100)]

pdff(paste0(outputdir,'/testresult1'))# open a pdf for plotting
## plot samples first
tplot(x=0:5, y=samplefreqs,
      xlab='Truth', ylab='probability',
      ylim=c(0,1),
      lty=1, # solid lines
      lwd=1, # thin
      col=7, # grey
      alpha=0.8) # quite transparent
##
## plot P(Y|X)
tplot(x=0:5, y=probdistr,
      type='b', # line+points
      lty=1, # solid
      lwd=3, # thicker
      col=1, # blue
      add=T) # add to previous plot
dev.off() # close pdf


## Now let's choose 100 datapoints from the full dataset
## excluding those used for training
ntest <- 100
testpoints <- sort(sample(setdiff(1:nrow(alldata), trainpoints), ntest))

## For each of these points
## we want to examine the posterior probabilities
## P(Truth=y | Output_class_0,..., ..., Output_class_5=...)
## for y=0,...,5
## and compare them with the prediction that was made by the ML algorithm
##
## We also want to compare this probability with the softmax
## obtained from the 'Output_class' weights
##
## And we compute the gain.
## Utility matrix gives +1 for correct choice, 0 otherwise
##
## We save all a plot for each datapoint as a pdf page
##
## We reuse the Y matrix from before,
## and for X we can use directly the data object 'alldata'
## but must select the variate columns first:
Xcols <- paste0('Output_class_', 0:5)

## Utility matrix
utilities <- diag(6)

## Preselect samples of population-frequencies to display
postsamples <- sample(1:ncol(testposterior), 100)

pdff(paste0(outputdir,'/testresult2'))
##
gain <- 0
gainsoftmax <- 0
##
surprise <- 0
surprisesoftmax <- 0
##
for(i in testpoints){
    ## get ML weights
    X <- alldata[i, ..Xcols] # note the tricky ".." syntax
    ## calculate softmax
    softmax <- exp(as.numeric(X))
    softmax <- softmax/sum(softmax)
    ## calculate posterior probability & samples
    posterior <- samplesFDistribution(Y=Y, X=X, mcoutput=outputdir, silent=TRUE, parallel=TRUE)
    probdistr <- rowMeans(posterior)
    samplefreqs <- posterior[, postsamples]
    ##
    ## make decision according to probability
    choice <- which.max(utilities%*%probdistr)
    choicesoftmax <- which.max(utilities%*%softmax)
    ##
    predictedvalue <- alldata[i,Predictions]+1
    truevalue <- alldata[i,Truth]+1
    gain <- gain + utilities[truevalue,choice]
    gainsoftmax <- gainsoftmax + utilities[truevalue,choicesoftmax]
    ##
    surprise <- surprise - log10(probdistr[truevalue])
    surprisesoftmax <- surprisesoftmax - log10(softmax[truevalue])
    ##
    ## plot samples first
    tplot(x=0:5, y=samplefreqs,
          xlab='Truth', ylab='probability', ylim=c(0,1),
          lty=1, lwd=1, col=7, alpha=0.8)
    ## plot P(Y|X)
    tplot(x=0:5, y=probdistr,
          type='b', lty=1, lwd=3, col=1, add=T)
    ## plot softmax
    tplot(x=0:5, y=softmax,
          type='b', lty=2, lwd=2, col=2, add=T)
    text(x=2.5, y=1, adj=0.5,
         labels=paste0('Datapoint ', i,
                       ' - Predicted: ', predictedvalue,
                       ', True: ',truevalue) )
    tplot(x=alldata[i,c('Truth','Predictions')], y=0.95,
          type='p', pch=c('T','P'), add=T)
}
dev.off()
##
cat('\nGain per sample:',gain/ntest,'\n')
cat('Softmax gain per sample:',gainsoftmax/ntest,'\n')
cat('\nSurprise per sample:',surprise/ntest,'\n')
cat('Softmax surprise per sample:',surprisesoftmax/ntest,'\n')

