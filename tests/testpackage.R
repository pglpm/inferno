startdir <- getwd()

#### Check and change working directory if necessary
if(basename(startdir) == 'tests'){
  setwd('../package/')
  cat('\nSwitching to "package" directory\n')
} else if(basename(startdir) == 'bayes_nonparametric_inference'){
  setwd('package/')
  cat('\nSwitching to "package" directory\n')
} else if(basename(startdir) != 'package'){
  stop('Please run this script from directory "package" or "tests" or base')
}

source('bnpi.R')

testdir <- '../tests/'

refdir <- paste0(testdir,
                 'reference_packagetest_seed16_240601T0912-V8-D15-K64-S120/')

seed <- 16

outputdirPrefix <- paste0(paste0(testdir,'__packagetest'))

currenttestdir <- inferpopulation(data = paste0(testdir, 'testdata.csv'),
                        metadata = paste0(testdir, 'metatestdata.csv'),
                        outputdir = outputdirPrefix,
                        appendtimestamp = TRUE, appendinfo = TRUE,
                        nsamples = 120, nchains = 12,
                        cleanup = FALSE, parallel = 4,
                        seed = seed)

#### Test whether Fdistribution output is identical
cat('\nVerifying equality of "Fdistribution.rds" (TRUE = passed):\n')
print(identical(
  readRDS(paste0(currenttestdir,'Fdistribution.rds')),
  readRDS(paste0(refdir,'Fdistribution.rds'))
))

#### Test whether computation log-1 is identical
## remove lines containing diagnostic times, as they can vary
cat('\nVerifying equality of "log-1" (TRUE = passed):\n')
reffile <- readLines(paste0(refdir,'log-1.txt'))
reffile <- reffile[!grepl('time', reffile, fixed=TRUE)]
currentfile <- readLines(paste0(currenttestdir,'log-1.log'))
currentfile <- currentfile[!grepl('time', currentfile, fixed=TRUE)]

print(identical(currentfile, reffile))

#### return to original directory, in case called with 'source'
cat('\nSwitching to original directory\n')
setwd(startdir)

