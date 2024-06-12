startdir <- getwd()

#### Check and change working directory if necessary
if(basename(startdir) == 'tests'){
  setwd(file.path('..','R'))
  cat('\nSwitching to "R" directory\n')
} else if(basename(startdir) == 'bayes_nonparametric_inference'){
  setwd('R')
  cat('\nSwitching to "R" directory\n')
} else if(basename(startdir) != 'R'){
  stop('Please run this script from directory "R" or "tests" or base')
}

source('bnpi.R')

testdir <- file.path('..', 'tests')

refdir <- file.path(testdir,
                 'reference_packagetest_seed16_240601T0912-V8-D15-K64-S120')

seed <- 16

outputdirPrefix <- file.path(testdir,'__packagetest')

currenttestdir <- inferpopulation(data = file.path(testdir, 'testdata.csv'),
                        metadata = file.path(testdir, 'metatestdata.csv'),
                        outputdir = outputdirPrefix,
                        output = 'directory',
                        appendtimestamp = TRUE, appendinfo = TRUE,
                        nsamples = 120, nchains = 12,
                        cleanup = FALSE, parallel = 4,
                        seed = seed)

#### Test whether Fdistribution output is identical
cat('\nVerifying equality of "Fdistribution.rds" (TRUE = passed):\n')
print(identical(
  readRDS(file.path(currenttestdir,'Fdistribution.rds')),
  readRDS(file.path(refdir,'Fdistribution.rds'))
))

#### Test whether computation log-1 is identical
## remove lines containing diagnostic times, as they can vary
cat('\nVerifying equality of "log-1" (TRUE = passed):\n')
reffile <- readLines(file.path(refdir,'log-1.txt'))
reffile <- reffile[!grepl('time', reffile, fixed=TRUE)]
currentfile <- readLines(file.path(currenttestdir,'log-1.log'))
currentfile <- currentfile[!grepl('time', currentfile, fixed=TRUE)]

print(identical(currentfile, reffile))

#### return to original directory, in case called with 'source'
cat('\nSwitching to original directory\n')
setwd(startdir)

