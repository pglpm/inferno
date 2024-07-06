startdir <- getwd()

#### Check and change working directory if necessary
if(basename(startdir) != 'tests'){
  cat('\nAre you in the correct folder?\n')
}

library('devtools')
load_all()

refdir <- 'reference_seed16_packagetest-240705T064934-vrt8_dat15_smp120'

seed <- 16

outputdirPrefix <- file.path('__packagetest')

currenttestdir <- inferpopulation(data = 'testdata.csv',
                        metadata = 'metatestdata.csv',
                        outputdir = outputdirPrefix,
                        output = 'directory',
                        appendtimestamp = FALSE, appendinfo = FALSE,
                        nsamples = 120, nchains = 12,
                        cleanup = FALSE, parallel = 4,
                        ## prior = TRUE,
                        ## lldata = 12,
                        useOquantiles = FALSE,
                        seed = seed)

#### Test whether Fdistribution output is identical
cat('\nVerifying equality of "Fdistribution.rds" (TRUE = passed):\n')
print(identical(
  readRDS(file.path(currenttestdir,'Fdistribution.rds')),
  readRDS(file.path(refdir,'Fdistribution.rds'))
))

#### Test whether MCtraces output is identical
cat('\nVerifying equality of "MCtraces.rds" (TRUE = passed):\n')
print(identical(
  readRDS(file.path(currenttestdir,'MCtraces.rds')),
  readRDS(file.path(refdir,'MCtraces.rds'))
))

#### Test whether computation log-1 is identical
## remove lines containing diagnostic times, as they can vary
cat('\nVerifying equality of "log-1" (TRUE = passed):\n')
reffile <- readLines(file.path(refdir,'log-1.txt'))
reffile <- reffile[!grepl('time', reffile, fixed=TRUE)]
currentfile <- readLines(file.path(currenttestdir,'log-1.log'))
currentfile <- currentfile[!grepl('time', currentfile, fixed=TRUE)]

print(identical(currentfile, reffile))


