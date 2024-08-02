startdir <- getwd()

#### Check and change working directory if necessary
if(basename(startdir) != 'test_custom'){
  cat('\nAre you in the correct folder?\n')
}

cat('\nInstalling local package "modelfreeinference" in local library\n\n')

if(!(Sys.getenv("R_LIBS_USER") %in% .libPaths())) {
  stop('Make sure your local installation directory,\n',
       Sys.getenv("R_LIBS_USER"),
       '\nexists.\n')
}

devtools::load_all()
## library('modelfreeinference')


seed <- 16

outputdirPrefix <- file.path('__packagetest')

currenttestdir <- inferpopulation(data = NULL,
                        metadata = 'metadata_test_custom.csv',
                        outputdir = outputdirPrefix,
                        output = 'directory',
                        appendtimestamp = TRUE,
                        appendinfo = TRUE,
                        nsamples = 120, nchains = 1,
                        cleanup = FALSE, parallel = 1,
                        prior = TRUE,
                        ## lldata = 12,
                        showKtraces = T,
                        showAlphatraces = T,
                        seed = seed)

if(FALSE){
  ## The maths of the new package version has changed a little,
  ## so a comparison with old-version results are not meaningful
refdir <- 'reference_seed16-vrt9_dat15_smp120'

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
}

