library('inferno')
# Load internal function
cat(getwd())
currentDir <- getwd()
utilPath <- file.path(currentDir, '../../../R/util.R')
source(utilPath)

# Stop cores
stopCores <- function(workers, silent) {
    if (!is.logical(workers$cluster)) {
        on.exit(closecoresonexit(workers$cluster, silent))
    }
}

# Test that we can open cores and parse arguments correctly
testNcores <- function(args, check, silent) {
    N <- length(args)
    for (i in 1:N) {
        parallel <- args[i]
        correctNcores <- check[i]
        message <- paste0('Testing setup_parallel() with argument ', parallel)
        print_pretty(message, silent)
        workers <- setup_parallel(parallel, silent)
        ncores <- workers$ncores
        stopCores(workers, silent)
        if (!(ncores == correctNcores)) {
            message <- paste0('WARNING: ', ncores, ' cores were set up, but ',
                              correctNcores, ' were requested.')
            print_pretty(message, silent = FALSE)
        }
    }
}

# Test that we can open additional workers
testOpenWorkers <- function(silent) {
    # Open 2 workers
    workers <- setup_parallel(2, silent)
    # Find workers again
    workers <- setup_parallel(TRUE, silent)
    # Open 2 more workers
    workers <- setup_parallel(4, silent)
    # Check if there are 4 workers
    ncores <- workers$ncores
    if ((workers$ncores) != 2) {
        message <- paste0('WARNING: ', ncores, ' cores were set up, but 
                          2 were expected.')
        print_pretty(message, silent = FALSE)
    }
    stopCores(workers, silent)

}

silent <- FALSE
print_pretty('', silent)
# Test integer input
argList <- c(1, 2, 4)
ncoresCheck <- c(1, 2, 4)
testNcores(argList, ncoresCheck, silent)
# Test boolean input
argList <- c(FALSE)
ncoresCheck <- c(1)
testNcores(argList, ncoresCheck, silent)
# Test opening nodes sequentially
testOpenWorkers(silent)
#