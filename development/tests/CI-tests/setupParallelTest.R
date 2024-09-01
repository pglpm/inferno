library('inferno')
# Load internal function
cat(getwd())
currentDir <- getwd()
utilPath <- file.path(currentDir, '../../../R/util.R')
source(utilPath)

# Stop cores
stopCores <- function(workers, silent) {
    if (!is.logical(workers$cluster)) {
        on.exit(closeCoresOnExit(workers$cluster, silent))
    }
}

# Test that we can open cores and parse arguments correctly
testNcores <- function(args, check, silent) {
    N <- length(args)
    for (i in 1:N) {
        parallel <- args[i]
        correctNcores <- check[i]
        message <- paste0('Testing setupParallel() with argument ', parallel)
        printPretty(message, silent)
        workers <- setupParallel(parallel, silent)
        ncores <- workers$ncores
        stopCores(workers, silent)
        if (!(ncores == correctNcores)) {
            message <- paste0('WARNING: ', ncores, ' cores were set up, but ',
                              correctNcores, ' were requested.')
            printPretty(message, silent = FALSE)
        }
    }
}

# Test that we can open additional workers
testOpenWorkers <- function(silent) {
    # Open 2 workers
    workers <- setupParallel(2, silent)
    # Find workers again
    workers <- setupParallel(TRUE, silent)
    # Open 2 more workers
    workers <- setupParallel(4, silent)
    # Check if there are 4 workers
    ncores <- workers$ncores
    if ((workers$ncores) != 2) {
        message <- paste0('WARNING: ', ncores, ' cores were set up, but 
                          2 were expected.')
        printPretty(message, silent = FALSE)
    }
    stopCores(workers, silent)

}

silent <- FALSE
printPretty('', silent)
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