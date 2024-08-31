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
            print_pretty(message, silent)
        }
    }
}

# Test that we can open additional workers
testOpenWorkers <- function(silent) {
    # Open 2 workers
    workers <- setup_parallel(2, silent)
    # Open 2 more workers
    workers <- setup_parallel(2, silent)
    # Check if there are 4 workers
    ncores <- workers$ncores
    if (!(workers$ncores) != 4) {
        message <- paste0('WARNING: ', ncores, ' cores were set up, but 
                          4 were requested.')
        print_pretty(message, silent)
    }
    stopCores()

}

argList <- c(0, 1, 2, 4, TRUE, FALSE)
ncoresCheck <- c(1, 1, 2, 4, 1, 1)
silent <- TRUE
testNcores(argList, ncoresCheck, silent)
