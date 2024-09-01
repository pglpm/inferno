#' Determine the status of parallel processing. If cores are requested,
#' set up the cores. If cores are registered, keep using those.
#' @param parallel logical or integer: whether to use pre-existing parallel
#'   workers, or how many to create and use. Default `TRUE`.
#' @param silent logical: give warnings or updates in the computation?
#'   Default `FALSE`.
#' @keywords internal
#' @return list (ncores, cluster) ncores is the number of cores and
#'   cluster object is a parallel::makeCluster(parallel) object or FALSE
#'   if there is no cluster.
setupParallel <- function(parallel, silent = FALSE) {
    # Set cluster object to false to start
    cluster <- FALSE
    # Check if 'parallel' argument is either logical or numeric
    if (!(is.logical(parallel) || (is.numeric(parallel) &&
        parallel > 0 && (round(parallel) == parallel)))) {
        stop('Argument parallel must be TRUE, FALSE or a positive integer.')
    }

    ## ## Alternative way to register cores;
    ## ## might need to be used for portability to Windows?
    ## registerDoSEQ()
    ## cluster <- makePSOCKcluster(ncores)

    if (isTRUE(parallel)) {
        if (foreach::getDoParRegistered()) {
            printPretty(paste0('Using already registered ', foreach::getDoParName(),
                           ' with ', foreach::getDoParWorkers(), ' workers'), silent)
            ncores <- foreach::getDoParWorkers()
        } else {
            printPretty('No parallel backend registered.', silent)
            ncores <- 1
        }
    } else if (parallel >= 2) {
        if (foreach::getDoParRegistered()) {
            ncores <- min(foreach::getDoParWorkers(), parallel)
            printPretty(paste0('Using already registered ', foreach::getDoParName(),
                        ' with ', foreach::getDoParWorkers(), ' workers'), silent)
            if (parallel > ncores) {
                message <- paste0('NOTE: there are pre-registered cores,
                    but less than requested in the "parallel" argument. 
                    Running with ', ncores, ' cores')
                printPretty(message, silent)
            }
        } else {
            cluster <- parallel::makeCluster(parallel)
            doParallel::registerDoParallel(cluster)
            printPretty(paste0('Registered ', foreach::getDoParName(),
                    ' with ', foreach::getDoParWorkers(), ' workers'), silent)
            ncores <- parallel
        }
    } else {
        printPretty('No parallel backend registered.', silent)
        ncores <- 1
    }
    workers <- list('ncores' = ncores, 'cluster' = cluster)
    return(workers)
}

closeCoresOnExit <- function(cluster, silent = FALSE) {
    printPretty('Closing connections to cores.', silent)
    foreach::registerDoSEQ()
    parallel::stopCluster(cluster)
    env <- foreach:::.foreachGlobals
    rm(list = ls(name = env), pos = env)
}

# Print a either a string/value on one line after removing leading or
# trailing whitespace.
# Possible to add automatic formatting of maximum amount of characters per line
# This function can be modified to enable optional logging of console output
printPretty <- function(message, silent = FALSE, addNewline = TRUE) {
    # Trim whitespace
    trim <- function (x) gsub('^\\s+|\\s+$', '', strsplit(x, '\n')[[1]])
    if (!silent) {
        cat(trim(message), '\n')
    }
}

#funStart <- function(fun = TRUE) {
#
#}