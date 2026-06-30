lll <- function(
    x = 0,
    sinkfile = 'testlog.log',
    verbose = TRUE
){
    cl <- parallel::makeCluster(2)
    closecoresonexit <- function(){
        message('Closing connections to cores.')
        parallel::stopCluster(cl)
        sink()
        gc(full = TRUE)
    }
    on.exit(closecoresonexit())

    sink(file = sinkfile, append=TRUE, split=verbose)
    cat('This is main cat 1.\n')
    message('This is message 1.\n')

    parallel::parLapply(
        cl = cl,
        X = 1:2,
        fun = lllfun,
        sinkfile = sinkfile,
        verbose = verbose
    )

    cat('This is maincat 2.\n')

    x
}

lllfun <- function(
    acore,
    sinkfile,
    verbose
){
    file2 <- paste0('testlog_', acore, '.log')

    sink()
    sink(file = sinkfile, append = TRUE, split = verbose)
    cat('\nThis is outcat 1, core', acore, '\n')

    sink(file = file2, append = TRUE, split = FALSE)
    cat('\nThis is intcat 1, core', acore, '\n')

    sink()
    sink(file = sinkfile, append = TRUE, split = verbose)
    cat('\nThis is outcat 2, core', acore, '\n')

    sink(file = file2, append = TRUE, split = FALSE)
    cat('\nThis is intcat 2, core', acore, '\n')

    closecons <- function(){
        ## Close output to log files
        ## flush(outcon)
        sink(file = NULL, type = 'output')
        sink(file = NULL, type = 'message')
        ## close(outcon)
    }
    on.exit(closecons())

    acore
}
