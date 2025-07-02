#### Cleanup a learn()-output directory
#' @keywords internal
util_cleanup <- function(path){
    if(all(file.exists(file.path(path, c(
        'learnt.rds',
        'MCtraces.rds',
        'MCtraces.pdf',
        'log-1.log',
        'plotquantiles_learnt.pdf',
        'plotsamples_learnt.pdf'
    ))))){
        invisible(file.remove(dir(path,
            pattern = paste0('^___.*\\..*$'),
            full.names = TRUE
        )))
    } else {
        stop("This doesn't look like a 'learn()' output directory")
    }
}
