#' Create a grid of values for a variate
#'
#' This function create a set of values for a variate `vrt`, based on the metadata stored in a `learnt` object created by the [learn()] function). The set of values depends on the type of variate (nominal or continuous, rounded, and so on, see [`metadata`]). The range of values is chosen to include, and extend slightly beyond, the range observed in the data used in the [learn()] function. Variate domains are always respected.
#'
#' @param vrt character: name of the variate, must match one of the names in the `metadata` file provided to the [learn()] function.
#' @param learnt Either a string with the name of a directory or full path for a 'learnt.rds' object, produced by the [learn()] function, or such an object itself.
#' @param length.out numeric, positive (default 129): number of values to be created; used only for continuous, non-rounded variates (see [`metadata`]).
#'
#' @return A numeric or character vector of values.
#'
#' @export
vrtgrid <- function(
    vrt,
    learnt,
    length.out = 129
){

    ## Extract auxmetadata
    ## If learnt is a string, check if it's a folder name or file name
    if (is.character(learnt)) {
        ## Check if 'learnt' is a folder containing learnt.rds
        if (file_test('-d', learnt) &&
                file.exists(file.path(learnt, 'learnt.rds'))) {
            learnt <- readRDS(file.path(learnt, 'learnt.rds'))
        } else {
            ## Assume 'learnt' the full path of learnt.rds
            ## possibly without the file extension '.rds'
            learnt <- paste0(sub('.rds$', '', learnt), '.rds')
            if (file.exists(learnt)) {
                learnt <- readRDS(learnt)
            } else {
                stop('The argument "learnt" must be a folder containing learnt.rds, or the path to an rds-file containing the output from "learn()".')
            }
        }
    }
    ## Add check to see that learnt is correct type of object?
    auxmetadata <- learnt$auxmetadata
    rm(learnt)

    ## Consistency checks
    if(!(vrt %in% auxmetadata$name)){
        stop('Variate "', vrt, '" not present in the list of variates.')
    }

    adata <- as.list(auxmetadata[auxmetadata$name == vrt, ])

        if(adata$mcmctype %in% c('R', 'C')){
            out <- seq(adata$plotmin, adata$plotmax, length.out = length.out)
        } else if(adata$mcmctype == 'D'){
            out <- seq(adata$plotmin, adata$plotmax, by = 2 * adata$halfstep)
            ## step <- as.integer(ceiling((adata$plotmax - adata$plotmin) /
            ##                                (2 * adata$halfstep)))
            ## seqstep <- seq_len(step)
            ## factors <- seqstep[!(step %% seqstep)]
            ## step <- factors[which.min(abs(factors - length.out + 1))] *
            ##     2 * adata$halfstep
            ## out <- seq(adata$plotmin, adata$plotmax,
            ##     length.out = 1 + round((adata$plotmax - adata$plotmin)/step))
        } else if(adata$mcmctype %in% c('B', 'O', 'N')){
            out <- unname(unlist(adata[paste0('V', seq_len(adata$Nvalues))]))
        } else {
            stop('Unknown variate type for "', vrt, '".')
        }
    out
}
