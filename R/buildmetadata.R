#' Build preliminary metadata flie
#'
#' @param data data.frame object or filepath
#' @param file string: name of output metadata file;
#'   NULL: output metadata as data.frame
#' @param includevrt string or NULL: name variates in data to be included
#' @param excludevrt string or NULL: name variates in data to be excluded
#' @param addsummary2metadata Bool: also output some diagnostic statistics
#'    in the metadata?
#' @param backupfiles Bool: rename previous metadata file if it exists?
#' @param verbose Bool: output heuristics for each variate, default TRUE
#'
#' @return nothing or data.table object
#'
#' @export
buildmetadata <- function(
    data,
    file = NULL,
    includevrt = NULL,
    excludevrt = NULL,
    addsummary2metadata = FALSE,
    backupfiles = FALSE,
    verbose = TRUE
) {

    ## Greater Common Denominator function
    ## used to calculate multiplicity of values of continuous variates
    ## and guess whether they are rounded
    gcd <- function(...) {
        suppressWarnings(Reduce(function(a, b) {
            if (b == 0) a else Recall(b, a %% b)
        }, c(...)))
    }

    ## Read datafile if it exists
    datafile <- NULL
    if (is.character(data) && file.exists(data)) {
        datafile <- data
        data <- read.csv(datafile, na.strings = '')
    }
    data <- as.data.frame(data)

#### Prepare empty metadata frame
    metadata <- as.data.frame(
        c(list(name = NA,
            type = NA,
            Nvalues = NA,
            rounding = NA,
            domainmin = NA,
            domainmax = NA,
            minincluded = NA,
            maxincluded = NA,
            ## lowvalue = NA,
            ## centralvalue = NA,
            ## highvalue = NA,
            plotmin = NA,
            plotmax = NA),
            if (addsummary2metadata) {
                list(datamin = NA,
                    datamax = NA,
                    datamaxrep = NA,
                    dataNvalues = NA)
            }
        )
    )[-1,]

    cat('\n')
    variatelist <- colnames(data)
    if(!is.null(includevrt)) {
        variatelist <- intersect(variatelist, includevrt)
    }
    if(!is.null(excludevrt)) {
        variatelist <- setdiff(variatelist, excludevrt)
    }
    ## Loop over variates (columns) in data frame
    for (name in variatelist) {
        if(verbose){ cat('\n* Variate', paste0('"', name, '"'), ':\n') }
        ## remove missing values
        x <- data[[name]]
        x <- x[!is.na(x)]
        uniquex <- length(unique(x))
        ## if this variate has only one value, then it bears no information
        if (uniquex <= 1) {
            message('WARNING: variate "', name, '"',
                ' does not have at least two distinct values.',
                '\nDiscarded because non-informative\n')
            next
        }

        ## If the variate is numeric,
        ## calculate several characteristics used later
        if (is.numeric(x)) {
            ## these two are used for diagnostic purposes
            datamin <- min(x)
            datamax <- max(x)
            ## values strictly within the interior of the domain
            ix <- x[!(x %in% range(x))]
            maxrep <- max(table(ix)) # max of repeated inner values
            meanrep <- mean(table(ix)) # average of repeated inner values
            ## differences between consecutive unique values
            jumpsx <- unique(signif(diff(sort(unique(x))), 3))
            rangex <- diff(range(x))
            multi <- 10^(-min(floor(log10(jumpsx))))
            ## jumps between unique values are integer multiples
            ## of 'jumpquantum'
            jumpquantum <- gcd(round(jumpsx * multi)) / multi
        }

#### Now make educated guess for the type of variate

        if (uniquex == 2) {
            ## Binary variate? (has only two unique values)
            type <- 'binary'
            Nvalues <- 2
            rounding <- NA
            domainmin <- NA
            domainmax <- NA
            minincluded <- NA
            maxincluded <- NA
            ## lowvalue <- NA
            ## centralvalue <- NA
            ## highvalue <- NA
            plotmin <- NA
            plotmax <- NA
            datavalues <- sort(as.character(unique(x)))
            names(datavalues) <- paste0('V', 1:2)
            ##
            if(verbose){
                cat('  Two different values detected:\n',
                    paste0('"', datavalues, '"', collapse=', '),
                    '\n')
                cat('  Assuming variate to be BINARY.\n')
            }

        } else if (!is.numeric(x)) {
            ## Nominal variate? (non-numeric values)
            type <- 'nominal'
            Nvalues <- uniquex
            rounding <- NA
            domainmin <- NA
            domainmax <- NA
            minincluded <- NA
            maxincluded <- NA
            ## lowvalue <- NA
            ## centralvalue <- NA
            ## highvalue <- NA
            plotmin <- NA
            plotmax <- NA
            datavalues <- sort(as.character(unique(x)))
            names(datavalues) <- paste0('V', 1:Nvalues)
            ##
            if(verbose){
                cat('  - ', Nvalues, 'different non-numeric values detected:\n',
                    paste0('"', datavalues, '"', collapse=', '),
                    '\n')
                cat('  Assuming variate to be NOMINAL.\n')
            }

        } else if (uniquex <= 10) {
            ## Ordinal variate with few numeric values?
            type <- 'ordinal'
            Nvalues <- uniquex
            rounding <- NA
            domainmin <- NA
            domainmax <- NA
            minincluded <- NA
            maxincluded <- NA
            ## lowvalue <- NA
            ## centralvalue <- NA
            ## highvalue <- NA
            plotmin <- NA
            plotmax <- NA
            datavalues <- sort(as.character(unique(x)))
            names(datavalues) <- paste0('V', 1:Nvalues)
            ##
            if(verbose){
                cat('  - Only', Nvalues, 'different numeric values detected:\n',
                    paste0('"', datavalues, '"', collapse=', '),
                    '\n')
                cat('  Assuming variate to be ORDINAL.\n')
            }

        } else if (jumpquantum >= 1) {
            ## Ordinal variate with many numeric values?
            type <- 'ordinal'
            Nvalues <- uniquex
            rounding <- jumpquantum
            domainmin <- datamin
            domainmax <- datamax
            minincluded <- NA
            maxincluded <- NA
            ## lowvalue <- NA
            ## centralvalue <- NA
            ## highvalue <- NA
            plotmin <- NA
            plotmax <- NA
            datavalues <- NULL
            ##
            if(verbose){
                cat(' - ', Nvalues, 'different numeric values detected\n')
                cat('  there is a large minimum distance of', jumpquantum,
                    'between datapoints.\n')
                cat('  Assuming variate to be ORDINAL.\n')
            }
            message('WARNING: please check variate "', name, '":',
                '\nit seems ordinal, but it could also be',
                ' a rounded continous variate\n')


        } else {
            ## The variate seems continuous
            type <- 'continuous'
            Nvalues <- Inf
            ## preliminary values, possibly modified below
            rounding <- NA
            domainmin <- signif(datamax - 3 * rangex, 1)
            domainmax <- signif(datamax + 3 * rangex, 1)
            minincluded <- FALSE
            maxincluded <- FALSE

            Q1 <- quantile(x, probs = 0.25, type = 6)
            ## centralvalue <- quantile(x, probs = 0.5, type = 6)
            Q3 <- quantile(x, probs = 0.75, type = 6)
            ## ## Borderline case if the first and second quartile have the same value
            ## if (lowvalue == highvalue) {
            ##   lowvalue <- if (sum(x < Q1) > 0) {
            ##                 max(x[x < Q1])
            ##               } else {
            ##                 max(x[x <= Q1])
            ##               }
            ##   highvalue <- if (sum(x > Q3) > 0) {
            ##                  min(x[x > Q3])
            ##                } else {
            ##                  min(x[x >= Q3])
            ##                }
            ## }
            plotmin <- datamin - (Q3 - Q1) / 2
            plotmax <- datamax + (Q3 - Q1) / 2
            datavalues <- NULL
            ##
            if(verbose){
                cat('  - Numeric values between',
                    datamin, 'and', datamax, '\n')
                cat('  Assuming variate to be CONTINUOUS.\n')
            }

#### Consider subcases for openness of domain
            if (sum(x == min(x)) > meanrep) {
                ## seems to have left-closed domain
                domainmin <- datamin
                plotmin <- datamin
                minincluded <- TRUE
                ##
                if(verbose){
                    cat('  - Many datapoints have minimum value\n')
                    cat('  Assuming "domainmin" to be the minimum observed value\n')
                    cat('  and to be included in the domain',
                        '(singular probabilities there).\n')
                }
            } else if (all(x > 0)) {
                ## seems to be strictly positive
                domainmin <- 0
                plotmin <- max((domainmin + datamin) / 2, plotmin)
                ##
                if(verbose){
                    cat('  - All values are non-negative\n')
                    cat('  Assuming "domainmin" to be 0\n')
                    cat('  but value "0" not to be part of domain.\n')
                }
            }

            if (sum(x == datamax) > meanrep) {
                ## seems to right-closed domain
                domainmax <- datamax
                plotmax <- datamax
                maxincluded <- TRUE
                ##
                if(verbose){
                    cat('  - Many datapoints have maximum value\n')
                    cat('  Assuming "domainmax" to be the maximum observed value\n')
                    cat('  and to be included in the domain',
                        '(singular probabilities there).\n')
                }
            }

#### Consider subcases for rounding
            ## ## Variate has integer values, it's possibly ordinal
            ## if (jumpquantum >= 1) {
            ##     message('\nNOTE: variate ', name, ' might be ordinal,\n',
            ##         'but is treated as continuous and rounded\n',
            ##         'owing to its large range of values.\n')
            ## }

            if (!(
                jumpquantum / rangex < 1e-5 || # no visible gaps between values
                (jumpquantum / rangex < 1e-3 &&
                 maxrep <= 100)  # visible gaps, but not many value repetitions
            )) {
                ## The variate seems to have quantized differences between unique values
                ## hence it might be rounded or integer/ordinal

                ## if(rangex / rounding > 256) {
                ##   ## The variate seems rounded,
                ##   ## but no latent-variable representation is needed
                ##   message('\nNOTE: variate ', name, ' is probably rounded,\n',
                ##           'but is treated as not rounded ',
                ##           'owing to its large range of values.\n')
                ##   rounding <- 0
                ## } else {
                rounding <- jumpquantum
                ##
                if(verbose){
                    cat('  - There is a minimum distance of', jumpquantum,
                        'between datapoints.\n')
                    cat('  Assuming variate to be ROUNDED.\n')
                }
                ## domainmin <- signif(datamin - 4 * rangex, 1)
                ## domainmax <- signif(datamax + 4 * rangex, 1)
            } # End rounded case

        } # End continuous-variate case

#### Create metadata object
        metadata <- merge(metadata,
            c(
                list(name = name,
                    type = type,
                    Nvalues = Nvalues,
                    rounding = rounding,
                    domainmin = domainmin,
                    domainmax = domainmax,
                    minincluded = minincluded,
                    maxincluded = maxincluded,
                    ## lowvalue = signif(lowvalue,3),
                    ## centralvalue = signif(centralvalue,3),
                    ## highvalue = signif(highvalue,3),
                    plotmin = signif(plotmin,2),
                    plotmax = signif(plotmax,2)),
                if (addsummary2metadata) {
                    list(datamin = datamin,
                        datamax = datamax,
                        datamaxrep = max(table(x)),
                        dataNvalues = uniquex)
                },
                as.list(datavalues)
            ),
            sort = FALSE, all = TRUE)
    } # End loop over variates

                                        # Save to file if the file parameter is set
    if (!is.null(file)) {
                                        # If file is a string, name the metadata file that
        if (is.character(file)) {
            file <- paste0(sub('.csv$', '', file), '.csv')
        }
                                        # If file is not a string, and the data is read from a file,
                                        # save the metadata as metadata_<datafile>.csv
        else if (!is.null(datafile)) {
            file <- paste0(sub('.csv$', '', file), '_metadata.csv')
        }
                                        # If the file already exists, rename the old file as a backup file
                                        # and name this metadata file as intended
        if (backupfiles) {
            if (file.exists(file)) {
                file.rename(from = file, to = paste0(sub('.csv$', '', file), '_bak', format(Sys.time(), '%y%m%dT%H%M%S'),'.csv'))
            }
        }
                                        # Save the file
        write.csv(metadata, file, row.names = FALSE, quote = FALSE, na = '')
        cat('\nSaved proposal metadata file as', paste0('"', file, '"'), '\n')

    } else {
                                        # Else just print to console
        metadata
    }
}
