#' Metadata and helper function for metadata
#'
#' Metadata and helper function to build a template metadata file or object.
#'
#' The \code{\link{learn}} function needs metadata about the variates present in the data. Such metadata can be provided either as a `csv` file or as a \code{\link[base]{data.frame}}. The function `buildmetadata` creates a template metadata csv-file, or outputs a metadata data.frame, by trying to guess metadata information from the dataset. The user *must* then modify and correct this template, using it as a starting point to prepare the correct metadata information.
#'
#' @param data A dataset, given as a \code{\link[base]{data.frame}}
#' or as a file path to a csv file.
#' @param file String: name of csv file where the metadata should be saved;
#'   if `NULL`: output metadata as `VALUE`.
#' @param includevrt Character or `NULL`: name of variates in dataset to be included.
#' @param excludevrt Character or `NULL`: name of variates in dataset to be excluded.
#' @param addsummary2metadata Logical: also output some diagnostic statistics
#'    in the metadata? Default `FALSE`.
#' @param backupfiles Logical: rename previous metadata file if it exists?
#' Default `TRUE`.
#' @param verbose Logical: output heuristics for each variate? Default `TRUE`.
#'
#' @return If `file = NULL`, a preliminary metadata file is created
#'   and `VALUE` is `NULL`;
#'   otherwise `VALUE` is a \code{\link[base]{data.frame}} containing the metadata.
#'
#' @section Metadata information and format:
#'
#' In order to correctly learn from a dataset, the \code{\link{learn}} function needs information that is not contained in the data themeselves; that is, it needs *meta*data. Metadata are provided either as a `csv` file or as a \code{\link[base]{data.frame}}.
#'
#' A metadata file or data.frame must contain one row for each simple variate in the given inference problem, and the following fields (columns), even if some of them may be empty:
#'
#' `name`, `type`, `domainmin`, `domainmax`, `datastep`, `minincluded`, `maxinluded`, `V1`, `V2`, ...
#'
#' with possibly additional `V`-fields, sequentially numbered, and possibly the two fields `plotmin` and `plotmax`.
#'
#' The `type` field has three possible values: `nominal`, `ordinal`, `continuous`. The remaining fields that must be filled in depend on the `type` field. Here is a list of requirements:
#'
#' - **`nominal`** and **`ordinal`**: require *either* `V1`, `V2`, ... fields *or* `domainmin`, `domainmax`, `datastep` (all three) fields. No other fields are required, but `plotmin` and `plotmax` can optionally be specified.
#'
#' - **`continuous`**: requires `domainmin`, `domainmax`, `datastep`, `minincluded`, `maxincluded`, and optionally `plotmin` and `plotmax`.
#'
#' Here are the meanings and possible values of the fields:
#'
#' **`name`**: The name of the variate. This must be the same character string as it appears in the dataset (be careful about upper- and lower-case).
#'
#' **`type`**: The data type of variate `name`. Possible values are `nominal`, `ordinal`, `continuous`.
#'
#' - A *nominal* (also called *categorical*) variate has a discrete, finite number of possible values which have no intrinsic ordering. Examples could be a variate related to colour, with values "red", "green", "blue", and so on; or a variate related to cat breeds, with values "Siamese", "Abyssinian", "Persian", and so on. The possible valuesof the variate must be given in the fields `V1`, `V2`, and so on. It is important to include values that are possible but are *not* present in the dataset. A variate having only two possible values (binary variate), for example "yes" and "no", can be specified as nominal.
#'
#' - An *ordinal* variate has a discrete, finite number of possible values which do have an intrinsic ordering. Examples could be a Likert-scaled variate for the results of a survey, with values "very dissatisfied", "dissatisfied", "satisfied", "very satisfied"; or a variate related to the levels of some quantities, with values "low", "medium", "high"; or a variate having a numeric scale with values from 1 to 10. Whether a variate is nominal or ordinal often depends on the context. The possible values of the variate but be given in either one (but not both) or two ways: (1) in the fields `V1`, `V2`, ..., as for nominal variates; (2) as the fields `domainmin`, `domainmax`, `datastep`. Option (2) only works with numeric, equally spaced values: it assumes that the first value is `domainmin`, the second is `domainmin`+`datastep`, the third is `domainmin`+2*`datastep`, and so on up to the last value, `domainmax`.
#'
#' - A *continuous* variate has a continuum of values with an intrinsic ordering. Examples could be a variate related to the width of an object; or to the age of a person; or one coordinate of an object in a particular reference system. A continuous variate requires specification of the fields `domainmin`, `domainmax`, `datastep`, `minincluded`, `maxincluded`, `plotmin`, `plotmax`. Some naturally continuous variates are often rounded to a given precision; for instance, the age of a person might be reported as rounded to the nearest year (25 years, 26 years, and so on); or the length of an object might be reported to the nearest centimetre (1 m, 1.01 m, 1.02 m, and so on). The minimum distance between such rounded values **must** be reported in the `datastep` field; this would be `1` in the age example and `0.01` in the length example above. See below for further explanation of why reporting such rounding is important.
#'
#' **`domainmin`**: The minimum value that the variate (ordinal or continuous) can take on. Possible values are a real number, or `-Inf`, or an empty value, which is then interpreted as `-Inf`. Some continuous variates, like age or distance or temperature, are naturally positive, and therefore have `domainmin` equal `0`. But in other contexts the minimum value could be different. For instance, if a given inference problem only involves people of age 18 or more, then `domainmin` would be set to `18`.
#'
#' **`domainmax`**: The maximum value that the variate (ordinal or continuous) can take on. Possible values are a real number, or `+Inf`, or an empty value, which is then interpreted as `+Inf`. As with `domainmin`, the maximum value depends on the context. An age-related variate could theoretically have `domainmax` equal to `+Inf`; but if a given study categorizes some people as "90 years old or older", then `domainmax` should be set to `90`.
#'
#' **`datastep`**: The minimum distance between the values of a variate (ordinal or continuous). Possible values are a positive real number, or `0`, or an empty value, which is then interpreted as `0`. For a numeric ordinal variate, `datastep` is the step between consecutive values. For a continuous *rounded* variate, `datastep` is the minimum distance between different values that occurs because of rounding; see the examples given above. The function `buildmetadata` has some heuristics to determine whether the variate is rounded or not. See further details under the section Rounding below.
#'
#' **`minincluded`**, **`maxincluded`**: Whether the minimum (`domainmin`) and maximum(`domainmax`) values of a *continuous* variate can really appear in the data or not. Possible values are `TRUE`, `FALSE`, or an empty value, which is then interpreted as `FALSE`. Here are some examples about the meaning of these fields. (a) A continuous *unrounded* variate such as temperature has 0 as a minimum possible value `domainmin`, but this value itself is physically impossible and can never appear in data; in this case `minincluded` is set to `FALSE`. (b) A variate related to the *unrounded* length, in metres, of some objects may take on any positive real value; but suppose that all objects of length 5 or less are grouped together under the value `5`. It is then possible for several datapoints to have value `5`: one such datapoint could originally have the value 3.782341...; another the value 4.929673..., and so on. In this case `domainmin` is set to `5`, and `minincluded` is set to `TRUE`. Similarly for the maximum value of a variate and `maxincluded`. Note that if `domainmin` is `-Inf`, then `minincluded` is automatically set to `FALSE`, and similarly for `maxincluded` if `domainmax` is `+Inf`.
#'
#' **`plotmin`**, **`plotmax`**: The software draws some probability plots for each variate, after learning from the data with the \code{\link{learn}} function. `plotmin` and `plotmax` give the plotting range in these plots. They might be necessary because, for instance, although the minimum variate value is 0 and the maximum is 90, the variability of the data and the range of interest in the inference problem is between 20 and 40; these would then be the values of `plotmin` and `plotmax`. Possible values are real numbers, with `plotmin` strictly less than `plotmax`; one or both may be empty values. In case of an empty value, the software internally tries to choose an optimal value, typically so as to include all data given in the plot.
#'
#' @section Rounded continuous variates:
#'
#' To be written.
#'
#' @section Necessity of metadata:
#'
#' To be written.
#'
#' @aliases metadata buildmetadata
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

    ## keywords to recognize ordinal variates
    ordinalkeywords <- c('low', 'high', 'very', 'medium', 'strong',
        'poor', 'good', 'always')

    ## collect warning messages to be given at the end
    warninglist <- NULL

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
            domainmin = NA,
            domainmax = NA,
            datastep = NA,
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

    ## build list of variates, including only those in 'includevrt'
    ## and excluding those in 'excludevrt'
    variatelist <- colnames(data)
    if(!is.null(includevrt)) {
        variatelist <- intersect(variatelist, includevrt)
    }
    if(!is.null(excludevrt)) {
        variatelist <- setdiff(variatelist, excludevrt)
    }

    if(verbose){ cat('\nAnalyzing', length(variatelist), 'variates for',
        nrow(data), 'datapoints.\n') }

    ## Loop over variates (columns) in data frame
    ## cat('\n')
    for (name in variatelist) {
        if(verbose){ cat(paste0('\n* "', name, '" variate:\n')) }
        ## remove missing values
        x <- data[[name]]
        x <- x[!is.na(x)]

        ## check if it's a 'factor' object and transform
        if(is.factor(x)){
            x <- as.character(x)
        }

        unx <- sort(unique(x))
        uniquex <- length(unx)

        ## if this variate has only one value, then it bears no information
        if (uniquex <= 1) {
            if(verbose){
                cat('  Only one value present.\n')
            }
            warninglist <- c(warninglist,
                paste0('\n* "', name, '" variate',
                    ' does not have at least two distinct values.',
                    '\nDiscarded because non-informative.',
                    '\nPlease insert its characteristics by hand in the metadata file.'))
            next
        }

        ## If the values are strings,
        ## check whether they can be reinterpreted as numeric
        if(!is.numeric(unx) &&
               !any(is.na(suppressWarnings(as.numeric(unx))))) {
            x <- as.numeric(x)
            unx <- sort(as.numeric(unx))
        }

        ## If the variate is numeric,
        ## calculate several characteristics used later
        if (is.numeric(unx) && uniquex > 2) {
            ## these two are used for diagnostic purposes
            datamin <- min(x)
            datamax <- max(x)
            ## values strictly within the interior of the domain
            ix <- x[!(x %in% range(x))]
            maxrep <- max(table(ix)) # max of repeated inner values
            meanrep <- mean(table(ix)) # average of repeated inner values
            ## differences between consecutive unique values
            jumpsx <- unique(signif(diff(unx), 3))
            rangex <- diff(range(x))
            multi <- 10^(-min(floor(log10(jumpsx))))
            ## jumps between unique values are integer multiples
            ## of 'jumpquantum'
            jumpquantum <- gcd(jumpsx * multi) / multi
            datavalues <- NULL
        } else {
            datavalues <- as.character(unx)
        }

### Check whether the variate has been first rounded
### and then transformed to logarithmic scale.
### This may lead to problems
        if(is.numeric(unx) && datamin > 0 && uniquex > 2) {
            for(base in c(2, exp(1), 10)) {
                ljumpsx <- unique(signif(diff(log(unx, base)), 3))
                lrangex <- diff(log(range(x), base))
                lmulti <- 10^(-min(floor(log10(ljumpsx))))
                ## jumps between unique values are integer multiples
                ## of 'jumpquantum'
                ljumpquantum <- gcd(ljumpsx * lmulti) / lmulti
                if (!(
                    ljumpquantum / lrangex < 1e-5 || # no visible gaps between values
                        (ljumpquantum / lrangex < 1e-3 &&
                             maxrep <= 100)  # visible gaps, but not many value repetitions
                )) {
                    warninglist <- c(warninglist,
                        paste0('\n* "', name, '" variate',
                            ' appears to have been rounded\n',
                            'and then transformed to logarithmic scale.\n',
                            'This may lead to problems in the inference.\n',
                            'Preferably, transform it back to non-logarithmic scale.'))
                    break
                }
            }
        }

#### Now make educated guess for the type of variate

        ## if (uniquex == 2) {
        ##     ## Binary variate? (has only two unique values)
        ##     type <- 'binary'
        ##     Nvalues <- 2
        ##     datastep <- NA
        ##     domainmin <- NA
        ##     domainmax <- NA
        ##     minincluded <- NA
        ##     maxincluded <- NA
        ##     ## lowvalue <- NA
        ##     ## centralvalue <- NA
        ##     ## highvalue <- NA
        ##     plotmin <- NA
        ##     plotmax <- NA
        ##     names(datavalues) <- paste0('V', 1:2)
        ##     ##
        ##     if(verbose){
        ##         cat('  Two different values detected:\n',
        ##             paste0('"', datavalues, '"', collapse=', '),
        ##             '\n')
        ##         cat('  Assuming variate to be BINARY.\n')
        ##     }
        ##
        ## } else
            if (uniquex == 2 || (!is.numeric(x) &&
                       !any(sapply(ordinalkeywords, function(keyw){
                           grepl(keyw, datavalues, fixed = TRUE)
                       })))) {
            ## Nominal variate? (non-numeric values)
            type <- 'nominal'
            Nvalues <- uniquex
            domainmin <- NA
            domainmax <- NA
            datastep <- NA
            minincluded <- NA
            maxincluded <- NA
            ## lowvalue <- NA
            ## centralvalue <- NA
            ## highvalue <- NA
            plotmin <- NA
            plotmax <- NA
            names(datavalues) <- paste0('V', 1:Nvalues)
            ##
            if(verbose){
                cat('  - ', Nvalues, 'different',
                    if(uniquex > 2){' non-numeric'},
                    ' values detected:\n',
                    paste0('"', datavalues, '"', collapse=', '),
                    '\n')
                cat('  which do not seem to refer to an ordered scale.\n')
                cat('  Assuming variate to be NOMINAL.\n')
            }

        } else if (!is.numeric(x)) {
            ## Nominal variate? (non-numeric values)
            ## Ordinal variate with few numeric values?
            type <- 'ordinal'
            Nvalues <- uniquex
            domainmin <- NA
            domainmax <- NA
            datastep <- NA
            minincluded <- NA
            maxincluded <- NA
            ## lowvalue <- NA
            ## centralvalue <- NA
            ## highvalue <- NA
            plotmin <- NA
            plotmax <- NA
            names(datavalues) <- paste0('V', 1:Nvalues)
            ##
            if(verbose){
                cat('  - ', Nvalues, 'different non-numeric values detected:\n',
                    paste0('"', datavalues, '"', collapse=', '),
                    '\n')
                cat('  which seem to refer to an ordered scale.\n')
                cat('  Assuming variate to be ORDINAL.\n')
                cat('  Please appropriately reorder its values in metadata file.\n')
            }
            warninglist <- c(warninglist,
                paste0('\n* "', name, '" variate',
                    ' appears to be ordinal.',
                    '\nPlease appropriately reorder its values in metadata file.'))


#### Now we know that the variate is numeric

        } else if (uniquex <= 10 &&
                       (all(x == round(x)) || jumpquantum == round(jumpquantum))) {
            ## Ordinal variate with few numeric values?
            type <- 'ordinal'
            Nvalues <- uniquex
            ## if the values are spaced by an integer,
            ## no need to use the 'V' fields
            if(jumpquantum == round(jumpquantum)) {
                datastep <- jumpquantum
                domainmin <- datamin
                domainmax <- datamax
                datavalues <- NULL
            } else {
                domainmin <- NA
                domainmax <- NA
                datastep <- NA
                datavalues <- as.character(unx)
                names(datavalues) <- paste0('V', 1:Nvalues)
            }
            minincluded <- NA
            maxincluded <- NA
            ## lowvalue <- NA
            ## centralvalue <- NA
            ## highvalue <- NA
            plotmin <- NA
            plotmax <- NA
            ##
            if(verbose){
                cat('  - Only', Nvalues, 'different numeric values detected:\n')
                if(jumpquantum == round(jumpquantum)) {
                    cat('from', domainmin, 'to', domainmax,
                        'in steps of', jumpquantum, '\n')
                } else {
                    cat(paste0('"', datavalues, '"', collapse=', '), '\n')
                }
                cat('  Assuming variate to be ORDINAL.\n')
            }

            ## } else if (uniquex <= 10 && ) {
            ##     ## Ordinal variate with many numeric values?
            ##     type <- 'ordinal'
            ##     Nvalues <- uniquex
            ##     datastep <- jumpquantum
            ##     domainmin <- datamin
            ##     domainmax <- datamax
            ##     minincluded <- NA
            ##     maxincluded <- NA
            ##     ## lowvalue <- NA
            ##     ## centralvalue <- NA
            ##     ## highvalue <- NA
            ##     plotmin <- NA
            ##     plotmax <- NA
            ##     datavalues <- NULL
            ##     ##
            ##     if(verbose){
            ##         cat(' - ', Nvalues, 'different numeric values detected\n')
            ##         cat('  distance between datapoints is a multiple of integer',
            ##             jumpquantum, '\n')
            ##         cat('  Assuming variate to be ORDINAL.\n')
            ##     }
            ##     warninglist <- c(warninglist,
            ##         paste0('* "', name, '" variate',
            ##             ' appears to be ordinal,',
            ##             '\nbut it could also be a rounded continous variate'))


        } else {
            ## The variate seems continuous
            type <- 'continuous'
            Nvalues <- Inf
            ## preliminary values, possibly modified below
            datastep <- NA
            domainmin <- -Inf # signif(datamax - 3 * rangex, 1)
            domainmax <- +Inf # signif(datamax + 3 * rangex, 1)

            ## Q1 <- quantile(x, probs = 0.25, type = 6)
            ## centralvalue <- quantile(x, probs = 0.5, type = 6)
            ## Q3 <- quantile(x, probs = 0.75, type = 6)
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
            plotmin <- datamin - IQR(x, type=6) / 2
            plotmax <- datamax + IQR(x, type=6) / 2
            datavalues <- NULL
            ##
            if(verbose){
                cat('  - Numeric values between',
                    datamin, 'and', datamax, '\n')
                cat('  Assuming variate to be CONTINUOUS.\n')
            }

            roundedflag <- FALSE
#### Consider subcases for datastep
            ## for rounded variates it doesn't matter whether the boundaries
            ## are included in the domain or not
            if (!(
                jumpquantum / rangex < 1e-5 || # no visible gaps between values
                    (jumpquantum / rangex < 1e-3 &&
                         maxrep <= 100)  # visible gaps, but not many value repetitions
            )) {
                roundedflag <- TRUE
                ## The variate seems to have quantized differences
                ## between unique values
                ## hence it might be rounded or integer/ordinal

                ## if(rangex / datastep > 256) {
                ##   ## The variate seems rounded,
                ##   ## but no latent-variable representation is needed
                ##   message('\nNOTE: variate ', name, ' is probably rounded,\n',
                ##           'but is treated as not rounded ',
                ##           'owing to its large range of values.\n')
                ##   datastep <- 0
                ## } else {
                datastep <- jumpquantum
                minincluded <- maxincluded <- NA
                ##
                if(verbose){
                    cat('  - Distance between datapoints',
                        'is a multiple of', jumpquantum, '\n')
                    cat('  Assuming variate to be ROUNDED.\n')
                }
                if(jumpquantum >=1) {
                    warninglist <- c(warninglist,
                        paste0('\n* "', name, '" variate',
                            ' appears to be continuous and rounded,',
                            '\nbut it could also be an ordinal variate'))
                }
                ## domainmin <- signif(datamin - 4 * rangex, 1)
                ## domainmax <- signif(datamax + 4 * rangex, 1)

            } # End rounded case


#### Consider subcases for closedness of domain
            if(!roundedflag){
                minincluded <- FALSE
                maxincluded <- FALSE
            }
            if (sum(x == datamin) > 1) {
                ## seems to have left-closed domain
                domainmin <- datamin
                plotmin <- datamin
                if(!roundedflag){
                    minincluded <- TRUE
                }
                ##
                if(verbose){
                    cat('  - Several datapoints have minimum value', datamin, '\n')
                    cat('  Assuming "domainmin" to be this minimum observed value\n')
                    if(!roundedflag){
                        cat('  and to be included in the domain',
                            '(singular probabilities there).\n')
                    }
                }
            } else if (all(unx > 0)) {
                ## seems to be strictly positive
                domainmin <- 0
                plotmin <- max((domainmin + datamin) / 2, plotmin)
                if(!roundedflag){
                    minincluded <- FALSE
                }
                ##
                if(verbose){
                    cat('  - All values are positive\n')
                    cat('  Assuming "domainmin" to be 0\n')
                    if(!roundedflag){
                        cat('  with 0 excluded from domain.\n')
                    }
                }
            } else if (all(unx >= 0)) {
                ## seems to be non-negative
                domainmin <- 0
                plotmin <- max((domainmin + datamin) / 2, plotmin)
                if(!roundedflag){
                    minincluded <- TRUE
                }
                ##
                if(verbose){
                    cat('  - All values are non-negative\n')
                    cat('  Assuming "domainmin" to be 0\n')
                    if(!roundedflag){
                        cat('  with 0 included in the domain.\n')
                    }
                }
            }

            if (sum(x == datamax) > meanrep) {
                ## seems to right-closed domain
                domainmax <- datamax
                plotmax <- datamax
                if(!roundedflag){
                    maxincluded <- TRUE
                }
                ##
                if(verbose){
                    cat('  - Many datapoints have maximum value,', datamax, '\n')
                    cat('  Assuming "domainmax" to be this maximum observed value\n')
                    if(!roundedflag){
                        cat('  and to be included in the domain',
                            '(singular probabilities there).\n')
                    }
                }
            }

        } # End continuous-variate case

#### Create metadata object
        metadata <- merge(metadata,
            c(
                list(name = name,
                    type = type,
                    ## Nvalues = Nvalues,
                    domainmin = domainmin,
                    domainmax = domainmax,
                    datastep = datastep,
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
    if(verbose){cat('\n')}

    ## Print warnings
    if(!is.null(warninglist)){
        message(
            '=========',
            '\nWARNINGS',
            ' - please make sure to check these variates in the metadata file:')
        for(awarning in warninglist){message(awarning)}
        cat('=========\n')
    }

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
