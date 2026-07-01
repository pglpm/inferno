#' Example metadata file
#'
#' @description
#' A [data frame][base::data.frame()] containing the prior information about the variates `species` and `bill_len` of the [datasets::penguins] dataset.
#'
#' @format ## `metadataExample`
#' A [data frame][base::data.frame()] with 2 rows and 10 columns.
#'
#' @returns No return value.
#'
#' @seealso
#' [metadatatemplate()] which helps producing this kind of metadata files from a given dataset.
#'
#' [learn()] which needs this kind of metadata files to "learn" from data.
"metadataExample"


#' Example `learnt` object produced by learn()
#'
#' @description
#' An example `learnt` object obtained by means of the [learn()] function, using the [datasets::penguins] dataset and the metadata in [metadataExample], according to the call
#'
#' ```
#' learn(data = penguins, metadata = metadataExample,
#'   nsamples = 225, nchains = 15)
#' ```
#'
#' It is a list that essentially contains posterior hyperparameters for drawing statistical inferences about the variates `species` and `bill_len`.
#'
#' **Note** that the [learn()] function that produced `learntExample` was called with the option to create only a limited number (225) of Monte Carlo samples, in order to reduce its memory size. Thus the numerical error associated with the Monte Carlo approximation is relatively in inferences drawn from the posterior hyperparameters saved in `learntExample`. It is only meant to be used for illustration purposes of the package's capabilities.
#'
#' @format ## `learntExample`
#' A list containing results from Markov-chain Monte Carlo computation, including diagnostics and variate metadata.
#'
#' @returns No return value.
#'
#' @seealso
#' [learn()], which produces this kind of object.
#'
#' [Pr()], [qPr()], [rPr()], [mutualinfo()]: functions that require this kind of object in order to calculate probabilities and quantiles, generate data points, and calculate mutual information.
"learntExample"
