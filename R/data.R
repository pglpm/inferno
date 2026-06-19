#' Example metadata file
#'
#' @description
#' A data frame containing the prior information about the variates `species` and `bill_len` of the [datasets::penguins] dataset.
#'
#' @format ## `metadataExample`
#' A data frame with 2 rows and 10 columns.
#'
#' @seealso
#' [metadatatemplate()] which helps producing this metadata file from a given dataset, [learn()] which uses to produce hyperparameters for posterior infereces, [learntExample] as an (unrealistic) `learnt' object produced from these metadata.
"metadataExample"


#' Example `learnt` object produced by learn()
#'
#' @description
#' An example `learnt` object obtained by means of the [learn()] function, using the [datasets::penguins] dataset and the metadata in [metadataExample]. It is a list that essentially contains posterior hyperparameters for drawing statistical inferences about the variates `species` and `bill_len`.
#'
#' **Note** that the [learn()] function that produced `learntExample` was called with special, poor parameter values to reduce its computation time and memory size, at the cost of Monte Carlo convergence. Thus no real inferential or statistical meaning should be attached to probabilities obtained from `learntExample`. It is only used for illustration purposes of the package's capabilities.
#'
#' @format ## `learntExample`
#' A list containing results from Markov-chain Monte Carlo computation, including diagnostics and variate metadata.
#' 
#' @seealso
#' [learn()] which produces this kind of object, [Pr()] which calculates variousposterior probabilities based on this kind of object.
"learntExample"
