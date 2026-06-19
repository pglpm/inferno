#' Example metadata file
#'
#' @description
#' A data frame containing the prior information about the variates `species` and `bill_len` of the [datasets::penguins] dataset.
#'
#' @format ## `metadataExample`
#' A data frame with 2 rows and 10 columns.
#'
#' @seealso
#' [metadatatemplate()] which helps producing this metadata file from a given dataset, [learn()] which uses to produce hyperparameters for posterior infereces, [learntExample] as a `learnt' object produced from these metadata.
"metadataExample"


#' Example `learnt` object produced by learn()
#'
#' @description
#' An example `learnt` object obtained by means of the [learn()] function, using the [datasets::penguins] dataset and the metadata in [metadataExample]. It is a list that essentially contains posterior hyperparameters for drawing statistical inferences about the variates `species` and `bill_len`.
#'
#' **Note** that the [learn()] function that produced `learntExample` was called with the option to create only a limited number (360) of Monte Carlo samples, in order to reduce its memory size. Thus the numerical error associated with the Monte Carlo approximation is relatively in inferences drawn from the posterior hyperparameters saved in `learntExample`. It is only meant to be used for illustration purposes of the package's capabilities.
#'
#' @format ## `learntExample`
#' A list containing results from Markov-chain Monte Carlo computation, including diagnostics and variate metadata.
#' 
#' @seealso
#' [learn()] which produces this kind of object, [Pr()] which calculates variousposterior probabilities based on this kind of object.
"learntExample"
