# Monte Carlo computation of posterior probability distribution

Compute the posterior probability distribution of the variates
conditional on the given data.

## Usage

``` r
learn(
  data,
  metadata,
  auxdata = NULL,
  outputdir = NULL,
  nsamples = 3600,
  nchains = 8,
  nsamplesperchain = 450,
  parallel = TRUE,
  seed = NULL,
  cleanup = TRUE,
  appendtimestamp = TRUE,
  appendinfo = TRUE,
  outputvalue = "directory",
  subsampledata = NULL,
  prior = missing(data) || is.null(data),
  startupMCiterations = 3600,
  minMCiterations = 0,
  maxMCiterations = +Inf,
  maxhours = +Inf,
  ncheckpoints = 12,
  maxrelMCSE = +Inf,
  minESS = 450,
  initES = 2,
  thinning = NULL,
  plottraces = !cleanup,
  showKtraces = FALSE,
  showAlphatraces = FALSE,
  hyperparams = list(ncomponents = 64, minalpha = -4, maxalpha = 4, byalpha = 1, Rshapelo
    = 0.5, Rshapehi = 0.5, Rvarm1 = 3^2, Cshapelo = 0.5, Cshapehi = 0.5, Cvarm1 = 3^2,
    Dshapelo = 0.5, Dshapehi = 0.5, Dvarm1 = 3^2, Bshapelo = 1, Bshapehi = 1, Dthreshold
    = 1, tscalefactor = 4.266, Oprior = "Hadamard", Nprior = "Hadamard", avoidzeroW =
    NULL, initmethod = "datacentre", Qerror = pnorm(c(-1, 1)))
)
```

## Arguments

- data:

  A dataset, given as a
  [`base::data.frame()`](https://rdrr.io/r/base/data.frame.html) or as a
  file path to a CSV file.

- metadata:

  A
  [`metadata`](https://pglpm.github.io/inferno/reference/metadatatemplate.md)
  object, given either as a data.frame object, or as a file pa to a CSV
  file.

- auxdata:

  A larger dataset, given as a base::data.frame() or as a file path to a
  CSV file. Such a dataset would be too many to use in the Monte Carlo
  sampling, but can be used to calculate hyperparameters.

- outputdir:

  Character: path to folder where the output should be saved. If `NULL`
  (default), a directory is created that has the same name as the data
  file but with suffix "`_output_`". If `FALSE`, a directory is created
  in the temporary-directory space.

- nsamples:

  Integer: number of desired Monte Carlo samples. Default 3600.

- nchains:

  Integer: number of Monte Carlo chains. Default 4.

- nsamplesperchain:

  Integer: number of Monte Carlo samples per chain.

- parallel:

  Logical or positive integer or cluster object. `TRUE`: use roughly
  half of available cores; `FALSE`: use serial computation; integer: use
  this many cores. It can also be a cluster object previously created
  with
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html);
  in this case the parallel computation will use this object.

- seed:

  Integer: use this seed for the random number generator. If missing or
  `NULL` (default), do not set the seed.

- cleanup:

  Logical: remove diagnostic files at the end of the computation?
  Default `TRUE`.

- appendtimestamp:

  Logical: append a timestamp to the name of the output directory
  `outputdir`? Default `TRUE`.

- appendinfo:

  Logical: append information about dataset and Monte Carlo parameters
  to the name of the output directory `outputdir`? Default `TRUE`.

- outputvalue:

  Character: if `'directory'`, return the output directory name as
  `VALUE`; if character `'learnt'`, return the `learnt` object
  containing the parameters obtained from the Monte Carlo computation.
  Any other value: `VALUE` is `NULL`.

- subsampledata:

  Integer: use only a subset of this many datapoints for the Monte Carlo
  computation.

- prior:

  Logical: Calculate the prior distribution?

- startupMCiterations:

  Integer: number of initial Monte Carlo iterations. Default 3600.

- minMCiterations:

  Integer: minimum number of Monte Carlo iterations to be doneby a
  chain. Default 0.

- maxMCiterations:

  Integer: Do at most this many Monte Carlo iterations per chain.
  Default `Inf`.

- maxhours:

  Numeric: approximate time limit, in hours, for the Monte Carlo
  computation to last. Default `Inf`.

- ncheckpoints:

  Integer: number of datapoints to use for checking when the Monte Carlo
  computation should end. If `NULL`, this is equal to number of
  variates + 2. If Inf, use all datapoints. Default 12.

- maxrelMCSE:

  Numeric positive: desired maximal *relative Monte Carlo Standard
  Error* of calculated probabilities with respect to their variability
  with new data. Default `+Inf`, so that `minESS` is used instead.
  `maxrelMCSE` is related to `minESS` by
  `maxrelMCSE = 1/sqrt(minESS + initES)`.

- minESS:

  Numeric positive: desired minimal Monte Carlo *Expected Sample Size*.
  If `NULL`, it is equal to the final `nsamplesperchain`. Default 400.
  `minESS` is related to `maxrelMCSE` by
  `minESS = 1/maxrelMCSE^2 - initES`.

- initES:

  Numeric positive: number of initial *Expected Samples* to discard.

- thinning:

  Integer: thin out the Monte Carlo samples by this value. If `NULL`
  (default): let the diagnostics decide the thinning value.

- plottraces:

  Logical: save plots of the Monte Carlo traces of diagnostic values?
  Default `TRUE`.

- showKtraces:

  Logical: save plots of the Monte Carlo traces of the K parameter?
  Default `FALSE`.

- showAlphatraces:

  Logical: save plots of the Monte Carlo traces of the Alpha parameter?
  Default `FALSE`.

- hyperparams:

  List: hyperparameters of the prior.

## Value

Name of directory containing output files, or learnt object, or `NULL`,
depending on argument `outputvalue`.

## Details

This function takes as main inputs a set of data and metadata, and
computes the probability distribution for new data. Its computation can
also be interpreted as an estimation of the frequencies of the variates
in the *whole population*, beyond the sample data. The probability
distribution is not assumed to be Gaussian or of any other specific
shape. The computation is done via Markov-chain Monte Carlo.

This function creates an object, contained in a `learnt.rds` file, which
is used in all subsequent probabilistic computations. Other information
about the computation is provided in logs and plots, saved in a
directory specified by the user.

See
[`vignette('inferno_start')`](https://pglpm.github.io/inferno/articles/inferno_start.md)
for an introductory example.

## Examples

``` r
## Create dataset with 10 points of variate 'V' for demonstration
dataset <- data.frame(V = rnorm(n = 10))

## Create metadatafile
metadata <- data.frame(
    name = 'V',
    type = 'continuous'
)
```
