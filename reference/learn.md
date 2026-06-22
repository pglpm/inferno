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
  appendinfo = TRUE,
  valueislearnt = TRUE,
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

  [metadata](https://pglpm.github.io/prova/reference/metadatatemplate.md)
  about the dataset's variates, given either as a data frame or as a
  file path to a CSV file.

- auxdata:

  A larger dataset, given as a data frame or as a file path to a CSV
  file. Such a dataset would be too many to use in the Monte Carlo
  sampling, but can be used to calculate hyperparameters.

- outputdir:

  `NULL` or `NA` or character: path to folder where output information
  and diagnostics should be saved. If `NULL` (default), a directory is
  created in the temporary-directory space given by
  [`base::tempdir()`](https://rdrr.io/r/base/tempfile.html). If `NA`, a
  directory is created in the current working directory given by
  [`base::getwd()`](https://rdrr.io/r/base/getwd.html). If character,
  this is taken to be the output directory; it should of course be
  writable by the user.

- nsamples:

  Integer: number of desired Monte Carlo samples. Default 3600.

- nchains:

  Integer: number of Monte Carlo chains. Default 4.

- nsamplesperchain:

  Integer: number of Monte Carlo samples per chain.

- parallel:

  Logical or positive integer or cluster object. `TRUE` (default): use
  roughly half of available cores; `FALSE` (default): use serial
  computation; integer: use this many cores. It can also be a cluster
  object previously created with
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html);
  in this case the parallel computation will use this object.

- seed:

  Integer: use this seed for the random number generator. If missing or
  `NULL` (default), do not set the seed.

- cleanup:

  Logical: remove diagnostic files at the end of the computation?
  Default `TRUE`.

- appendinfo:

  Logical: append information about number of variates ('V'), number of
  data points ('D'), number of Monte Carlo samples ('S'), and timestamp,
  to the name of the output directory `outputdir`? The appended string
  has the format 'Vn_Dn_Sn_YYMMDDTHHMMSS'. Default `TRUE`.

- valueislearnt:

  Logical or `NULL`: should the `VALUE` returned be the `learnt` object
  containing the results from the Monte Carlo computation? If `FALSE`,
  then `VALUE` is the output directory name. If `NULL`, then `VALUE` is
  `NULL`. Default `TRUE`.

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

A "learnt" object, or name of directory containing output files, or
`NULL`, depending on argument `valueislearnt`.

## Details

This function takes as main inputs a set of data and metadata, and
computes the probability distribution for new data. Its computation can
also be interpreted as an estimation of the frequencies of the variates
in the *whole population*, beyond the sample data. The probability
distribution is not assumed to be Gaussian or of any other specific
shape. The computation is done via Markov-chain Monte Carlo.

This function creates an object, typically saved in a `learnt.rds` file,
which is used in all subsequent probabilistic computations. Other
information about the computation is provided in logs and plots, saved
in a directory specified by the user.

See
[`vignette('start')`](https://pglpm.github.io/prova/articles/start.md)
for an introductory example.

## See also

[`metadatatemplate()`](https://pglpm.github.io/prova/reference/metadatatemplate.md)
to help writing metadata files.

[`Pr()`](https://pglpm.github.io/prova/reference/Pr.md) to calculate
probabilities, and
[`qPr()`](https://pglpm.github.io/prova/reference/qPr.md) to calculate
quantiles, given the data processed by `learn()`.

[`rPr()`](https://pglpm.github.io/prova/reference/rPr.md) to generate
datapoints similar to the data processed by `learn()`.

[`mutualinfo()`](https://pglpm.github.io/prova/reference/mutualinfo.md)
to calculate mutual information given the data processed by `learn()`.

[`pread.csv()`](https://pglpm.github.io/prova/reference/prova.data.md)
and
[`pwrite.csv()`](https://pglpm.github.io/prova/reference/prova.data.md)
to read and write CSV files in the format used by `learn()`.

## Examples

``` r
## Create dataset with 5 points of variate 'V' for demonstration:
dataset <- data.frame(V = rnorm(n = 5))

## Create metadata file:
metadata <- data.frame(name = 'V', type = 'continuous')

## Learn from the data:
learnt <- learn(
  data = dataset, metadata = metadata,
      ## the following parameters are unrealistic
      ## only used to reduce computation time for this example
  nsamples = 10, nchains = 1, startupMCiterations = 10, maxhours = 0
)
#> 
#> Registered socket cluster with 1 nodes on host ‘localhost’.
#> 
#> Learning from 5 datapoints, 1 variates.
#> 
#>  Saving output in directory
#>  /tmp/RtmpXYdZSE/prova-V1_D5_S10_260622T090138_1a2869d4ceb7 
#> 
#> Starting Monte Carlo sampling of 10 samples by 1 chains
#> in a space of 191 (effectively 261) dimensions.
#> Using 1 cores: 10 samples per chain, max 1 chains per core.
#> Requested:   ESS 450   rel.MCSE 0.047.
#> Core logs are being saved in individual files.
#> C-compiling samplers appropriate to the variates (Nimble v1.4.2)
#> this can take tens of minutes. Please wait...
#>                                                                
#> Finished Monte Carlo sampling.
#> Highest number of Monte Carlo iterations across chains: 10.
#> Highest number of used mixture components: 3.
#> 
#> NOTE: 1 chains were stopped before reaching required precision
#> in order to meet the required time constraints.
#> 
#> Checking test data
#> (#1 #2 #3 #4 #5):
#> rel. CI error: 0.349 to 0.93
#> ESS: 4.73 to 10
#> needed thinning: 1.22 to 8.65
#> average: 0.0749 to 0.67
#> quantile width: 0.127 to 3.7
#> 
#> Plotting final Monte Carlo traces and marginal samples...
#> Total computation time: 41 secs
#> Average preparation & finalization time: 40 secs.
#> Average Monte Carlo time per chain: 0.7 secs.
#> Max total memory used: approx 360MB.
#> Max memory used per core: approx 360MB.
#> Removing temporary output files.
#> 
#> Finished.
#> **********************************************************
#>  Output saved in directory
#> /tmp/RtmpXYdZSE/prova-V1_D5_S10_260622T090138_1a2869d4ceb7
#> **********************************************************
#> Closing connections to cores.

## Check structure of `learnt` object:
str(learnt)
#> List of 6
#>  $ Rmean      : num [1, 1:64, 1:10] -3.34 2.52 2.83 -3.87 -2.56 ...
#>  $ Rsd        : num [1, 1:64, 1:10] 0.799 1.649 0.74 0.984 1.6 ...
#>  $ W          : num [1:64, 1:10] 6.31e-01 1.83e-11 4.20e-64 4.04e-285 1.36e-91 ...
#>  $ MCindex    : num [1:10(1d)] 1 2 3 4 5 6 7 8 9 10
#>  $ auxmetadata:'data.frame': 1 obs. of  24 variables:
#>   ..$ name             : chr "V"
#>   ..$ type             : chr "continuous"
#>   ..$ mcmctype         : chr "R"
#>   ..$ id               : int 1
#>   ..$ transform        : chr "identity"
#>   ..$ Nvalues          : int NA
#>   ..$ indexpos         : int NA
#>   ..$ halfstep         : num 0
#>   ..$ domainmin        : num -Inf
#>   ..$ domainmax        : num Inf
#>   ..$ minincluded      : logi FALSE
#>   ..$ maxincluded      : logi FALSE
#>   ..$ tdomainmin       : num -Inf
#>   ..$ tdomainmax       : num Inf
#>   ..$ domainminplushs  : num -Inf
#>   ..$ domainmaxminushs : num Inf
#>   ..$ tdomainminplushs : num -Inf
#>   ..$ tdomainmaxminushs: num Inf
#>   ..$ tlocation        : num 0.0796
#>   ..$ tscale           : num 0.459
#>   ..$ plotmin          : num -2.24
#>   ..$ plotmax          : num 2.01
#>   ..$ V1               : logi NA
#>   ..$ V2               : logi NA
#>  $ auxinfo    :List of 12
#>   ..$ nchains          : num 1
#>   ..$ npoints          : int 5
#>   ..$ hyperparams      :List of 21
#>   .. ..$ ncomponents : num 64
#>   .. ..$ minalpha    : num -4
#>   .. ..$ maxalpha    : num 4
#>   .. ..$ byalpha     : num 1
#>   .. ..$ Rshapelo    : num 0.5
#>   .. ..$ Rshapehi    : num 0.5
#>   .. ..$ Rvarm1      : num 9
#>   .. ..$ Cshapelo    : num 0.5
#>   .. ..$ Cshapehi    : num 0.5
#>   .. ..$ Cvarm1      : num 9
#>   .. ..$ Dshapelo    : num 0.5
#>   .. ..$ Dshapehi    : num 0.5
#>   .. ..$ Dvarm1      : num 9
#>   .. ..$ Bshapelo    : num 1
#>   .. ..$ Bshapehi    : num 1
#>   .. ..$ Dthreshold  : num 1
#>   .. ..$ tscalefactor: num 4.27
#>   .. ..$ Oprior      : chr "Hadamard"
#>   .. ..$ Nprior      : chr "Hadamard"
#>   .. ..$ initmethod  : chr "datacentre"
#>   .. ..$ Qerror      : num [1:2] 0.159 0.841
#>   ..$ maxiterations    : num 10
#>   ..$ maxusedcomponents: num 3
#>   ..$ nonfinitechains  : num 0
#>   ..$ stoppedchains    : num 1
#>   ..$ rel. CI error    : Named num [1:6] 0.397 0.462 0.902 0.628 0.93 ...
#>   .. ..- attr(*, "names")= chr [1:6] "gmean" "1" "2" "3" ...
#>   ..$ ESS              : Named num [1:6] 10 10 10 10 4.73 ...
#>   .. ..- attr(*, "names")= chr [1:6] "gmean" "1" "2" "3" ...
#>   ..$ needed thinning  : Named num [1:6] 1.58 2.13 8.14 3.94 8.65 ...
#>   .. ..- attr(*, "names")= chr [1:6] "gmean" "1" "2" "3" ...
#>   ..$ average          : Named num [1:6] 0.1262 0.0749 0.2593 0.1149 0.6699 ...
#>   .. ..- attr(*, "names")= chr [1:6] "gmean" "1" "2" "3" ...
#>   ..$ quantile width   : Named num [1:6] 0.127 0.129 1.106 0.239 3.701 ...
#>   .. ..- attr(*, "names")= chr [1:6] "gmean" "1" "2" "3" ...
```
