# Monte Carlo computation of posterior probability distribution

Compute the posterior joint probability distribution of the variates
conditional on the given data, by means of Markov-chain Monte Carlo,
using the package **Nimble**.

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

`learn()` saves several files in a directory. By default this output
directory is a temporary directory within the one used by
[`base::tempdir()`](https://rdrr.io/r/base/tempfile.html), but an
alternative one can be chosen with the argument `outputdir =`. The
output directory contain several diagnostic files for the Monte Carlo
computation; in particular:

- `MCtraces.pdf`: shows several trace plots of the Monte Carlo sampling;
  the correspondin data are in the file `MCtraces.rds`.

- `plotsamples_learnt.pdf`, `plotquantiles_learnt.pdf`: show the
  marginal posterior distributions of each individual variate, together
  with their variability (as samples or quantiles).

- `log-1.log`, `log-2.log`, ... one for each parallel core; report the
  progress of each parallel Monte Carlo computation and notes about it.

- `rng_seed.rds`: the state of the pseudorandom seed (see
  [base::Random](https://rdrr.io/r/base/Random.html)) when `learn()` was
  called.

- `metadata.csv`: a copy of the metadata.

It is recommended that you give an explicit argument `outputdir =` and
save the directory with the files above for future reference. In
particular, the `MCtraces.pdf` plot and `MCtraces.rds` data can be
useful to report Monte Carlo convergence in any work of yours that used
**Prova**.

## Details

This function takes as main inputs a set of data and metadata, and
computes the full joint probability distribution for new data, including
its variability; see
[`Pr()`](https://pglpm.github.io/prova/reference/Pr.md). This
computation can also be interpreted as an estimation of the full joint
frequency distribution of the variates in the *whole population*, beyond
the sample data, together with its uncertainty.

The computation is "non-parametric": probability or frequency
distributions are not assumed to be Gaussian or of any other specific
shape; no "model" is assumed. The mathematical representation of the
space of joint frequency distributions follows Dunson & Bhattacharya
(2011).

The computation is done via Markov-chain Monte Carlo, using the package
[**Nimble**](https://cran.r-project.org/package=nimble). "Convergence"
of the Monte Carlo computation is assessed with methods described in
Vehtari & al. (2021) and Kwon & al. (2025).

This function creates a "learnt" object, typically saved in a
`learnt.rds` file, which is used in all subsequent probabilistic
computations. Other information about the computation is provided in
logs and plots, saved in a directory specified by the user.

See
[`vignette('intro')`](https://pglpm.github.io/prova/articles/intro.md)
for an introductory example.

## References

For the mathematical representation of the frequency space:

- Dunson, Bhattacharya (2011): *Nonparametric Bayes regression and
  classification through mixtures of product kernels*
  <doi:10.1093/acprof:oso/9780199694587.003.0005>.

- Ishwaran, Zarepour (2002): *Exact and approximate sum representations
  for the Dirichlet process* <doi:10.2307/3315951>.

About Bayesian inference under exchangeability ("population inference"):

- Lindley, Novick (1981): *The role of exchangeability in inference*,
  <doi:10.1214/aos/1176345331>.

- Bernardo, Smith (2000): *Bayesian Theory*. Wiley
  <doi:10.1002/9780470316870>. About nonparametrics:

- Müller et al. (2015): *Nonparametric Bayesian inference*. IMS
  <doi:10.1007/978-3-319-18968-0>.

- Hjort et al. (2010): *Bayesian Nonparametrics*. Cambridge University
  Press <doi:10.1017/CBO9780511802478>.

About Markov-chain Monte Carlo and "convergence":

- de Valpine, Paciorek, Turek, & al. (2026): *NIMBLE: MCMC, Particle
  Filtering, and Programmable Hierarchical Modeling*
  <doi:10.5281/zenodo.1211190>,
  <https://cran.r-project.org/package=nimble>.

- Kwon & al. (2025): *MCMC stopping rules in latent variable modelling*
  <doi:10.1111/bmsp.12357>.

- Vehtari & al. (2021): *Rank-normalization, folding, and localization:
  an improved R-hat for assessing convergence of MCMC*
  <doi:10.1214/20-BA1221>.

- Gilks & al. (1998): *Markov Chain Monte Carlo in Practice*. Chapman &
  Hall/CRC <doi:10.1201/b14835>.

- D. J. C. MacKay (2005): *Information Theory, Inference, and Learning
  Algorithms*. Cambridge University Press
  <https://www.inference.org.uk/itila/book.html>.

Mathematical background of **Prova**:

- Porta Mana
  <https://github.com/pglpm/prova/raw/main/development/manual/pglpm2024-bayes_nonparam.pdf>.

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
## Create dataset with 3 points of variate 'V' for demonstration:
dataset <- data.frame(V = rnorm(n = 3))

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
#> Learning from 3 datapoints, 1 variates.
#> 
#>  Saving output in directory
#>  /tmp/Rtmpw0NX5F/prova-V1_D3_S10_260624T104950_1a7836bbfa 
#> 
#> Starting Monte Carlo sampling of 10 samples by 1 chains
#> in a space of 191 (effectively 259) dimensions.
#> Using 1 cores: 10 samples per chain, max 1 chains per core.
#> Requested:   ESS 450   rel.MCSE 0.047.
#> Core logs are being saved in individual files.
#> C-compiling samplers appropriate to the variates (Nimble v1.4.2)
#> this can take tens of minutes. Please wait...
#>                                                                
#> Finished Monte Carlo sampling.
#> Highest number of Monte Carlo iterations across chains: 10.
#> Highest number of used mixture components: 2.
#> 
#> NOTE: 1 chains were stopped before reaching required precision
#> in order to meet the required time constraints.
#> 
#> Checking test data
#> (#1 #2 #3):
#> rel. CI error: 0.277 to 0.757
#> ESS: 10 to 10
#> needed thinning: 1 to 5.72
#> average: 0.135 to 0.867
#> quantile width: 0.216 to 3.27
#> 
#> Plotting final Monte Carlo traces and marginal samples...
#> Total computation time: 34 secs
#> Average preparation & finalization time: 32 secs.
#> Average Monte Carlo time per chain: 0.54 secs.
#> Max total memory used: approx 360MB.
#> Max memory used per core: approx 360MB.
#> Removing temporary output files.
#> 
#> Finished.
#> ********************************************************
#>  Output saved in directory
#> /tmp/Rtmpw0NX5F/prova-V1_D3_S10_260624T104950_1a7836bbfa
#> ********************************************************
#> Closing connections to cores.

## Check structure of `learnt` object:
str(learnt)
#> List of 6
#>  $ Rmean      : num [1, 1:64, 1:10] 2.92 -1.05 3.25 2.55 -2.54 ...
#>  $ Rsd        : num [1, 1:64, 1:10] 0.346 7.822 1.841 1.145 0.638 ...
#>  $ W          : num [1:64, 1:10] 1.11e-12 1.17e-09 1.76e-27 3.55e-262 7.41e-28 ...
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
#>   ..$ tscale           : num 0.507
#>   ..$ plotmin          : num -2.34
#>   ..$ plotmax          : num 1.99
#>   ..$ V1               : logi NA
#>   ..$ V2               : logi NA
#>  $ auxinfo    :List of 12
#>   ..$ nchains          : num 1
#>   ..$ npoints          : int 3
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
#>   ..$ maxusedcomponents: num 2
#>   ..$ nonfinitechains  : num 0
#>   ..$ stoppedchains    : num 1
#>   ..$ rel. CI error    : Named num [1:4] 0.384 0.277 0.395 0.757
#>   .. ..- attr(*, "names")= chr [1:4] "gmean" "1" "2" "3"
#>   ..$ ESS              : Named num [1:4] 10 10 10 10
#>   .. ..- attr(*, "names")= chr [1:4] "gmean" "1" "2" "3"
#>   ..$ needed thinning  : Named num [1:4] 1.48 1 1.56 5.72
#>   .. ..- attr(*, "names")= chr [1:4] "gmean" "1" "2" "3"
#>   ..$ average          : Named num [1:4] 0.194 0.178 0.135 0.867
#>   .. ..- attr(*, "names")= chr [1:4] "gmean" "1" "2" "3"
#>   ..$ quantile width   : Named num [1:4] 0.216 0.26 0.229 3.267
#>   .. ..- attr(*, "names")= chr [1:4] "gmean" "1" "2" "3"
```
