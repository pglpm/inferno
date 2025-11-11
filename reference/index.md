# Package index

## Main functions

- [`flexiplot()`](https://pglpm.github.io/inferno/reference/flexiplot.md)
  : Plot numeric or character values

- [`hist(`*`<probability>`*`)`](https://pglpm.github.io/inferno/reference/hist.probability.md)
  : Plot the variability of an object of class "probability" as a
  histogram

- [`write.csvi()`](https://pglpm.github.io/inferno/reference/inferno.data.md)
  [`read.csvi()`](https://pglpm.github.io/inferno/reference/inferno.data.md)
  :

  Write and read data in **inferno**

- [`learn()`](https://pglpm.github.io/inferno/reference/learn.md) :
  Monte Carlo computation of posterior probability distribution

- [`metadatatemplate()`](https://pglpm.github.io/inferno/reference/metadatatemplate.md)
  : Metadata and helper function for metadata

- [`mutualinfo()`](https://pglpm.github.io/inferno/reference/mutualinfo.md)
  : Calculate mutual information between groups of joint variates

- [`plot(`*`<probability>`*`)`](https://pglpm.github.io/inferno/reference/plot.probability.md)
  : Plot an object of class "probability"

- [`plotFsamples()`](https://pglpm.github.io/inferno/reference/plotFsamples.md)
  : Plot one-dimensional posterior probabilities

- [`plotquantiles()`](https://pglpm.github.io/inferno/reference/plotquantiles.md)
  : Plot pairs of quantiles

- [`Pr()`](https://pglpm.github.io/inferno/reference/Pr.md) : Calculate
  posterior probabilities

- [`qPr()`](https://pglpm.github.io/inferno/reference/qPr.md) :
  Calculate quantiles

- [`rPr()`](https://pglpm.github.io/inferno/reference/rPr.md) : Generate
  datapoints

- [`subset(`*`<probability>`*`)`](https://pglpm.github.io/inferno/reference/subset.probability.md)
  :

  Subset variates of an object of class `probability`

- [`vrtgrid()`](https://pglpm.github.io/inferno/reference/vrtgrid.md) :
  Create a grid of values for a variate

## Internal functions

For developers.

- [`buildauxmetadata()`](https://pglpm.github.io/inferno/reference/buildauxmetadata.md)
  : Build preliminary metadata flie
- [`createQfunction()`](https://pglpm.github.io/inferno/reference/createQfunction.md)
  : Calculate and save transformation function for ordinal variates
- [`fftNGS()`](https://pglpm.github.io/inferno/reference/fftNGS.md) :
  Find optimal FFT size
- [`funAC()`](https://pglpm.github.io/inferno/reference/funAC.md) :
  Compute autocovariance
- [`funESS3()`](https://pglpm.github.io/inferno/reference/funESS3.md) :
  Compute ESS
- [`funESSLD()`](https://pglpm.github.io/inferno/reference/funESSLD.md)
  : Calculate MC effective sample size using LaplacesDemon's algorithm
- [`funMCEI()`](https://pglpm.github.io/inferno/reference/funMCEI.md) :
  Calculate quantile width through batches
- [`funMCEQ()`](https://pglpm.github.io/inferno/reference/funMCEQ.md) :
  Calculate credibility quantiles on estimated quantile
- [`funMCSELD()`](https://pglpm.github.io/inferno/reference/funMCSELD.md)
  : Calculate MC standard error using LaplacesDemon's batch means
- [`inferno`](https://pglpm.github.io/inferno/reference/inferno-package.md)
  [`inferno-package`](https://pglpm.github.io/inferno/reference/inferno-package.md)
  : inferno: Inference in R with Bayesian nonparametrics
- [`learnbind()`](https://pglpm.github.io/inferno/reference/learnbind.md)
  : Bind 3D arrays by first dimension
- [`mcjoin()`](https://pglpm.github.io/inferno/reference/mcjoin.md) :
  Concatenate mcsample objects
- [`mcsubset()`](https://pglpm.github.io/inferno/reference/mcsubset.md)
  : Eliminate samples from a 'learnt' object
- [`rowcumsum()`](https://pglpm.github.io/inferno/reference/rowcumsum.md)
  : Cumulative sum along first dimension
- [`rowinvcumsum()`](https://pglpm.github.io/inferno/reference/rowinvcumsum.md)
  : Inverse cumulative sum along first dimension
- [`util_Pcheckpoints()`](https://pglpm.github.io/inferno/reference/util_Pcheckpoints.md)
  : Calculate joint frequencies for checkpoints in learn()
- [`util_combineYX()`](https://pglpm.github.io/inferno/reference/util_combineYX.md)
  : Calculate probabilities, quantiles, etc, for all Y and X
  combinations
- [`util_denorm()`](https://pglpm.github.io/inferno/reference/util_denorm.md)
  : Utility function to avoid finite-precision accuracys
- [`util_joinPtraces()`](https://pglpm.github.io/inferno/reference/util_joinPtraces.md)
  : Join '\_\_\_\_tempPtraces-' files
- [`util_learntvar2sd()`](https://pglpm.github.io/inferno/reference/util_learntvar2sd.md)
  : Convert learnt with R/C/Dvar to learnt with R/C/Dsd
- [`util_lprobsargsyx()`](https://pglpm.github.io/inferno/reference/util_lprobsargsyx.md)
  : Prepare arguments for util_lprobsyx from data
- [`util_lprobsbase()`](https://pglpm.github.io/inferno/reference/util_lprobsbase.md)
  : Calculate collection of log-probabilities for different components
  and samples
- [`util_lprobsmi()`](https://pglpm.github.io/inferno/reference/util_lprobsmi.md)
  : Calculate pairs of log-probabilities for mutualinfo()
- [`util_prepPcheckpoints()`](https://pglpm.github.io/inferno/reference/util_prepPcheckpoints.md)
  : Format datapoints for testing of MCMC progress
- [`util_qYXcont()`](https://pglpm.github.io/inferno/reference/util_qYXcont.md)
  : Calculate quantiles for continuous Y by bisection
- [`util_qYXdiscr()`](https://pglpm.github.io/inferno/reference/util_qYXdiscr.md)
  : Calculate quantiles for discrete Y by bisection
- [`vtransform()`](https://pglpm.github.io/inferno/reference/vtransform.md)
  : Transforms variates to different representations
- [`workerfun()`](https://pglpm.github.io/inferno/reference/workerfun.md)
  : Worker function called by learn()
