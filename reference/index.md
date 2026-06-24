# Package index

## Main functions

- [`Pr()`](https://pglpm.github.io/prova/reference/Pr.md) : Calculate
  posterior probabilities

- [`flexiplot()`](https://pglpm.github.io/prova/reference/flexiplot.md)
  : Plot numeric or character values

- [`hist(`*`<probability>`*`)`](https://pglpm.github.io/prova/reference/hist.probability.md)
  : Plot the variability of an object of class "probability" as a
  histogram

- [`learn()`](https://pglpm.github.io/prova/reference/learn.md) : Monte
  Carlo computation of posterior probability distribution

- [`learntExample`](https://pglpm.github.io/prova/reference/learntExample.md)
  :

  Example `learnt` object produced by learn()

- [`metadataExample`](https://pglpm.github.io/prova/reference/metadataExample.md)
  : Example metadata file

- [`metadatatemplate()`](https://pglpm.github.io/prova/reference/metadatatemplate.md)
  : Metadata and helper function for metadata

- [`mutualinfo()`](https://pglpm.github.io/prova/reference/mutualinfo.md)
  : Calculate mutual information between groups of joint variates

- [`plot(`*`<probability>`*`)`](https://pglpm.github.io/prova/reference/plot.probability.md)
  : Plot an object of class "probability"

- [`plotquantiles()`](https://pglpm.github.io/prova/reference/plotquantiles.md)
  : Plot pairs of quantiles

- [`print(`*`<probability>`*`)`](https://pglpm.github.io/prova/reference/print.probability.md)
  : Print an object of class "probability"

- [`pwrite.csv()`](https://pglpm.github.io/prova/reference/prova.data.md)
  [`pread.csv()`](https://pglpm.github.io/prova/reference/prova.data.md)
  :

  Write and read CSV files in **Prova**

- [`qPr()`](https://pglpm.github.io/prova/reference/qPr.md) : Calculate
  quantiles

- [`rPr()`](https://pglpm.github.io/prova/reference/rPr.md) : Generate
  datapoints

- [`vrtgrid()`](https://pglpm.github.io/prova/reference/vrtgrid.md) :
  Create a grid of values for a variate

## Internal functions

For developers.

- [`buildauxmetadata()`](https://pglpm.github.io/prova/reference/buildauxmetadata.md)
  : Build preliminary metadata flie

- [`createQfunction()`](https://pglpm.github.io/prova/reference/createQfunction.md)
  : Calculate and save transformation function for ordinal variates

- [`fftNGS()`](https://pglpm.github.io/prova/reference/fftNGS.md) : Find
  optimal FFT size

- [`funAC()`](https://pglpm.github.io/prova/reference/funAC.md) :
  Compute autocovariance

- [`funESS3()`](https://pglpm.github.io/prova/reference/funESS3.md) :
  Compute ESS

- [`funESSLD()`](https://pglpm.github.io/prova/reference/funESSLD.md) :
  Calculate MC effective sample size using LaplacesDemon's algorithm

- [`funMCEI()`](https://pglpm.github.io/prova/reference/funMCEI.md) :
  Calculate quantile width through batches

- [`funMCEQ()`](https://pglpm.github.io/prova/reference/funMCEQ.md) :
  Calculate credibility quantiles on estimated quantile

- [`funMCSELD()`](https://pglpm.github.io/prova/reference/funMCSELD.md)
  : Calculate MC standard error using LaplacesDemon's batch means

- [`learnbind()`](https://pglpm.github.io/prova/reference/learnbind.md)
  : Bind 3D arrays by first dimension

- [`mcjoin()`](https://pglpm.github.io/prova/reference/mcjoin.md) :
  Concatenate mcsample objects

- [`mcsubset()`](https://pglpm.github.io/prova/reference/mcsubset.md) :
  Eliminate samples from mcsamples object

- [`plotFsamples()`](https://pglpm.github.io/prova/reference/plotFsamples.md)
  : Plot one-dimensional posterior probabilities

- [`prova`](https://pglpm.github.io/prova/reference/prova-package.md)
  [`prova-package`](https://pglpm.github.io/prova/reference/prova-package.md)
  : prova: Nonparametric Probabilistic-Statistical Variate Analysis with
  Automated Markov-Chain Monte Carlo

- [`prsubset()`](https://pglpm.github.io/prova/reference/prsubset.md) :

  Subset variates of an object of class `probability`

- [`rowcumsum()`](https://pglpm.github.io/prova/reference/rowcumsum.md)
  : Cumulative sum along first dimension

- [`rowinvcumsum()`](https://pglpm.github.io/prova/reference/rowinvcumsum.md)
  : Inverse cumulative sum along first dimension

- [`util_Pcheckpoints()`](https://pglpm.github.io/prova/reference/util_Pcheckpoints.md)
  : Calculate joint frequencies for checkpoints in learn()

- [`util_combineYX()`](https://pglpm.github.io/prova/reference/util_combineYX.md)
  : Calculate probabilities, quantiles, etc, for all Y and X
  combinations

- [`util_denorm()`](https://pglpm.github.io/prova/reference/util_denorm.md)
  : Utility function to avoid finite-precision accuracys

- [`util_joinPtraces()`](https://pglpm.github.io/prova/reference/util_joinPtraces.md)
  : Join '\_\_\_\_tempPtraces-' files

- [`util_learntvar2sd()`](https://pglpm.github.io/prova/reference/util_learntvar2sd.md)
  : Convert learnt with R/C/Dvar to learnt with R/C/Dsd

- [`util_lprobsargsyx()`](https://pglpm.github.io/prova/reference/util_lprobsargsyx.md)
  : Prepare arguments for util_lprobsyx from data

- [`util_lprobsbase()`](https://pglpm.github.io/prova/reference/util_lprobsbase.md)
  : Calculate collection of log-probabilities for different components
  and samples

- [`util_lprobsmi()`](https://pglpm.github.io/prova/reference/util_lprobsmi.md)
  : Calculate pairs of log-probabilities for mutualinfo()

- [`util_prepPcheckpoints()`](https://pglpm.github.io/prova/reference/util_prepPcheckpoints.md)
  : Format datapoints for testing of MCMC progress

- [`util_qYXcont()`](https://pglpm.github.io/prova/reference/util_qYXcont.md)
  : Calculate quantiles for continuous Y by bisection

- [`util_qYXdiscr()`](https://pglpm.github.io/prova/reference/util_qYXdiscr.md)
  : Calculate quantiles for discrete Y by bisection

- [`vtransform()`](https://pglpm.github.io/prova/reference/vtransform.md)
  : Transforms variates to different representations

- [`workerfun()`](https://pglpm.github.io/prova/reference/workerfun.md)
  : Worker function called by learn()
