# Prova v1.0.0

* Improved documentation, also of internal functions.
* Introduced argument `verbose =` in several function, so as to show information to console only if the user so desires. This argument is `FALSE` for all functions except `learn()`.
* Added example of imputation in vignette.
* Rewritten README with description of the package's main features.
* Cleaned up code.

# Prova v0.9.0

* Improved documentation, also of internal functions.
* Removed unused functions.
* Cleaned up code.


# Prova v0.7.0

* Initial CRAN submission.
* Introduced the `print()` method for probability objects.
* Slightly modified output of `Pr()` and related functions.
* Added references.
* Modified vignettes.


# Prova v0.6.5

* Added examples to all functions visible to the user.
* Fixed some bugs (in particular a bug that would allow the input of a non-existent nominal value giving erroneous probabilities).
* Corrected typos.

# Prova v0.6.0

Renamed the package from **inferno** to **Prova**. It turns out, unsurprisingly, that the name "inferno" was already in use for several inference-related software packages of various kinds.

# Inferno v0.5.5

* Removed a bug that broke `learn()` with datasets consisting of only one nominal variate.
* Fix bug in handling metadata.


# Inferno v0.5.0

* Added function `qPr()` for the computation of quantiles and their variability.

* Corrected some bugs.

* **NB: learnt objects created with versions < 0.5.0 are incompatible with the functions of version 0.5.0.** To convert to the new version, use the function `util_learntvar2sd(file)`, where file is the path of the learnt object to be converted.



# Inferno v0.3.2

* Improved graphic output of `flexiplot()`.


# Inferno v0.3.1

* Added the function `rPr()`, which generates datapoints for any desired set of joint variates, according to the posterior probability calculated with the `learn()` function. See documentation.
* Added arguments xjitter and yjitter to `flexiplot()`, useful for scatterplots of discrete variates.
* Updated the `mutualinfo()` function, which should also be a little faster.
* Updated and corrected documentation.


# Inferno v0.3.0

**NB: this release makes all relevant functions incompatible with objects obtained with previous releases.** Please submit an issue if you'd like to convert your previous results in a format compatible with the new release. A conversion utility will be made available soon if there are enough requests.

* The `Pr()` function has a new argument `tails =`, and now accepts arbitrary combinations with point-value arguments `(Y = y)` and left- or right-open interval arguments (`Y <= y` and `Y >= y`), the latter for ordinal and continuous variates only. Thus it covers and extends the use of the now-obsolete function `tailPr()`. See documentation, especially about the new argument tails.
* The `Pr()` function now outputs two new elements: values.MCerror and quantiles.MCerror, quantifying the accuracy of the Monte Carlo calculation of the values and quantiles elements. See documentation.
* New handling of ordinal and nominal variates, which should be faster and use slightly less memory.
* More precise calculation of probabilities for rounded and discrete variates.


# Inferno v0.2.3

* Improved (hopefully) stopping rule of the Markov-chain Monte Carlo computation. Now partly based on the "bulk ESS" function from Vehtari & al.
* Mainly for debugging purposes, `learn()` now continuously updates the Monte Carlo trace plot during calculations.
* A couple more internal functions used for debugging and Monte Carlo monitoring.


# Inferno v0.2.1

Updates to GitHub:
Added GitHub Actions workflow for automatic testing of the software.

Updates to code:

* New logical argument "verbose" (def. TRUE) in `buildmetadata()`. When TRUE, messages are given for each variate, explaining the internal heuristics and guessing to determine the various metadata values.

* Modified handling of rounded continuous variates, now more consistent according to discussion in issue #50.

* Elimination of type-"L" variates in Monte Carlo sampling. The type "D" handles both rounded continuous variates and ordinal variates having domain with more than 10 values. `samplesFdistribution()` and other functions have been updated accordingly.

* Rewritten `plotFsamples()`. Now it goes through every variate type in turn, and should be easier to understand.

* Modified the information contained in the internal "auxmetadata" object, and accordingly modified all functions that use this object.

Performed tests:

* Performed a battery of tests against many datasets available in base-R. This lead to the unveiling and fixing of several small bugs.

The tests were performed to check the working of `buildmetadata()`, `buildauxmetadata()`, `samplesFdistribution()`, `plotFsamples()`.

* With the mentioned datasets, `samplesFdistribution()` has been checked against a clearer (but much slower), for-loop-based script -- written from scratch -- to calculate the various probabilities. This script also uses mathematical formulae that are theoretically identical but numerically different when it comes to finite-precision arithmetic. Some bugs have been fixed

* The latter test also shows that errors coming from finite-precision arithmetic are all below 10^-15.


# Inferno v0.2.0

Major changes:

* New initialization procedure for Monte Carlo sampling (this led to great improvements)

* New stopping rule for the Monte Carlo sampling (partly based on the ideas in doi.org/10.1080/10618600.2015.1044092)

* New handling of rounded and ordinal variates

* New argument 'auxdata' in `inferpopulation()`: the user can here give a much larger dataset (of which 'data' argument is presumably a subset), which is used to calculate some general statistics to improve the inference. The idea is that 'auxdata' cannot be used for the Monte Carlo proper, owing to memory or time limitations, but at least we can squeeze some other useful information out of it.

* New argument 'relerror' in `inferpopulation()`: an (approximate) upper bound to the desired numerical error. It's the numerical error relative to the width of the probability distribution.

* New argument 'ncheckpoints' in `inferpopulation()`: number of datapoints to be used to check Monte Carlo convergence. NULL value (default) is equal to the number of variates + 1.

Minor changes:

* Elimination of fields "centralvalue", "lowvalue", "highvalue" in metadata - they aren't required anymore

* Elimination of not-used diagnostics in log files

* Elimination of dependence on LaplacesDemon

* Corrections and improvements of plots and log information

* Correction in how parallel cores were closed (was leading to errors)

* Various bug fixes

# Inferno v0.1.0 - First package release

The code is now available as an R package on GitHub.

This release includes many changes to the code.

* Changes in file and folder structures
* Documentation has been added
* `mutualinfo()` has been added. This function calculates mutual information between groups of joint variates.
