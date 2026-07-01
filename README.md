<!-- badges: start -->
  [![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.19136582.svg)](https://doi.org/10.5281/zenodo.17226082)
<!-- badges: end -->

*[Looking for the R package "inferno"? you came to the right place! It's been renamed "Prova"]*

<img src="https://github.com/pglpm/prova/raw/main/man/figures/prova_symbol.jpg" alt="Ensemble of densities" width="100%"/>

"prova" /'prɔva/ (Italian)

- [(noun)](https://dictionary.cambridge.org/dictionary/italian-english/prova): test, trial, assessment, proof, evidence, sign, indication, try, attempt.
- [(verb)](https://dictionary.cambridge.org/dictionary/italian-english/provare): test!, try out!, assess!, attempt!, prove!, demonstrate!, show!


# ***Prova***: probabilistic-statistical variate analysis and inference, nonparametric and with automated Markov-chain Monte Carlo

This repository provides an R package and some theoretical background to perform probabilistic and statistical analysis of data. These are the main features:

- Any combination of **binary**, **nominal**, **ordinal**, **continuous** variates. Continuous variates can be bounded or unbounded, and also rounded or discretized.
- **No modelling assumptions** such as gaussianity, linearity, or any other kind of model. The analysis and inferences are fully *non-parametric*.
- **No assumptions about functional dependence** between variates. The analysis and inferences are therefore more general than those by neural networks, random forests, or similar machine-learning algorithms.
- **Automatic imputation** of missing data: all sample data are used, even those that lacks some variate values. The imputation is done with a principled method (the marginalization rule of probability theory), rather than ad-hoc procedures.
- Easy and straightforward **subgroup analyses** and **stratified analyses**, for any division of variates, with full statistical details.
- **Quantification of generalization** beyond the finite sample size. In other words, quantification of uncertainty of results regarding the whole, unsampled, population.
- Straightforward use within **decision theory**, such as **clinical decision-making**. Users can immediately combine the probabilistic results with any measures of utilities, such as [quality-adjusted life years](https://toolbox.eupati.eu/glossary/quality-adjusted-life-year/).
- **Quantification of associations** between any kinds of variates, without modelling assumptions (gaussianity, linearity, etc.), thanks to the use of [*mutual information*](https://electropedia.org/iev/iev.nsf/display?openform&ievref=171-07-26).
- **Automated Markov-chain Monte Carlo** computation. Users unfamiliar with Monte Carlo methods don't have to worry, because the computations are handled automatically.

The package at bottom does Bayesian nonparametric inference (also called "density inference" or "inference under exchangeability"), which makes all features above possible.

The [introductory vignette](https://pglpm.github.io/prova/articles/intro.html) explains, with a guided example, most of the features above, as well as the main ideas and functions. It can be particularly useful for researchers who are more familiar with traditional "frequentist" statistics but would like to try the Bayesian approach. See the [post](https://www.apadivisions.org/division-7/publications/newsletters/developmental/2018/07/bayesian-statistics) by Barbara W. Sarnecka, frequentist statistician turned Bayesian, for a brilliant overview of the Bayesian advantages. The [vignette about mutual information](https://pglpm.github.io/prova/articles/mutualinfo.html) explains the use of this powerful measure of association.

The package is under continuous development, but the core functionalities work and have been tested in concrete research questions; see [example applications](#example-applications) below.

The package internally does the computations necessary for Bayesian inference by means of Monte Carlo methods thanks to the R package [**Nimble**](https://r-nimble.org/). As already mentioned, this computation is automated. Users familiar with Monte Carlo methods can still access computational details and can even change some of the computation hyperparameters.

## Installation

**You need to have installed the package [**Nimble**](https://r-nimble.org/), *at least version 1.4.2*.** Please follow its [installation instructions](https://r-nimble.org/manual/cha-installing-nimble.html) for your operating system.

Newer versions of **Prova** can be installed with
```
remotes::install_github('pglpm/prova')
```


## Documentation

The vignette [*Bayesian nonparametric inference with **Prova***](https://pglpm.github.io/prova/articles/intro.html) is a step-by-step introduction to **Prova** and also to Bayesian nonparametrics. It guides you through a concrete example with various kinds of inferences. You may also try to follow it using a dataset of your own.

Other tutorials are available at [pglpm.github.io/prova](https://pglpm.github.io/prova/), or can be accessed in an R session with `browseVignettes('prova')`.

A summary of the theoretical foundations, including further references, is available in [this draft](https://github.com/pglpm/prova/raw/main/development/manual/pglpm2024-bayes_nonparam.pdf). The main idea for the internal mathematical representation comes from [Dunson & Bhattacharya](https://doi.org/10.1093/acprof:oso/9780199694587.003.0005) and [Ishwaran & Zarepour](https://doi.org/10.2307/3315951).

For a low-level course on Bayesian nonparametric inference and Decision Theory see [Data Science and AI Prototyping](https://pglpm.github.io/ADA511/).


## Example applications

- [*Personalized prognosis & treatment using an optimal predictor machine: An example study on conversion from Mild Cognitive Impairment to Alzheimer's Disease*](https://doi.org/10.31219/osf.io/8nr56).

- [*Don't guess what's true: choose what's optimal. A probability transducer for machine-learning classifiers*](https://doi.org/10.31219/osf.io/vct9y)

- [*Does the evaluation stand up to evaluation? A first-principle approach to the evaluation of classifiers*](https://doi.org/10.31219/osf.io/7rz8t)

- [*Calibrated and uncertain? Evaluating uncertainty estimates in
binary classification models*](https://doi.org/10.1088/2632-2153/ae45ed)


## Projects using **Prova**:

- [InfernoCalibNet](https://m4siko.github.io/InfernoCalibNet/): uncertainty-aware predictions for medical AI using CNN and Bayesian nonparametrics framework
- [parkinsonbayes](https://github.com/pglpm/parkinsonbayes/): Examples of Bayesian nonparametric inference for studies of Parkinson's Disease
- [Inferno-App](https://github.com/Myddis/Inferno-App/): PySide6 application that integrates Python and R functionality using the **Inferno** (old version of **Prova**) R package.


## Contact

Please report bugs and request features or specific documentation on [GitHub Issues](https://github.com/pglpm/prova/issues).
If you have other questions about application, theory, technical implementation, feel free to contact Luca <pglXYZ@portamanaXYZ.org> (remove 'XYZ' for anti-spam purposes).


## Disclaimer

No large language models were used in the production of this software and of its documents.
