*[Looking for the R package "inferno"? you came to the right place! It's been renamed "Prova"]*

<img src="https://github.com/pglpm/prova/raw/main/man/figures/prova_symbol.jpg" alt="Ensemble of densities" width="100%"/>

"prova" /'prɔva/ (Italian)

- [(noun)](https://dictionary.cambridge.org/dictionary/italian-english/prova): test, trial, assessment, proof, evidence, sign, indication, try, attempt.
- [(verb)](https://dictionary.cambridge.org/dictionary/italian-english/provare): test!, try out!, assess!, attempt!, prove!, demonstrate!, show!


# ***Prova***: Probabilistic-statistical variate analysis and inference with Bayesian nonparametrics

This repository provides an R package and some theoretical background for *Bayesian nonparametric population inference*, which can also be called "inference under exchangeability" or "density inference". The package is especially apt for the study of statistics and associations of subpopulations or subgroups. The [introductory vignette](https://pglpm.github.io/prova/articles/start.html) explains with an example the main ideas and functions, in particular for researchers who are more familiar with traditional "frequentist" statistics but would like to try the Bayesian approach; see the [post](https://www.apadivisions.org/division-7/publications/newsletters/developmental/2018/07/bayesian-statistics) by Barbara W. Sarnecka, frequentist statistician turned Bayesian, for a brilliant overview of the Bayesian advantages.

The package is under development and has not yet reached a stable phase: function names and arguments may still change, new functions will be added. More tutorials will be prepared. Also the package name is still under consideration.

But the core functionalities do work, and have been tested in concrete research questions; see [example applications](#example-applications) below.

The package internally does the computations necessary for Bayesian inference by means of Monte Carlo methods, thanks to the R package [**Nimble**](https://r-nimble.org/). Users unfamiliar with Monte Carlo methods don't have to worry, because the computations are handled automatically. Users familiar with Monte Carlo methods can easily have access to computational details and can even change some of the computation hyperparameters.

If you want to test the package we'd be very happy to help in resolving possible issues and in understanding the functionalities.

## Installation

You need to have the package [**Nimble**](https://r-nimble.org/), *at least version 1.4.2*, installed. Please follow its [installation instructions](https://r-nimble.org/manual/cha-installing-nimble.html) for your operating system.

You can then install **Prova** in R by using the `remotes` package:
```
remotes::install_github('pglpm/prova')
```

To install a tagged version:
```
remotes::install_github('pglpm/prova@vx.y.z')
```

To install from source, first clone the repo:
```
git clone https://github.com/pglpm/prova.git
```

then install the package in R:

```
install.packages(pkgs='path/to/prova', repos=NULL)
```

the installation will also automatically install all required R-dependencies.


## Documentation

The vignette [*Bayesian nonparametric inference with **Prova***](https://pglpm.github.io/prova/articles/start.html) is a step-by-step introduction to **Prova** and also to Bayesian nonparametrics. It guides you through a concrete example with various kinds of inferences. You may also try to follow it using a dataset of your own.

Other tutorials are available at [pglpm.github.io/prova](https://pglpm.github.io/prova/)

A summary of the theoretical foundations, including further references, is available in [this draft](https://github.com/pglpm/prova/raw/main/development/manual/optimal_predictor_machine.pdf). The main idea for the internal mathematical representation comes from [Dunson & Bhattacharya](https://doi.org/10.1093/acprof:oso/9780199694587.003.0005) and [Ishwaran & Zarepour](https://doi.org/10.2307/3315951).

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
