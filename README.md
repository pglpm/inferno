*[If you're looking for the package "inferno", you've come to the right place! It has been renamed "prova"]*

<img src="https://github.com/pglpm/prova/raw/main/development/manual/prova_symbol.png" alt="Ensemble of densities" width="100%"/>

# ***prova***: Probabilistic-statistical analysis and inference with Bayesian nonparametrics

This repository provides an R package and some theoretical background for *Bayesian nonparametric population inference*, which can also be called "inference under exchangeability" or "density inference". The package is especially apt for the study of statistics and associations of subpopulations or subgroups. The [introductory vignette](https://pglpm.github.io/prova/articles/start.html) explains with an example the main ideas and functions, in particular for researchers who are more familiar with traditional "frequentist" statistics but would like to try the Bayesian approach; see the [post](https://www.apadivisions.org/division-7/publications/newsletters/developmental/2018/07/bayesian-statistics) by Barbara W. Sarnecka, frequentist statistician turned Bayesian, for a brilliant overview of the Bayesian advantages.

The package is under development and has not yet reached a stable phase: function names and arguments may still change, new functions will be added. More tutorials will be prepared. Also the package name is still under consideration.

But the core functionalities do work, and have been tested in concrete research questions; see [example applications](#example-applications) below.

The package internally does the computations necessary for Bayesian inference by means of Monte Carlo methods, thanks to the R package [**Nimble**](https://r-nimble.org/). Users unfamiliar with Monte Carlo methods don't have to worry, because the computations are handled automatically. Users familiar with Monte Carlo methods can easily have access to computational details and can even change some of the computation hyperparameters.

If you want to test the package we'd be very happy to help in resolving possible issues and in understanding the functionalities.

## Installation

You need to have the package [**Nimble**](https://r-nimble.org/), *at least version 1.4.2*, installed. Please follow its [installation instructions](https://r-nimble.org/manual/cha-installing-nimble.html) for your operating system.

You can then install **prova** in R by using the `remotes` package:
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

The vignette [*Bayesian nonparametric inference with **prova***](https://pglpm.github.io/prova/articles/start.html) is a step-by-step introduction to **prova** and also to Bayesian nonparametrics. It guides you through a concrete example with various kinds of inferences. You may also try to follow it using a dataset of your own.

Other tutorials, still drafts, are available at [pglpm.github.io/prova](https://pglpm.github.io/prova)

A summary of the theoretical foundations, including further references, is available in [this draft](https://github.com/pglpm/prova/raw/main/development/manual/optimal_predictor_machine.pdf). The main idea for the internal mathematical representation comes from [Dunson & Bhattacharya](https://doi.org/10.1093/acprof:oso/9780199694587.003.0005) and [Ishwaran & Zarepour](https://doi.org/10.2307/3315951).

For a low-level course on Bayesian nonparametric inference and Decision Theory see [Foundations of data science](https://pglpm.github.io/ADA511).


## **prova** App
An application has been built upon **prova**. This app can be used for testing out the features of **prova** with just a few button clicks, without having to write any code in R yourself.

### Desktop Application
* Currently available for Windows and MacOS. Download and install the desktop application by following this: [Installation Guide](https://github.com/h587916/Inferno-App/releases/tag/1.0.2).

### Cross-Platform Open Source Version
* For Windows, macOS, and Linux, you can run the PySide6 app locally using Python by following this: [Setup Guide](https://github.com/h587916/Inferno-App?tab=readme-ov-file#inferno-app).

## Example applications

- [*Personalized prognosis & treatment using an optimal predictor machine: An example study on conversion from Mild Cognitive Impairment to Alzheimer's Disease*](https://doi.org/10.31219/osf.io/8nr56).

- [*Don't guess what's true: choose what's optimal. A probability transducer for machine-learning classifiers*](https://doi.org/10.31219/osf.io/vct9y)

- [*Does the evaluation stand up to evaluation? A first-principle approach to the evaluation of classifiers*](https://doi.org/10.31219/osf.io/7rz8t)

- [*Calibrated and uncertain? Evaluating uncertainty estimates in
binary classification models*](https://doi.org/10.1088/2632-2153/ae45ed)


Projects using **prova**:

- [InfernoCalibNet](https://m4siko.github.io/InfernoCalibNet).
- [parkinsonbayes](https://github.com/pglpm/parkinsonbayes).



## Contact

Please report bugs and request features or specific documentation on [GitHub Issues](https://github.com/pglpm/prova/issues).
If you have other questions about application, theory, technical implementation, feel free to contact Luca <pglXYZ@portamanaXYZ.org> (remove 'XYZ' for anti-spam purposes).


## Disclaimer

No large language models were used in the production of this software and of its documents.
