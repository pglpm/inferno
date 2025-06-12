<img src="https://github.com/pglpm/inferno/raw/main/development/manual/inferno_symbol.png" alt="Ensemble of densities" width="100%"/>

# Inferno: ***Infer***ence in ***R*** with Bayesian ***no***nparametrics


This repository provides an R package and some theoretical background for Bayesian nonparametric inference under exchangeability, or "inference about populations".
The package is under rapid development and has not reached a stable phase. This means that function names and arguments may still change. More tutorials will be added. The package name is also still under consideration.

The core functionalities do work, however, and have been tested on concrete research questions; see references below. If you want to test the package, he developer would be very happy to help in resolving possible issues and in better understanding the functionalities.

## Documentation
Documentation and some tutorials are available at [pglpm.github.io/inferno](https://pglpm.github.io/inferno). These are also work in progress.

A summary of the theoretical foundations, including further references, is available in [this draft](https://github.com/pglpm/inferno/raw/main/development/manual/optimal_predictor_machine.pdf). The main idea for the internal mathematical representation comes from Dunson & Bhattacharya: [*Nonparametric Bayes regression and classification through mixtures of product kernels*](https://doi.org/10.1093/acprof:oso/9780199694587.003.0005).

For a low-level course on Bayesian nonparametric inference see [Foundations of data science](https://pglpm.github.io/ADA511).


## Example applications

- [*Personalized prognosis & treatment using an optimal predictor machine: An example study on conversion from Mild Cognitive Impairment to Alzheimer's Disease*](https://doi.org/10.31219/osf.io/8nr56).

- [*Don't guess what's true: choose what's optimal. A probability transducer for machine-learning classifiers*](https://doi.org/10.31219/osf.io/vct9y)

- [*Does the evaluation stand up to evaluation? A first-principle approach to the evaluation of classifiers*](https://doi.org/10.31219/osf.io/7rz8t)


Projects using ***Inferno***:

- [InfernoCalibNet](https://m4siko.github.io/InfernoCalibNet).
- [parkinsonbayes](https://github.com/pglpm/parkinsonbayes).

## Installation
Install the package with R by using the `remotes` package:
```
remotes::install_github('pglpm/inferno')
```
To install a tagged version:
```
remotes::install_github('pglpm/inferno@vx.y.z')
```
To install from source, first clone the repo:
```
git clone https://github.com/pglpm/inferno.git
```
then install the package in R:
```
install.packages(pkgs='path/to/inferno', repos=NULL)
```
the installation will automatically also install all required R-dependencies.

## Inferno App
An application has been built upon ***Inferno***. This app can be used for testing out the features of ***Inferno*** with just a few button clicks, without having to write any code in R yourself. 

### Desktop Application
* Currently available for Windows and MacOS. Download and install the desktop application by following this: [Installation Guide](https://github.com/h587916/Inferno-App/releases/tag/1.0.2).

### Cross-Platform Open Source Version
* For Windows, macOS, and Linux, you can run the PySide6 app locally using Python by following this: [Setup Guide](https://github.com/h587916/Inferno-App?tab=readme-ov-file#inferno-app).


## Contact

Please report bugs and request features or specific documentation on [GitHub Issues](https://github.com/pglpm/inferno/issues).
If you have other questions about application, theory, technical implementation, feel free to contact Luca <pglXYZ@portamanaXYZ.org> (remove 'XYZ' for anti-spam purposes).
