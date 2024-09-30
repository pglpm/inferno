<img src="https://github.com/pglpm/inferno/raw/main/development/manual/inferno_symbol.png" alt="Ensemble of densities" width="100%"/>

# Bayesian nonparametric inference

[![Documentation](https://img.shields.io/badge/Documentation-blue)](https://pglpm.github.io/inferno/)

This repository provides an R package and some theoretical background for Bayesian nonparametric inference under exchangeability, or "inference about populations".
The package is under rapid development and has not reached a stable phase. This means that function names and arguments may still change. The package name is also still under consideration. However, the core functionalities and probability calculations work. While the code is still in its '0.X' phase, we recommend contacting the developers if you want to start using the package for a research project. We would love more "beta-testers"!

## Contact

Please report bugs and request features or specific documentation on [GitHub Issues](https://github.com/pglpm/inferno/issues).
If you have other questions feel free to contact the developers:
* For application/theory: PierGianLuca Porta Mana <pgl@portamana.org>
* For technical implementation: Aurora Grefsrud <agre@hvl.no>

## Documentation
View the documentation at [pglpm.github.io/inferno](https://pglpm.github.io/inferno/). As the code is still very much in the development phase the documentation is also a work in progress.
Tutorials will be posted as soon as possible. Feel free to take a look at the [draft on the theoretical foundation](https://github.com/pglpm/inferno/raw/main/development/manual/optimal_predictor_machine.pdf). A concrete example application in medicine is given in [this paper](https://doi.org/10.31219/osf.io/8nr56); example applications in machine learning are [this paper](https://doi.org/10.31219/osf.io/7rz8t) and [this paper](https://doi.org/10.31219/osf.io/vct9y).

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

## Further reading
For a low-level course on Bayesian nonparametric population inference see [Foundations of data science](https://pglpm.github.io/ADA511/).
