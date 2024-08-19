<img src="development/manual/opm_symbol.png" alt="Ensemble of densities" width="100%"/>

# Bayesian nonparametric inference

[![Documentation](https://img.shields.io/badge/Documentation-blue)](https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/SusyFitter)

This repository provides an R package and some theoretical background for Bayesian nonparametric inference under exchangeability, or "inference about populations".
The package itself is under rapid development, and has not reached a stable phase yet. However, the core functionalities are working, and the package can be used as is. While the code is still in it's '0.X' phase, we recommend that you contact the developers if you want to start using the package for a research project.

## Contact

Please report bugs, request features and ask for specific documentation on [GitHub Issues](https://github.com/pglpm/bayes_nonparametric_inference/issues).
If you have other questions feel free to contact the developers:
* For application/theory: PierGianLuca Porta Mana <pgl@portemana.org>
* For technical implementation: Aurora Grefsrud <agre@hvl.no>

## Documentation
View the documentation at [pglpm.github.io/bayes_nonparametric_inference](https://pglpm.github.io/bayes_nonparametric_inference/index.html). As the code is still very much in the development phase the documentation is also a work in progress.
Tutorials will be posted as soon as possible. Feel free to take a look at the [draft on the theoretical foundation](https://github.com/pglpm/bayes_nonparametric_inference/blob/main/development/manual/optimal_predictor_machine.pdf). A concrete example application, in medicine, is given in [this paper](https://doi.org/10.31219/osf.io/8nr56).

## Installation
Install the package with R by using the `remotes` package:
```
remotes::install_github('pglpm/bayes_nonparametric_inference')
```


## Further reading
For a low-level course on Bayesian nonparametric population inference see [Foundations of data science](https://pglpm.github.io/ADA511/).
