% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_mcmcstop.R
\name{funMCEI}
\alias{funMCEI}
\title{Calculate quantile width through batches}
\usage{
funMCEI(x, fn, p = c(0.055, 0.945), ...)
}
\arguments{
\item{x}{A matrix, rows being MC samples and columns being quantities whose MCSE is to be estimated.}
}
\value{
Estimates of the MC standard error for each trace. Division by sqrt(N) is already performed.
}
\description{
Modified from
from https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/ESS.R
}
\keyword{internal}
