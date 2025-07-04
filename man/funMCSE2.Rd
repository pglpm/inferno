% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_mcmcstop.R
\name{funMCSE2}
\alias{funMCSE2}
\title{Function for calculating MC standard error, from Geyer's mcmc package}
\usage{
funMCSE2(x)
}
\arguments{
\item{x:}{a matrix, rows being MC samples and columns being quantities whose MCSE is to be estimated.}
}
\value{
Estimates of the MC standard error for each trace. Division by sqrt(N) is already performed.
}
\description{
Function for calculating MC standard error, from Geyer's mcmc package
}
\keyword{internal}
