% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_mcmcstop.R
\name{funMCSELD}
\alias{funMCSELD}
\title{Calculate MC standard error using LaplacesDemon's batch means}
\usage{
funMCSELD(x)
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
