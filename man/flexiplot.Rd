% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_tplotfunctions.R
\name{flexiplot}
\alias{flexiplot}
\title{Plot numeric or character values}
\usage{
flexiplot(
  x,
  y,
  xdomain = NULL,
  ydomain = NULL,
  xlim = NULL,
  ylim = NULL,
  type = "l",
  pch = c(1, 2, 0, 5, 6, 3),
  lty = c(1, 2, 4, 3, 6, 5),
  lwd = 2,
  col = palette(),
  alpha.f = 1,
  xjitter = NULL,
  yjitter = NULL,
  grid = TRUE,
  cex.main = 1,
  add = FALSE,
  ...
)
}
\arguments{
\item{x}{Numeric or character: vector of x-coordinates. If missing, a numeric vector \code{1:...} is created having as many values as the rows of \code{y}.}

\item{y}{Numeric or character: vector of y coordinates. If missing, a numeric vector \code{1:...} is created having as many values as the rows of \code{x}.}

\item{xdomain, ydomain}{Character or numeric or \code{NULL} (default): vector of possible values of the variables represented in the \code{x}- and \code{y}-axes, in case the \code{x} or \code{y} argument is a character vector. The ordering of the values is respected. If \code{NULL}, then \code{unique(x)} or \code{unique(y)} is used.}

\item{xlim, ylim}{\code{NULL} (default) or a vector of two values. In the latter case, if any of the two values is not finite (including \code{NA} or \code{NULL}), then the \code{min} or \code{max} \code{x}- or \code{y}-coordinates of the plotted points are used.}

\item{alpha.f}{Numeric, default 1: opacity of the colours, \code{0} being completely invisible and \code{1} completely opaque.}

\item{xjitter, yjitter}{Logical or \code{NULL} (default): add \code{\link[base:jitter]{base::jitter()}} to \code{x}- or \code{y}-values? Useful when plotting discrete variates. If \code{NULL}, jitter is added if the values are of character class.}

\item{grid}{Logical: whether to plot a light grid. Default \code{TRUE}.}

\item{...}{Other parameters to be passed to \code{\link[graphics:matplot]{graphics::matplot()}}.}
}
\description{
Plot function that modifies and expands the \strong{graphics} package's \code{\link[graphics:matplot]{graphics::matplot()}} function in several ways.
}
\details{
Some of the additional features provided by \code{flexiplot} are the following. First, either or both \code{x} and \code{y} arguments can be of class \code{\link[base:character]{base::character}}. In this case, axes labels corresponding to the unique values are used (see arguments \code{xdomain} and \code{ydomain}). A jitter can also be added to the generated points, via the \code{xjitter} and \code{yjitter} switches. Second, it allows for the specification of only a lower or upper limit in the \code{xlim} and \code{ylim} arguments. Third, it uses a cleaner plotting style and a default argument \code{type = 'l'} (line plot) rather than \code{type = 'p'} (point plot).
}
