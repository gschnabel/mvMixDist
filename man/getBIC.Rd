% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist_generic.R
\name{getBIC}
\alias{getBIC}
\title{Evaluate Bayesian information criterion (BIC)}
\usage{
getBIC(dist, x, weights = rep(1, ncol(x)))
}
\arguments{
\item{dist}{a distribution object}

\item{x}{a \code{p x n} matrix where each column is an observation vector}

\item{weights}{importance of each observation}
}
\value{
the value of the Bayesian information criterion
}
\description{
Evaluate Bayesian information criterion (BIC)
}
\note{
If not all weights are equal one, the effective sample size (ESS)
      is used as the number of data points
}
