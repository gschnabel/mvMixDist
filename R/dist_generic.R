setClass("probdist",
         contains="VIRTUAL")

#' Draw a sample
#'
#' @param dist a distribution object 
#' @param n number of samples to draw from distribution
#'
#' @return a \code{p x n} matrix, each column contains an observation vector
#' @export
#' 
setGeneric("getSample", 
           function(dist,n)
             standardGeneric("getSample"))

#' Get density
#'
#' @param dist a distribution object
#' @param x contains \code{p x n} matrix where columns are observation vectors
#' @param log if TRUE return logarithm of density
#'
#' @return a vector with the (log) densities for the observation vectors
#' @export
#'
setGeneric("getDensity", 
           function(dist, x, log=TRUE, vars=NULL)
             standardGeneric("getDensity"))

#' Maximum Likelihood Estimate
#'
#' @param dist a distribution object
#' @param x contains \code{p x n} matrix where columns are observation vectors
#' @param weights importance of each observation
#' @param ... distribution specific tuning parameters
#'
#' @return a distribution object with optimized distribution parameters
#' @export
#'
setGeneric("getMLEstimate",
           function(dist, x, weights=rep(1,ncol(x)), ...)
             standardGeneric("getMLEstimate")
)


#' Get Prior Density
#'
#' @param dist a distribution object
#' @param n number of samples
#' @param ... distribution specific tuning parameters
#'
#' @export
#'
setGeneric("getPriorDensity",
           function(dist, ...)
               standardGeneric("getPriorDensity")
)

#' Get Posterior Sample
#'
#' @param dist a distribution object
#' @param x contains \code{p x n} matrix where columns are observation vectors
#' @param n number of samples
#' @param ... distribution specific tuning parameters
#'
#' @return a sample from the posterior distribution (via Gibbs sampling)
#' @export
#'
setGeneric("getPosteriorSample",
           function(dist, x, n, ...)
             standardGeneric("getPosteriorSample")
)


#' Evaluate Bayesian information criterion (BIC)
#'
#' @param dist a distribution object
#' @param x a \code{p x n} matrix where each column is an observation vector
#' @param weights importance of each observation
#'
#' @return the value of the Bayesian information criterion
#' @note If not all weights are equal one, the effective sample size (ESS)
#'       is used as the number of data points
#' @export
#'
setGeneric("getBIC",
           function(dist, x, weights=rep(1,ncol(x)))
             standardGeneric("getBIC")
)


#' Get membership probability
#'
#' @param dist a distribution object
#' @param x contains \code{p x n} matrix where columns are observation vectors
#' @param log if TRUE returns logarithm of membership probabilities
#' 
#' @return a \code{m x n} matrix where each column is associated with an observation vector
#'         and contains for each of the \code{m} components the probability of membership 
#'      
#' @note Only distributions which are mixtures of distributions implement this function.
#' @export
#'
setGeneric("getMembership",
           function(dist, x, log=TRUE) 
             standardGeneric("getMembership")
)








