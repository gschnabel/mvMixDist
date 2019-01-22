#' Create Multivariate Normal Distribution
#'
#' @param mean the center vector
#' @param sigma the covariance matrix
#'
#' @return a distribution object of a multivariate normal distribution
#' @import MCMCpack mvtnorm
#' @export
#'
createDist_MVN <- function(mean, sigma, prior = list()) {
  L <- chol(sigma)
  logDetSigma <- 2*sum(log(diag(L)))
  invL <- backsolve(chol(sigma),diag(ncol(sigma)))
  new("mvndist",mean=mean, sigma=sigma, L=L, invL=invL,
      logDetSigma=logDetSigma, prior=prior)
}


setClass("mvndist",
         slots=list(
           mean="numeric",
           sigma="matrix",
           L="matrix",
           invL="matrix",
           logDetSigma="numeric",
           prior = 'list'
         ),
         prototype=list(
           mean=0,
           sigma=matrix(1),
           L=matrix(1),
           invL=matrix(1),
           logDetSigma=1,
           prior = list()
         ),
         validity=function(object) {
           all(dim(object@sigma)==length(object@mean))
         },
         contains="probdist")



setMethod("getSample",
          signature=list(
            dist = "mvndist",
            n = "numeric"
          ),
          definition=function(dist, n) {
            numDim <- length(dist@mean)
            randVecs <- matrix(rnorm(numDim*n),ncol=n)
            t(dist@L) %*% randVecs + dist@mean
          })


setMethod("getDensity",
          signature=list(
            dist = "mvndist",
            x = "matrix"
          ),
          definition=function(dist, x, log=TRUE, vars=NULL) {
            if (is.null(vars)) {
              numDim <- length(dist@mean)
              stopifnot(nrow(x)==numDim)
              logNormConst <- (-numDim/2)*log(2*pi) - 0.5*dist@logDetSigma
              res <- logNormConst - 0.5*colSums((t(dist@invL) %*% (x-dist@mean))^2)
              if (!log) res <- exp(res)
              res
            }
            else
            {
              subX <- x
              subMean <- dist@mean[vars]
              subSigma <- dist@sigma[vars,vars,drop=FALSE]
              marDist <- createDist_MVN(subMean,subSigma)
              getDensity(marDist, subX, log, vars=NULL)
            }
          })



setMethod("getMLEstimate",
          signature=list(
            dist = "mvndist",
            x = "matrix"
          ),
          definition=function(dist, x, weights) {
            pars <- cov.wt(t(x), wt=weights)
            createDist_MVN(pars$center, pars$cov)
          })

setMethod("getPosteriorSample",
          signature=list(
            dist = "mvndist",
            x = "matrix",
            n = "numeric"
          ),
          definition=function(dist, x, n) {
            
            # prior
            prior <- dist@prior
            mu0 <- prior$mu0
            kappa0 <- prior$kappa0
            nu0 <- prior$nu0
            phi <- prior$phi
            # data
            numObs <- ncol(x)
            D <- cov.wt(t(x), method = "ML")
            # posterior
            mu1 <- (kappa0*mu0 + numObs*D$center) / (kappa0 + numObs)
            kappa1 <- kappa0 + numObs
            nu1 <- nu0 + numObs
            phi1 <- phi + numObs*D$cov + (kappa0*numObs) / (kappa0+numObs) * 
              tcrossprod(D$center - mu0)
            # create samples
            distList <- replicate(n, NULL, simplify = FALSE)
            for (i in seq(n)) {
              sampleCov <- riwish(nu1, phi1)
              sampleMean <- as.vector(rmvnorm(1, mu1, 1/kappa1 * sampleCov))
              distList[[i]] <- createDist_MVN(sampleMean, sampleCov, prior = dist@prior)
            }
            distList
          })

