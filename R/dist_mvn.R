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
           prior = 'list',
           extra = 'list'
         ),
         prototype=list(
           mean=0,
           sigma=matrix(1),
           L=matrix(1),
           invL=matrix(1),
           logDetSigma=1,
           prior = list(),
           extra = list()
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


setMethod("getPriorDensity",
          signature=list(
            dist = "mvndist"
          ),
          definition=function(dist, log=TRUE) {

              prior <- dist@prior
              priorDens <- NA_real_
              if (prior$type == "Jeffrey") {
                 priorDens <- (-(nrow(dist@sigma) + 1)) * dist@logDetSigma 
              }
              else if (prior$type == "normal-invWish") {
                priorDens <- log(diwish(dist@sigma, prior@nu0, prior@phi)) + 
                      dmvnorm(dist@mean, prior@mu0, 1/prior@kappa0 * dist@sigma, log=TRUE)
              } else 
                  stop(paste0("Unknown prior type ", prior$type))

              if (log) priorDens else exp(priorDens)
          })


setMethod("getPosteriorSample",
          signature=list(
            dist = "mvndist",
            x = "matrix",
            n = "numeric"
          ),
          definition=function(dist, x, n) {
            
            # storage for results
            distList <- replicate(n, NULL, simplify = FALSE)
            # prior
            prior <- dist@prior
            if (prior$type == "normal-invWish") {
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
                for (i in seq(n)) {
                  # for the normal-invere Wishart prior
                  # we get a normal-inverse Wishart posterior
                  # with some new parameters mu1, kappa1, nu1, phi1
                  # sampling now works by first drawing a covariance matrix
                  # from the inverse Wishart distribution and then
                  # drawing a new mean from a multivariate normal distribution
                  # using the sampled covariance matrix
                  sampleCov <- riwish(nu1, phi1)
                  sampleMean <- as.vector(rmvnorm(1, mu1, 1/kappa1 * sampleCov))
                  distList[[i]] <- createDist_MVN(sampleMean, sampleCov, prior = prior)
                }
            }
            else if (prior$type=="Jeffrey") {
              
               # for the Jeffrey's prior, we need Gibbs sampling to alternate
               # between sampling a mean vector and a covariance matrix
               # given a mean vector, the conditional posterior for the
               # covariance matrix is an inverse Wishart distribution
               # given a covariance matrix, the mean is sampled from 
               # a multivariate normal distribution
               numObs <- ncol(x)
               obsMean <- rowMeans(x)
               p <- length(obsMean)

               sampleCov <- dist@sigma
               sampleMean <- dist@mean

               # FIXME this is a cheat...
               #if (numObs==0)
               #    obsMean <- rep(0, p)
               #numObs <- max(numObs, p)

               for (i in seq(n)) {
                   sampleMean <- as.vector(rmvnorm(1, obsMean, sampleCov / numObs))
                   Smat <- numObs * cov.wt(t(x), center = sampleMean, method = "ML")$cov
                   # FIXME another cheat
                   #diag(Smat) <- diag(Smat) + 1e-6
                   sampleCov <- riwish(numObs, Smat)
                   distList[[i]] <- createDist_MVN(sampleMean, sampleCov, prior = prior)
               }
            }
            distList
          })


