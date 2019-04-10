#' Create Multivariate Exponential Distribution
#'
#' @return a distribution object of an exponential distribution
#' @export
#'
createDist_exp <- function(coefs, box, prior=list()) {
  new("expdist", coefs=coefs, box=box, prior=prior)
}


setClass("expdist",
         slots=list(
           coefs = 'numeric',
           box = 'numeric',
           prior = 'list',
           extra = 'list'
        ),
         prototype=list(
           coefs = numeric(0),
           box = numeric(0),
           extra = list()
         ),
         validity = function(object) {
           length(object@box) %% 2 == 0 &&
           2*length(object@coefs) == length(object@box)
         },
         contains="probdist")



setMethod("getSample",
          signature=list(
            dist = "expdist",
            n = "numeric"
          ),
          definition=function(dist, n) {
              a <- dist@coefs
              # TODO: catastrophic cancelling because 0 / 0 if a ~ 0 but who cares...
              a[abs(a) < 1e-8] <- 1e-8
              numDim <- length(dist@box) / 2
              xlo <- dist@box[1:numDim]
              xhi <- dist@box[(numDim+1):(2*numDim)]

              p <- runif(n*length(a))
              res <- log(p*(exp(a*xhi)-exp(a*xlo)) + exp(a*xlo)) / a
              dim(res) <- c(length(a), n) 
              res 
          })


setMethod("getDensity",
          signature=list(
            dist = "expdist",
            x = "matrix"
          ),
          probdens <- function(dist, x, log=TRUE) {
              a <- dist@coefs
              # TODO: catastrophic cancelling because 0 / 0 if a ~ 0 but who cares...
              a[abs(a) < 1e-8] <- 1e-8
              numDim <- length(dist@box) / 2
              xlo <- dist@box[1:numDim]
              xhi <- dist@box[(numDim+1):(2*numDim)]
              stopifnot(nrow(x) == numDim)
              const <- sum(log(abs(a))) - sum(log(abs(exp((a*xhi)) - exp(a*xlo)))) 
              res <- const + colSums(a*x) 
              # filter points outside
              sel <- apply(x >= xlo & x <= xhi, 2, all) 
              res[!sel] <- (-Inf)
              if (log) res else exp(res)
          })


setMethod("getMLEstimate",
          signature=list(
            dist = "expdist",
            x = "matrix"
          ),
          definition=function(dist, x, weights) {
              a <- dist@coefs
              numDim <- length(dist@box) / 2
              xlo <- dist@box[1:numDim]
              xhi <- dist@box[(numDim+1):(2*numDim)]
              stopifnot(nrow(x)==length(a))
              # filter points outside box
              sel <- apply(x >= xlo & x <= xhi, 2, all)
              x <- x[, sel]
              a_mle <- rowMeans(x)
              createDist_exp(a_mle, dist@box)
          })


setMethod("getPriorDensity",
          signature=list(
            dist = "expdist"
          ),
          definition=function(dist, log=TRUE) {
              # uniform prior for the moment
              prior <- dist@prior
              alo <- prior$coefs_lo
              ahi <- prior$coefs_hi
              res <- (-1)*sum(log(abs(alo-ahi)))
              if (log) res else exp(res)  
          })


setMethod("getPosteriorSample",
          signature=list(
            dist = "expdist",
            x = "matrix",
            n = "numeric"
          ),
          definition=function(dist, x, n) {

              prior <- dist@prior
              alo <- prior$coefs_lo
              ahi <- prior$coefs_hi
              numCoefs <- length(alo)
              numDim <- length(dist@box) / 2
              xlo <- dist@box[1:numDim]
              xhi <- dist@box[(numDim+1):(2*numDim)]

              xmeans <- rowMeans(x)
              N <- ncol(x)

              agrid <- alo + (ahi-alo) * t(matrix(seq(0, 1, length=1e4), 1e4, numCoefs))
              # TODO: catastrophic cancelling because 0 / 0 if a ~ 0 but who cares...
              agrid[abs(agrid) < 1e-8] <- 1e-8

              const <- N*log(abs(agrid)) - N*log(abs(exp((agrid*xhi)) - exp(agrid*xlo))) 
              logdens <- const + N*agrid*xmeans
              logdens <- t(apply(logdens, 1, function(x) { x-max(x) } ))
              dens <- exp(logdens)
              cdf <- t(apply(dens, 1, function(x) { x[-length(x)] + diff(x) / 2 }))
              diffa <- t(apply(agrid, 1, diff))
              cdf <- cdf * diffa
              cdf <- t(apply(cdf, 1, cumsum))
              cdf <- cdf / cdf[,ncol(cdf)]

              res <- matrix(0, numCoefs, n)
              for (i in seq_len(numCoefs)) {
                  res[i,] <- agrid[i, findInterval(runif(n), cdf[i,]) + 1]
              } 
              apply(res, 2, function(a) {
                  curDist <- dist
                  curDist@coefs <- a
                  curDist
              })
          }) 


