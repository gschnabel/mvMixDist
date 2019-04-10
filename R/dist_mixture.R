

#' Create mixture distribution
#'
#' @param prop vector of probabilities for the components
#' @param compList list of distributions
#'
#' @return a distribution object containing a mixture distribution
#' @export
#'
createDist_Mix <- function(prop=numeric(0),compList=list(),
                           prior = list()) {
  
  if (length(compList)==0) {
    compDim <- 0L
    prop <- 0
  }
  else {
    compDim <- nrow(getSample(compList[[1]],1L))
    prop <- prop / sum(prop)
  }
  new("mixdist",prop=prop, comp=compList, dim=compDim,
      prior = prior)
}


setClass("mixdist",
         slots=list(
           prop="numeric",
           comp="list",
           dim="numeric",
           prior = 'list',
           extra = 'list'
         ),
         prototype=list(
           prop=numeric(0),
           comp=list(),
           dim=0L,
           prior = list(),
           extra = list()
         ),
         validity=function(object) {
           isTRUE(length(object@prop)==length(object@comp) &&
                    all(diff(sapply(object@comp,function(x) nrow(getSample(x,1L))))==0))
         },
         contains = "probdist")


setMethod("getSample",
          signature=list(
            dist = "mixdist",
            n = "numeric"
          ),
          definition=function(dist, n) {
            
            numComp <- length(dist@prop)
            compIdx <- sample(numComp,n,prob=dist@prop,replace=TRUE)
            res <- matrix(0,nrow=dist@dim, ncol=n)
            for (curComp in seq(numComp)) {
              curIdx <- which(compIdx==curComp)
              if (length(curIdx)>0)
                res[,curIdx] <- getSample(dist@comp[[curComp]],length(curIdx))
            }
            res
          })

setMethod("getDensity",
          signature=list(
            dist = "mixdist",
            x = "matrix"
          ),
          definition=function(dist, x, log=TRUE, vars=NULL) {
            numComp <- length(dist@prop)
            logProp <- log(dist@prop)
            res <- matrix(0,nrow=numComp,ncol=ncol(x))
            for (curComp in seq(numComp)) {
              res[curComp,] <- logProp[curComp] + getDensity(dist@comp[[curComp]], x, log=TRUE, vars)
            }
            res <- apply(res,2,function(x) {
              xmax <- max(x)
              log(sum(exp(x - xmax))) + xmax
            })
            if (!isTRUE(log)) res <- exp(res)
            res
          })


setMethod("getPriorDensity",
          signature=list(
            dist = "mixdist"
          ),
          definition=function(dist, log=TRUE) {
            
            numComp <- length(dist@prop)
            priorDens <- log(ddirichlet(dist@prop, dist@prior$alpha))
            for (k in seq_len(numComp)) {
              priorDens <- priorDens + getPriorDensity(dist@comp[[k]], log = TRUE)
            } 
            if (log) priorDens else exp(priorDens) 
          })


setMethod("getMLEstimate",
          signature=list(
            dist = "mixdist",
            x = "matrix"
          ),
          definition=function(dist, x, weights, maxIter=50) {
            
            numComp <- length(dist@prop)
            for (curIter in seq(maxIter)) {
              # expectation
              memShip <- getMembership(dist, x, log=FALSE)
              memShip <- t(weights*t(memShip))
              dist@prop <- rowSums(memShip) / sum(weights)
              # maximization
              for (curComp in seq(numComp))
                dist@comp[[curComp]] <- getMLEstimate(dist@comp[[curComp]], x, memShip[curComp,])
            }
            dist
          })


setMethod("getMembership",
          signature=list(
            dist = "mixdist",
            x = "matrix"
          ),
          definition=function(dist, x, log=TRUE) {
            numComp <- length(dist@prop)
            logProp <- log(dist@prop)
            res <- matrix(0,nrow=numComp,ncol=ncol(x))
            for (curComp in seq(numComp))
              res[curComp,] <- logProp[curComp] + getDensity(dist@comp[[curComp]], x, log=TRUE)
            res <- apply(res,2,function(x) x-max(x))
            dim(res) <- c(numComp, ncol(x))
            
            if (isTRUE(log)) {
              res <- apply(res,2,function(x) x - log(sum(exp(x))))
            }
            else
              res <- apply(exp(res),2,function(x) x/sum(x))
            dim(res) <- c(numComp, ncol(x))
            res
          })


setMethod("getPosteriorSample",
          signature=list(
            dist = "mixdist",
            x = "matrix",
            n = "numeric"
          ),
          definition=function(dist, x, n) {
            
            alpha <- dist@prior$alpha
            numComp <- length(dist@prop)
            distList <- replicate(n, NULL, simplify = FALSE)
            for (i in seq(n)) {
              # sample membership
              memShip <- getMembership(dist, x, log=FALSE)
              numComp <- nrow(memShip)
              z <- rep(0, ncol(x))
              for (j in seq_along(z)) {
                z[j] <- sample.int(numComp, 1, prob = memShip[,j])
              }
              # sample components
              for (curComp in seq(numComp)) {
                selectedPoints <- x[,z==curComp,drop=FALSE]
                dist@comp[[curComp]] <- getPosteriorSample(dist@comp[[curComp]], 
                                                           selectedPoints, 1)[[1]]
              }
              # sample from dirichlet posterior
              counts <- colSums(outer(z, 1:numComp, `==`))
              dist@prop <- as.vector(rdirichlet(1, alpha + counts))
              # density for marginal likelihood
              likeDens <- sum(getDensity(dist, x, log=TRUE))
              priorDens <- getPriorDensity(dist, log=TRUE)
              scaledPostDens <- likeDens + priorDens
              dist@extra <- list(priorDens = priorDens,
                                 likeDens = likeDens,
                                 scaledPostDens = scaledPostDens)
              # save
              distList[[i]] <- dist

            }
            distList
          })
