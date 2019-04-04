#' Create Constant Distribution
#'
#' @return a distribution object of a constant distribution
#' @export
#'
createDist_Const <- function(box) {
  new("constdist",box=box)
}


setClass("constdist",
         slots=list(
           box = 'numeric',
           extra = 'list'
        ),
         prototype=list(
           box = numeric(0),
           extra = list()
         ),
         validity = function(object) {
           length(object@box) %% 2 == 0
         },
         contains="probdist")



setMethod("getSample",
          signature=list(
            dist = "constdist",
            n = "numeric"
          ),
          definition=function(dist, n) {
            numDim <- length(dist@box) / 2
            matrix(runif(n*numDim,
                      dist@box[1:numDim], 
                      dist@box[(numDim+1):(2*numDim)]),
                   nrow = numDim)
          })


setMethod("getDensity",
          signature=list(
            dist = "constdist",
            x = "matrix"
          ),
          definition=function(dist, x, log=TRUE, vars=NULL) {
            numDim <- length(dist@box) / 2
            mins <- dist@box[1:numDim]
            maxs <- dist@box[(numDim+1):(2*numDim)]
            lens <- abs(maxs - mins)
            if (is.null(vars)) vars <- seq_len(numDim)
            logfun <- get("log", pos = "package:base")
            if (isTRUE(log))
              -sum(logfun(lens[vars]))
            else
              1/prod(lens[vars])
          })



setMethod("getMLEstimate",
          signature=list(
            dist = "constdist",
            x = "matrix"
          ),
          definition=function(dist, x, weights) {
            createDist_Const(dist@box)
          })


setMethod("getPriorDensity",
          signature=list(
            dist = "constdist"
          ),
          definition=function(dist, log=TRUE) {
           if (log) 0 else 1  
          })


setMethod("getPosteriorSample",
          signature=list(
            dist = "constdist",
            x = "matrix",
            n = "numeric"
          ),
          definition=function(dist, x, n) {
            
            # create samples
            distList <- replicate(n, NULL, simplify = FALSE)
            for (i in seq(n)) {
              distList[[i]] <- dist
            }
            distList
          })



