dontrun <- function() {

  # starting values for Gibbs sampling
  library(mvMixDist)

  #myPrior <- list(type = "normal-invWish", mu0 = 1, kappa0 = 1,
  #                nu0 = 3, phi = rbind(c(1,0), c(0,1)))
  trueMean1 <- c(20,15)
  trueMean2 <- c(5,5)
  trueSigma1 <- diag(c(10,40))
  trueSigma2 <- diag(c(5,30))

  myPrior <- list(type = "Jeffrey")
  myDist1 <- createDist_MVN(mean=trueMean1,sigma=trueSigma1, prior = myPrior)
  myDist2 <- createDist_MVN(mean=trueMean2,sigma=trueSigma2, prior = myPrior)
  myDist3 <- createDist_Const(box=c(-20,-20,40,40))
  myDist4 <- createDist_MVN(mean=trueMean1+5, sigma=trueSigma1, prior = myPrior)

  myDist <- createDist_Mix(c(0.02, 0.38, 0.6), list(myDist1, myDist2, myDist3),
                           prior = list(alpha = rep(10,3)))

  myDist_b <- createDist_Mix(c(0.4,0.6,0.3,0.3), list(myDist1,myDist2,myDist3, myDist4), prior = list(alpha=c(10,10,10,10)))
  
  smpl <- getSample(myDist, 3000)
  getPriorDensity(myDist, log=TRUE)
  

  getMLEstimate(myDist, smpl)
  myDist@comp[[1]]@mean <- c(-10,30)
  myDist@comp[[1]]@sigma <- diag(c(1000,1000))
  myDist@comp[[2]]@mean <- c(15,35)
  myDist@comp[[2]]@sigma <- diag(c(1000,1000))
  postSmpl <- getPosteriorSample(myDist, smpl, 4000)
  postSmpl <- getPosteriorSample(myDist_b, smpl, 4000)
  
  cs <- sapply(postSmpl[1:4000], function(x) x@extra$scaledPostDens)
  plot(seq_along(cs), cs)
  cs <- log(mean(exp(cs-max(cs)))) + max(cs)
  cs

  y <- sapply(postSmpl, function(x) x@prop)
  y <- sapply(postSmpl, function(x) as.vector(x@comp[[4]]@sigma))
  y <- sapply(postSmpl, function(x) as.vector(x@comp[[4]]@mean))
  x <- seq_along(y)
  plot(x,y,type="l")
  mean(y)
  sd(y)
  hist(y)

  library(ellipse)
  for (i in 1:1000) {

    plot(smpl[1,], smpl[2,])
    for (k in seq_along(myDist_b@comp)) {
      tmp <- postSmpl[[i]]@comp[[k]]
      if (! "mvndist" %in% class(tmp)) next
      lines(ellipse(x = tmp@sigma, centre = tmp@mean, level = 0.68), col = 'red')
      lines(ellipse(x = tmp@sigma, centre = tmp@mean, level = 0.95), col = 'red')
    }
    Sys.sleep(0.2)
  }


}



