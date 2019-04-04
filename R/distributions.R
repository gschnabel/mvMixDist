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
  myDist <- createDist_Mix(c(0.1, 0.3, 0.6), list(myDist1, myDist2, myDist3),
                           prior = list(alpha = rep(1,3)))
  
  smpl <- getSample(myDist, 1000)
  getPriorDensity(myDist, log=TRUE)
  

  getMLEstimate(myDist, smpl)
  myDist@comp[[1]]@mean <- c(-10,30)
  myDist@comp[[1]]@sigma <- diag(c(1000,1000))
  myDist@comp[[2]]@mean <- c(15,35)
  myDist@comp[[2]]@sigma <- diag(c(1000,1000))
  postSmpl <- getPosteriorSample(myDist, smpl, 4000)
  
  y <- sapply(postSmpl, function(x) x@prop[3])
  y <- sapply(postSmpl, function(x) x@comp[[1]]@sigma[2,2])
  y <- sapply(postSmpl, function(x) x@comp[[2]]@mean[2])
  x <- seq_along(y)
  plot(x,y,type="l")
  mean(y)
  sd(y)
  hist(y)

  library(ellipse)
  for (i in 1:1000) {

    plot(smpl[1,], smpl[2,])
    tmp <- postSmpl[[i]]@comp[[2]]
    lines(ellipse(x = diag(c(5,30)), centre = c(5,5), level = 0.68), col = 'green', lwd = 2)
    lines(ellipse(x = diag(c(5,30)), centre = c(5,5), level = 0.95), col = 'green', lwd = 2)
    lines(ellipse(x = tmp@sigma, centre = tmp@mean, level = 0.68), col = 'red')
    lines(ellipse(x = tmp@sigma, centre = tmp@mean, level = 0.95), col = 'red')

    tmp <- postSmpl[[i]]@comp[[1]]
    lines(ellipse(x = diag(c(10,40)), centre = c(20,15), level = 0.68), col = 'green', lwd = 2)
    lines(ellipse(x = diag(c(10,40)), centre = c(20,15), level = 0.95), col = 'green', lwd = 2)
    lines(ellipse(x = tmp@sigma, centre = tmp@mean, level = 0.68), col = 'red')
    lines(ellipse(x = tmp@sigma, centre = tmp@mean, level = 0.95), col = 'red')
    Sys.sleep(0.2)
  }


}



