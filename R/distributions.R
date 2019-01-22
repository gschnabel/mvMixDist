dontrun <- function() {
  
  myDist1 <- createDist_MVN(mean=c(0,10),sigma=diag(c(1,100)), prior = list())
  myDist2 <- createDist_MVN(mean=c(20,10),sigma=diag(c(25,25)))
  myDist3 <- createDist_Const(c(-30,-30,60,60))
  myDist <- createDist_Mix(c(0.5,0.3,0.2),
                           list(myDist1,myDist2,myDist3))

  smpl <- getSample(myDist,100)
  plot(t(smpl),pty="s",asp=1)
  
  getMLEstimate(myDist,smpl, maxIter = 1000)
  getMembership(myDist,smpl)
  
  myDist3 <- createDist_Const(c(5,10,10,50))
  smpl <- getSample(myDist3, 100)
  plot(t(smpl))
  getDensity(myDist3, cbind(5,5), vars=1)

  myPrior <- list(mu0 = 1, kappa0 = 1,
                  nu0 = 3, phi = rbind(c(1,0), c(0,1)))
  
  myDist1 <- createDist_MVN(mean=c(0,10),sigma=diag(c(1,100)),
                            prior = myPrior)
  
  smpl <- getSample(myDist1, 1000)
  plot(t(smpl))  
  
  postSmpl <- getPosteriorSample(myDist1, smpl, 100)

  # true distribution
  myDist1 <- createDist_MVN(mean=c(0,10),sigma=diag(c(1,100)), prior = myPrior)
  myDist2 <- createDist_MVN(mean=c(20,10),sigma=diag(c(25,25)), prior = myPrior)
  myDist3 <- createDist_Const(box=c(-20,-20,40,40))
  myDist <- createDist_Mix(c(0.4, 0.4, 0.2), list(myDist1, myDist2, myDist3),
                                             prior = list(alpha = 1))
  
  smpl <- getSample(myDist, 1000)
  plot(t(smpl))
  
  # starting values for Gibbs sampling
  myDist1 <- createDist_MVN(mean=c(20,15),sigma=diag(c(10,40)), prior = myPrior)
  myDist2 <- createDist_MVN(mean=c(5,5),sigma=diag(c(5,30)), prior = myPrior)
  myDist3 <- createDist_Const(box=c(-20,-20,40,40))
  myDist <- createDist_Mix(c(0.1, 0.3, 0.6), list(myDist1, myDist2, myDist3),
                           prior = list(alpha = 1))
  
  
  getMLEstimate(myDist, smpl)
  
  postSmpl <- getPosteriorSample(myDist, smpl, 1000)
  
  y <- sapply(postSmpl, function(x) x@prop[1])
  x <- seq_along(y)
  plot(x,y,type="l")
  mean(y)
  sd(y)

}



