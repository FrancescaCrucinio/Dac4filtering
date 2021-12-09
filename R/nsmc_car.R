nsmc_car <- function(xOld, obs, sigmaX, sigmaY, M){
  # dimension, number islands and number of particles
  d <- ncol(xOld)
  Nparticles <- nrow(xOld)

  x <- xOld
  lZ <- rep(0, times = Nparticles)
  # loop over outer level particles
  q <- vector(mode = "list", length = Nparticles)
  for(i in 1:Nparticles){
    q[[i]] <- nsmc_inner_car(x[i, ], obs, sigmaX, sigmaY, M)
    lZ[i] <- q[[i]]$lZ
  }
  # resampling outer particles
  Wouter <- exp(lZ - max(lZ))
  Wouter <- Wouter/sum(Wouter)
  ancestors_outer <- stratified_resample(Wouter, Nparticles)
  q <- q[ancestors_outer]
  # outer particles
  x <- matrix(0, nrow = Nparticles, ncol = d)
  # backward simulation
  for(i in 1:Nparticles){
    x[i, ] <- q[[i]]$xa[sample.int(M, 1), ]
  }
  return(x)
}


nsmc_inner_car <- function(xinnerOld, obs, sigmaX, sigmaY, M){
  xinner <- matrix(0, nrow = M, ncol = d)
  xinnera <- matrix(0, nrow = M, ncol = d)
  # loop over dimension
  # j=1
  j <- 1
  xinner[, j] <- sum(xinnerOld) + sqrt(sigmaX) * rnorm(M)
  # weights
  lW <- -0.5*(obs[j] - xinner[, j])^2/sigmaY - 0.5*log(2*pi*sigmaY)
  max.lW <- max(lW)
  W <- exp(lW - max(lW))
  lZ <- log(mean(W)) + max.lW
  # resampling
  W <- W/sum(W)
  ancestors <- stratified_resample(W, M)
  xinner[, j] <- xinner[ancestors, j]
  xinnera[, j] <- xinner[, j]
  for(j in 2:d){
    # propose
    xinner[, j] <- (sum(xinnerOld[j:d]) + sum(xinner[, 1:j-1]))/d + sqrt(sigmaX) * rnorm(M)
    # weights
    lW <- -0.5*(obs[j] - xinner[, j])^2/sigmaY - 0.5*log(2*pi*sigmaY)
    max.lW <- max(lW)
    W <- exp(lW - max(lW))
    lZ <- lZ + log(mean(W)) + max.lW
    # resampling
    W <- W/sum(W)
    ancestors <- stratified_resample(W, M)
    xinner[, j] <- xinner[ancestors, j]
    xinnera[,1:j] <- xinnera[ancestors, 1:j]
    xinnera[,j] <- xinner[, j]
  }
  return(list("lZ" = lZ, "x" = xinner, "xa" = xinnera))
}
