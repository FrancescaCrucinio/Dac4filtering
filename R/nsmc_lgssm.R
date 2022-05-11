# Nested SMC for linear Gaussian SSM
nsmc_lgssm <- function(xOld, obs, tau, lambda, sigmaY, M){
  d <- ncol(xOld)
  Nparticles <- nrow(xOld)

  x <- xOld
  lZ <- rep(0, times = Nparticles)
  # loop over outer level particles
  q <- vector(mode = "list", length = Nparticles)
  for(i in 1:Nparticles){
    q[[i]] <- nsmc_inner_lgssm(x[i, ], obs, tau, lambda, sigmaY, M)
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
    x[i, d] <- q[[i]]$x[sample.int(M, 1), d]
    for (j in (d-1):1){
      lW <- -0.5*lambda*(x[i, j+1] - q[[i]]$x[, j])^2
      W <- exp(lW - max(lW))
      # resampling
      x[i, j] <- q[[i]]$x[stratified_resample(W/sum(W), 1), j]
    }
  }
  return(x)
}

# Inner SMC for linear Gaussian SSM
nsmc_inner_lgssm <- function(xinnerOld, obs, tau, lambda, sigmaY, M){
  xinner <- matrix(0, nrow = M, ncol = d)
  # loop over dimension
  # j=1
  j <- 1
  xinner[, j] <- 0.5*xinnerOld[j] + rnorm(M)/sqrt(tau)
  # weights
  lW <- -0.5*(obs[j] - xinner[, j])^2/sigmaY - 0.5*log(2*pi*sigmaY)
  max.lW <- max(lW)
  W <- exp(lW - max(lW))
  lZ <- log(mean(W)) + max.lW
  # resampling
  W <- W/sum(W)
  ancestors <- stratified_resample(W, M)
  xinner[ , j] <- xinner[ancestors, j]
  for(j in 2:d){
    # propose
    xinner[, j] <- (0.5*tau*xinnerOld[j] + lambda*xinner[, j-1])/(tau+lambda) + rnorm(M)/sqrt(tau+lambda)
    # weights
    lW <- -0.5*(obs[j] - xinner[, j])^2/sigmaY - 0.5*log(2*pi*sigmaY)
    max.lW <- max(lW)
    W <- exp(lW - max(lW))
    lZ <- lZ + log(mean(W)) + max.lW
    # resampling
    W <- W/sum(W)
    ancestors <- stratified_resample(W, M)
    xinner[, j] <- xinner[ancestors, j]
  }
  return(list("lZ" = lZ, "x" = xinner))
}
