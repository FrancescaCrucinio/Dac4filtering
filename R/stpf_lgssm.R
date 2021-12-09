stpf_lgssm <- function(xOld, obs, tau, lambda, sigmaY){
  # dimension, number islands and number of particles
  d <- dim(xOld)[3]
  Nparticles <- nrow(xOld)
  M <- ncol(xOld)

  x <- array(0, dim = c(Nparticles, M, d))
  lZ <- rep(0, times = Nparticles)
  # loop over islands
  for(i in 1:Nparticles){
    # loop over dimension
    # j=1
    j <- 1
    x[i, , j] <- 0.5*xOld[i, , j] + rnorm(M)/sqrt(tau)
    # weights
    lW <- -0.5*(obs[j] - x[i, , j])^2/sigmaY - 0.5*log(2*pi*sigmaY)
    max.lW <- max(lW)
    W <- exp(lW - max(lW))
    lZ[i] <- log(mean(W)) + max.lW
    # resampling
    W <- W/sum(W)
    ancestors <- mult_resample(W, M)
    xOld[i, , ] <- xOld[i, ancestors, ]
    x[i, , j] <- x[i, ancestors, j]
    for(j in 2:d){
      # propose
      x[i, , j] <- (0.5*tau*xOld[i, , j] + lambda*x[i, , j-1])/(tau+lambda) + rnorm(M)/sqrt(tau+lambda)
      # weights
      lW <- -0.5*(obs[j] - x[i, , j])^2/sigmaY - 0.5*log(2*pi*sigmaY)
      max.lW <- max(lW)
      W <- exp(lW - max(lW))
      lZ[i] <- lZ[i] + log(mean(W)) + max.lW
      # resampling
      W <- W/sum(W)
      ancestors <- mult_resample(W, M)
      xOld[i, , ] <- xOld[i, ancestors, ]
      x[i, , 1:j] <- x[i, ancestors, 1:j]
    }
  }
  # resampling islands
  Wisland <- exp(lZ - max(lZ))
  Wisland <- Wisland/sum(Wisland)
  ancestors_island <- mult_resample(Wisland, Nparticles)
  x <- x[ancestors_island, ,]
  return(x)
}
