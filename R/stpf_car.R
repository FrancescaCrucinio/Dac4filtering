stpf_car <- function(xOld, obs, sigmaX, sigmaY){
  # dimension, number islands and number of particles
  d <- dim(xOld)[3]
  Nparticles <- nrow(xOld)
  M <- ncol(xOld)

  x <- array(0, dim = c(Nparticles, M, d))
  lZ <- rep(0, times = Nparticles)
  # loop over islands
  for(i in 1:Nparticles){
    # loop over dimension
    for(j in 1:d){
      # propose
      if(j == d){
        x[i, , j] <- sum(xOld[i, , d])/d +
          rowSums(x[i, , seq(length.out = j-1), drop = FALSE], dims = 2)/d +
          sqrt(sigmaX) * rnorm(M)
      } else{
        x[i, , j] <- rowSums(xOld[i, , j:d])/d +
          rowSums(x[i, , seq(length.out = j-1), drop = FALSE], dims = 2)/d +
          sqrt(sigmaX) * rnorm(M)
      }
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
