nsmc_nl <- function(xOld, obs, nu, sigmaX, M){
  # dimension, number islands and number of particles
  # dimension and number of particles
  d <- nrow(xOld)
  Nparticles <- dim(xOld)[3]

  x <- xOld
  lZ <- rep(0, times = Nparticles)
  # loop over outer level particles
  q <- vector(mode = "list", length = Nparticles)
  for(i in 1:Nparticles){
    q[[i]] <- nsmc_inner_nl(x[, , i], obs, nu, sigmaX, M)
    lZ[i] <- q[[i]]$lZ
  }
  # resampling outer particles
  Wouter <- exp(lZ - max(lZ))
  Wouter <- Wouter/sum(Wouter)
  ancestors_outer <- stratified_resample(Wouter, Nparticles)
  q <- q[ancestors_outer]
  # outer particles
  x <- array(0, dim = c(d, d, Nparticles))
  # backward simulation
  for(i in 1:Nparticles){
    for (row in 1:d) {
      for (col in 1:d) {
        x[row, col, i] <- q[[i]]$x[row, col, sample.int(M, 1)]
      }
    }
  }
  return(x)
}




nsmc_inner_nl <- function(xinnerOld, obs, nu, sigmaX, M){
  xinner <- array(0, dim = c(d, d, M))
  lZ <- 0
  for(row in 1:d){
    for (col in 1:d){
      out_neighbours <- get_neighbours_weights(row, col, d)
      xMean <- sapply(rep(1, times = M), sample_mixture, out_neighbours$mixture_weights,
                      out_neighbours$current_x_neighbours, array(xinnerOld, dim = c(dim(xinnerOld), 1)), simplify = TRUE)
      xinner[row, col, ] <- xMean + sqrt(sigmaX)*rnorm(M)
      lW <- -0.5*(nu+1)*log(1+(xinner[row, col, ] - obs[row, col])^2/nu)
      max.lW <- max(lW)
      W <- exp(lW - max(lW))
      lZ <- lZ + log(mean(W)) + max.lW
      W <- W/sum(W)
      ancestors <- stratified_resample(W, M)
      xinner[row, col, ] <- xinner[row, col, ancestors]
    }
  }
  return(list("lZ" = lZ, "x" = xinner))
}
