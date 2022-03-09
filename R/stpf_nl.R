stpf_nl <- function(xOld, obs, nu, sigmaX){
  # dimension, number islands and number of particles
  d <- dim(xOld)[2]
  Nparticles <- dim(xOld)[3]
  M <- dim(xOld)[4]

  x <- array(0, dim = c(d, d, Nparticles, M))
  lZ <- rep(0, times = Nparticles)
  # loop over islands
  for(i in 1:Nparticles){
    for (col in 1:d){
      for (row in 1:d){
        out_neighbours <- get_neighbours_weights(row, col, d)
        xMean <- sapply(1:M, sample_mixture, out_neighbours$mixture_weights,
                        out_neighbours$current_x_neighbours, xOld[, , i, ], simplify = TRUE)
        x[row, col, i, ] <- xMean + sqrt(sigmaX)*rnorm(M)
        # weights
        lW <- -0.5*(nu+1)*log(1+(x[row, col, i, ] - obs[row, col])^2/nu)
        max.lW <- max(lW)
        W <- exp(lW - max(lW))
        lZ[i] <- log(mean(W)) + max.lW
        # resampling
        W <- W/sum(W)
        ancestors <- stratified_resample(W, M)
        xOld[, , i, ] <- xOld[, , i, ancestors]
        x[1:(row-1), col, i, ] <- x[1:(row-1), col, i, ancestors]
        x[ , 1:(col-1), i, ] <- x[ , 1:(col-1), i, ancestors]
        x[row, col, i, ] <- x[row, col, i, ancestors]
      }
    }
  }
  # resampling islands
  Wisland <- exp(lZ - max(lZ))
  Wisland <- Wisland/sum(Wisland)
  ancestors_island <- stratified_resample(Wisland, Nparticles)
  x <- x[, , ancestors_island, ]
  return(x)
}
