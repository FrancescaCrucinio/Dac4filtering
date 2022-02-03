stpf_nl <- function(xOld, obs, eta){
  # dimension, number islands and number of particles
  d <- dim(xOld)[3]
  Nparticles <- nrow(xOld)
  M <- ncol(xOld)

  x <- array(0, dim = c(Nparticles, M, d, d))
  lZ <- rep(0, times = Nparticles)
  # loop over islands
  for(i in 1:Nparticles){
    for (col in 1:d){
      for (row in 1:d){
        mixture_weights <- rep(0, times = 5)
        mixture_weights[1] <- 1
        if (row > 1) mixture_weights[3] <- 0.5
        if (row < d-1) mixture_weights[5] <- 0.5
        if (col > 1) mixture_weights[2] <- 0.5
        if (col < d-1) mixture_weights[4] <- 0.5
        mixture_weights <- mixture_weights/sum(mixture_weights)
        xMean <- rep(0, times = M)
        for (j in 1:M){
          # sample component of mixture
          mixture_component <- sample.int(5, size = 1, prob = mixture_weights)
          switch(mixture_component,
               1={xMean[j] <- xOld[i, j, row, col]},
               2={xMean[j] <- xOld[i, j, row, col-1]},
               3={xMean[j] <- xOld[i, j, row-1, col]},
               4={xMean[j] <- xOld[i, j, row, col+1]},
               5={xMean[j] <- xOld[i, j+1, row, col]}
          )
        }
        x[i, , row, col] <- xMean + rnorm(M)
        # weights
        lW <- -0.5*(eta+1)*log(1+(x[i, row, col] - obs[row, col])^2/eta)
        max.lW <- max(lW)
        W <- exp(lW - max(lW))
        lZ[i] <- log(mean(W)) + max.lW
        # resampling
        W <- W/sum(W)
        ancestors <- stratified_resample(W, M)
        xOld[i, , , ] <- xOld[i, ancestors, , ]
        x[i, , 1:(row-1), col] <- x[i, ancestors, 1:(row-1), col]
        x[i, , , 1:(col-1)] <- x[i, ancestors, , 1:(col-1)]
        x[i, , row, col] <- x[i, ancestors, row, col]
      }
    }
  }
  # resampling islands
  Wisland <- exp(lZ - max(lZ))
  Wisland <- Wisland/sum(Wisland)
  ancestors_island <- stratified_resample(Wisland, Nparticles)
  x <- x[ancestors_island, , ,]
  return(x)
}
