dac_nl_lightweight <- function(history, obs, sigmaX, eta){
  # dimension and number of particles
  d <- ncol(history[, , , 1])
  Nparticles <- nrow(history[, , , 1])
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # leaves
  # number of variables
  nv <- 1
  x <- array(0, dim = c(Nparticles, d, d))
  lW <- array(0, dim = c(Nparticles, d, d))
  # history indices
  historyIndex <- array(1:Nparticles, dim = c(Nparticles, d, d))
  # leaves
  for (col in 1:d){
    for (row in 1:d){
      mixture_weights <- rep(0, times = 5)
      mixture_weights[1] <- 1
      if (row > 1) mixture_weights[3] <- 0.5
      if (row < d-1) mixture_weights[5] <- 0.5
      if (col > 1) mixture_weights[2] <- 0.5
      if (col < d-1) mixture_weights[4] <- 0.5
      mixture_weights <- mixture_weights/sum(mixture_weights)
    }
    xMean <- sapply(1:Nparticles, sample_mixture, history[, , , 1], row, col, mixture_weights, simplify = TRUE)
    x[, row, col] <- xMean + sqrt(sigmaX)*rnorm(Nparticles)
    # weights
    lW <- -0.5*(eta+1)*log(1+sweep(x[, row, col], 2, obs[row, col])^2/eta)
  }
  # loop over tree levels excluding leaves
  for (u in 1:nlevels){
    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    # number of variables in each node
    nvNew <- nchild^u

    # updated particles
    xNew <- matrix(0, nrow = Nparticles, ncol = d)
    # updated history
    historyIndexNew <- array(0, dim = c(Nparticles, d, nodes))
    for (i in 1:nodes) {
      # get children indices
      ci <- child_indices_lattice(d, u, i, nvNew, nv)
      if(u > 1){
        # mutation

      }
      else{ # at the leaf level all histories are the same
        historyIndexNew[, , i] <- historyIndex[, , i]
      }
      # adaptive lightweight mixture resampling

    }
    historyIndex <- historyIndexNew
    x <- xNew
    nv <- nvNew
  }
  return(x)
}

sample_mixture <- function(n, xOld, row, col, mixture_weights){
  # sample component of mixture
  mixture_component <- sample.int(5, size = 1, prob = mixture_weights)
  switch(mixture_component,
         {xMean <- xOld[n, row, col]},
         {xMean <- xOld[n, row, col-1]},
         {xMean <- xOld[n, row-1, col]},
         {xMean <- xOld[n, row, col+1]},
         {xMean <- xOld[n, row+1, col]}
  )
  return(xMean)
}
