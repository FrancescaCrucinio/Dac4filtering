dac_nl_lightweight <- function(history, obs, sigmaX, nu, M = NULL){
  if(is.null(M)) {
    # number of samples for lightweight mixture (no adaptation)
    M <- ceiling(sqrt(Nparticles))
  }
  # dimension and number of particles
  d <- ncol(history[, , , 1])
  Nparticles <- nrow(history[, , , 1])
  # tree topology
  nchild <- 2
  nlevels <- log2(d^2)
  # leaves
  # number of variables
  nv <- 1
  x <- array(0, dim = c(Nparticles, d, d))
  lW <- array(0, dim = c(Nparticles, d, d))
  W <- array(0, dim = c(Nparticles, d, d))
  # history indices
  historyIndex <- array(1:Nparticles, dim = c(Nparticles, d^2, d))
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
      xMean <- sapply(1:Nparticles, sample_mixture, history[, , , 1], row, col, mixture_weights, simplify = TRUE)
      x[, row, col] <- xMean + sqrt(sigmaX)*rnorm(Nparticles)
      # weights
      lW[, row, col] <- -0.5*(nu+d)*log(1+(x[, row, col] - obs[row, col])^2/nu)
      max.lW <- max(lW[, row, col])
      W[, row, col] <- exp(lW[, row, col] - max.lW)
      W[, row, col] <- W[, row, col]/sum(W[, row, col])
    }
  }
  x <- matrix(x, nrow = Nparticles, ncol = d^2)
  lW <- matrix(lW, nrow = Nparticles, ncol = d^2)
  W <- matrix(W, nrow = Nparticles, ncol = d^2)
  obs <- c(t(obs))
  # loop over tree levels excluding leaves
  for (u in 1:nlevels){

    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    # number of variables in each node
    nvNew <- nchild^u

    # updated particles
    xNew <- matrix(0, nrow = Nparticles, ncol = d^2)
    # updated history
    # historyIndexNew <- array(0, dim = c(Nparticles, d^2, nodes))
    for (i in 1:nodes) {
      print(paste("u = ", u, "i = ", i))
      # get children indices
      ci <- child_indices_lattice(d, u, i, nv, nodes)
      print(paste(c(ci)))
      # if(u > 1){
      #   # mutation
      #
      # }
      # else{ # at the leaf level all histories are the same
        # historyIndexNew[, , i] <- historyIndex[, , i]
      # }
      # lightweight mixture resampling
      if(M == "adaptive") {

      }
      else{
        out <- nl_light(i, u, nv, ci, W, Nparticles, M, eta, x)
        indices <- out$resampled_indices
        xNew[, c(ci)] <- cbind(x[indices[, 1], ci[, 1]], x[indices[, 2], ci[, 2]])
        # historyIndexNew[, , i] <- historyIndexNew[indices[, 1], , i]
      }
    }
    # historyIndex <- historyIndexNew
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
