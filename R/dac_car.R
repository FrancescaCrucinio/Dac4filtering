dac_car_lightweight <- function(history, obs, sigmaX, sigmaY){
  # dimension and number of particles
  d <- ncol(history[, , 1])
  Nparticles <- nrow(history[, , 1])
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # leaves
  # number of variables
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  lW <- matrix(0, nrow = Nparticles, ncol = d)
  # history indices
  historyIndex <- array(1:Nparticles, dim = c(Nparticles, d, d))
  for (i in 1:nchild^nlevels){
    # propose
    x[, i] <- rowSums(history[, i:d, drop = FALSE, 1])/d + sqrt(sigmaX) * rnorm(Nparticles)
    # weights
    lW[, i] <- -0.5*(obs[i] - x[, i])^2/sigmaY - 0.5*log(2*pi*sigmaY)
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
      ci <- child_indices(i, nvNew)
      # mutation
      historyIndexNew[, , i] <- historyIndex[, , i]
      # adaptive lightweight mixture resampling
      xOld <- history[historyIndex[, , nchild*i], , 1]
      indices <- car_adaptive_light(Nparticles, i, u, nv, nvNew, ci, lW, Nparticles, sigmaX, x, xOld)
      # update particles
      xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv-1)], x[indices[, 2], (ci[1]+nv):ci[2]])
      # update xOld
      xOld[, ci[1]:ci[2]] <- xOld[indices[, 1], ci[1]:ci[2]]
      historyIndexNew[, , i] <- historyIndexNew[indices[, 1], , i]
    }
    x <- xNew
    nv <- nvNew
  }
  return(x)
}
