# Lightweight resampling DaC for linear Gaussian SSM
marginal_dac_lgssm_lightweight <- function(history, obs, tau, lambda, sigmaY, adaptive = FALSE){
  # number of samples for lightweight mixture
  theta <- ceiling(sqrt(Nparticles))
  # dimension and number of particles
  d <- ncol(history)
  Nparticles <- nrow(history)
  memory <- ifelse(Nparticles >= 10^3 & d == 2048, TRUE, FALSE)
  memory <- TRUE
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # leaves
  # number of variables
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  lW <- matrix(0, nrow = Nparticles, ncol = d)
  W <- matrix(0, nrow = Nparticles, ncol = d)
  for (i in 1:nchild^nlevels){
    sample_from_past <- sample.int(Nparticles, Nparticles, replace = TRUE)
    if (i == 1) {
      x[, i] <- 0.5*history[sample_from_past, i] + rnorm(Nparticles)/sqrt(tau)
    } else {
      x[, i] <- 0.5*tau*history[sample_from_past, i]/(tau+lambda) + rnorm(Nparticles)/sqrt(tau+lambda)
    }
    lW[, i] <- -0.5*(obs[i] - x[, i])^2/sigmaY - 0.5*log(2*pi*sigmaY)
    max.lW <- max(lW[, i])
    W[, i] <- exp(lW[, i] - max.lW)
    W[, i] <- W[, i]/sum(W[, i])
  }

  # loop over tree levels excluding leaves
  for (u in 1:nlevels){
    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    # number of variables in each node
    nvNew <- nchild^u
    # updated particles/normalizing constant
    xNew <- matrix(0, nrow = Nparticles, ncol = d)

    for (i in 1:nodes){
      # get children indices
      ci <- child_indices(i, nvNew)
      # lightweight mixture resampling
      if(adaptive){
        out <- marginal_lgssm_light_adaptive(Nparticles, i, u, nv, ci, lW, Nparticles, lambda, tau, x, history, memory)
      } else {
        out <- marginal_lgssm_light(i, u, nv, ci, W, Nparticles, theta, lambda, tau, x, history, memory)
      }
      # update after mixture resampling
      indices <- out$resampled_indices
      xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv-1)], x[indices[, 2], (ci[1]+nv):ci[2]])
     }
    x <- xNew
    nv <- nvNew
  }
  return(x)
}
