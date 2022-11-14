marginal_dac_msv_first_step <- function(obs, SigmaV, Sigma0, Nparticles, adaptive = FALSE){
  # dimension
  d <- ncol(Sigma0)
  # number of samples for lightweight mixture
  theta <- ceiling(sqrt(Nparticles))
  # tree topology
  nchild <- 2
  nlevels <- ceil(log2(d))
  # leaves
  # number of variables in each node
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  lW <- matrix(0, nrow = Nparticles, ncol = d)
  W <- matrix(0, nrow = Nparticles, ncol = d)
  for (i in 1:d){
    x[, i] <- sqrt(Sigma0[i, i])*rnorm(Nparticles)
    lW[, i] <- -0.5*(obs[i])^2/(exp(x[, i])*SigmaV[i,i])
    max.lW <- max(lW[, i])
    W[, i] <- exp(lW[, i] - max.lW)
    W[, i] <- W[, i]/sum(W[, i])
  }

  # loop over first two levels of tree excluding leaves
  for (u in 1:nlevels){
    nodes <- nchild^(nlevels-u)
    # number of variables in each node
    nvNew <- nchild^u
    # updated particles
    xNew <- matrix(0, nrow = Nparticles, ncol = d)
    for (i in 1:nodes){
      # get children indices
      ci <- child_indices(i, nvNew)
      # lightweight mixture resampling
      if(adaptive){
        out <- marginal_msv_light()
      } else {
        out <- marginal_msv_nofactor_light_fixed_theta_first_step(i, u, nv, ci, W, Nparticles, theta, SigmaV, Sigma0, x, obs)
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
