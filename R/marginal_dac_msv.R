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
        out <- marginal_msv_light_fixed_theta_first_step_vectorized(i, u, nv, ci, W, Nparticles, theta, SigmaV, Sigma0, x, obs)
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

marginal_dac_msv <- function(history, obs_current, obs_past, SigmaU, SigmaV, SigmaUV, SigmaX, phi, adaptive = FALSE){
  # dimension
  d <- ncol(history)
  # number of particles
  Nparticles <- nrow(history)
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
  # precompute mean vector
  mu <- t(apply(history, 1, function(x){ phi*x + SigmaUV %*% solve(SigmaV) %*% diag(1/sqrt(exp(x))) %*% obs_past}))
  for (i in 1:d){
    sample_from_past <- sample.int(Nparticles, Nparticles, replace = TRUE)
    x[, i] <- mu[sample_from_past, i] + sqrt(SigmaX[i, i])*rnorm(Nparticles)
    lW[, i] <- -0.5*(obs_current[i])^2/(exp(x[, i])*SigmaV[i,i])
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
        out <- marginal_msv_light_fixed_theta(i, u, nv, ci, W, Nparticles, theta, SigmaX, SigmaV, SigmaUV, phi,
                                                       x, obs_current, obs_past, history)
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
