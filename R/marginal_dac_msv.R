# Lightweight resampling DaC for linear Gaussian SSM
marginal_dac_msv <- function(history, obs, mu, Phi, Lambda, Sigma, adaptive = FALSE){
  # observation dimension
  p <- length(obs)
  # number of factors
  r <- 4
  # latent space dimension
  d <- p+r
  # number of samples for lightweight mixture
  theta <- ceiling(sqrt(Nparticles))
  # number of particles
  Nparticles <- nrow(history)
  # tree topology
  nchild <- 2
  nlevels <- ceil(log2(d))
  # leaves
  # number of variables in each node
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  lW <- matrix(0, nrow = Nparticles, ncol = d)
  W <- matrix(0, nrow = Nparticles, ncol = d)
  for (i in 1:p){
    sample_from_past <- sample.int(Nparticles, Nparticles, replace = TRUE)
    x[, i] <- mu[i] + Phi[i]*(history[sample_from_past, i] - mu[i]) + Sigma[i]*rnorm(Nparticles)
    lW[, i] <- -0.5*(obs[i])^2/exp(x[, i])
    max.lW <- max(lW[, i])
    W[, i] <- exp(lW[, i] - max.lW)
    W[, i] <- W[, i]/sum(W[, i])
  }
  for (i in (p+1):d){
    sample_from_past <- sample.int(Nparticles, Nparticles, replace = TRUE)
    x[, i] <- mu[i] + Phi[i]*(history[sample_from_past, i] - mu[i]) + Sigma[i]*rnorm(Nparticles)
    lW[, i] <- 0
    W[, i] <- 1
  }

  # loop over first two levels of tree excluding leaves
  for (u in 1:2){
    # number of nodes at this level
    nodes <- ifelse(u<2, nchild^(nlevels-u)- nchild + 1, nchild^(nlevels-u))
    # number of variables in each node
    nvNew <- nchild^u
    # node iteration starting point
    node_start <- ifelse(u<2, 1, 2)
    # updated particles
    xNew <- matrix(0, nrow = Nparticles, ncol = d)
    if(u == 2){
      xNew[, 1:2] <- x[, 1:2]
    }
    for (i in node_start:nodes){
      # get children indices
      if(u == 2){
        ci <- child_indices(i, nvNew) - 2
      } else {
        ci <- child_indices(i, nvNew)
      }
      # lightweight mixture resampling
      if(adaptive){
        out <- marginal_msv_light()
      } else {
        out <- marginal_msv_light_fixed_theta_nofactor(i, u, nv, ci, W, Nparticles, theta, mu, Phi, Sigma, x, history)
      }
      # update after mixture resampling
      indices <- out$resampled_indices
      xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv-1)], x[indices[, 2], (ci[1]+nv):ci[2]])
    }
    x <- xNew
    nv <- nvNew
  }
  # number of indices in first node
  nv_first <- 2
  # loop over remaining levels of tree
  for (u in 3:4){
    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    # number of variables in each node
    nvNew <- nchild^u
    # updated particles
    xNew <- matrix(0, nrow = Nparticles, ncol = d)

    for (i in seq_len(length.out = nodes-1)){ # first nodes have no factor contribution
      # get children indices
      if(i == 1){
        ci <- child_indices(i, nvNew-2)

      } else {
        ci <- child_indices(i, nvNew) - 2
      }
      # lightweight mixture resampling
      if(adaptive){
        out <- marginal_msv_light()
      } else {
        if(i == 1){
          out <- marginal_msv_light_fixed_theta_nofactor(i, u, nv_first, ci, W, Nparticles, theta, mu, Phi, Sigma, x, history)
        } else {
          out <- marginal_msv_light_fixed_theta_nofactor(i, u, nv, ci, W, Nparticles, theta, mu, Phi, Sigma, x, history)
        }
      }
      # update after mixture resampling
      indices <- out$resampled_indices
      if(i == 1){
        xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv_first-1)], x[indices[, 2], (ci[1]+nv_first):ci[2]])
        nv_first <- nvNew-2
      } else {
        xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv-1)], x[indices[, 2], (ci[1]+nv):ci[2]])
      }
    }
    # last node in level
    i <- nodes
    ci <- child_indices(i, nvNew) - 2
    # lightweight mixture resampling
    if(adaptive){
      out <- marginal_msv_light()
    } else {
      out <- marginal_msv_light_fixed_theta_factor(i, u, nv, ci, W, Nparticles, theta, mu, Phi, Sigma, Lambda, x, history, obs)
    }
    # update after mixture resampling
    indices <- out$resampled_indices
    xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv-1)], x[indices[, 2], (ci[1]+nv):ci[2]])
    x <- xNew
    nv <- nvNew
  }
  # root level
  i <- 1
  ci <- c(1, 30)
  # updated particles
  xNew <- matrix(0, nrow = Nparticles, ncol = d)
  # lightweight mixture resampling
  if(adaptive){
    out <- marginal_msv_light()
  } else {
    out <- marginal_msv_light_fixed_theta_factor(i, u, nv_first, ci, W, Nparticles, theta, mu, Phi, Sigma, Lambda, x, history, obs)
  }
  # update after mixture resampling
  indices <- out$resampled_indices
  xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv_first-1)], x[indices[, 2], (ci[1]+nv_first):ci[2]])
  x <- xNew
  return(x)
}
