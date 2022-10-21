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
  nv_last <- 0
  x <- matrix(0, nrow = Nparticles, ncol = d)
  lW <- matrix(0, nrow = Nparticles, ncol = d)
  W <- matrix(0, nrow = Nparticles, ncol = d)
  for (i in 1:p){
    sample_from_past <- sample.int(Nparticles, Nparticles, replace = TRUE)
    x[, i] <- mu[i] + Phi[i]*(history[sample_from_past, i] - mu[i]) + Sigma[i]*rnorm(Nparticles)
    lW[, i] <- -0.5*(obs[i])^2/exp(x[i, ])
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

  # loop over tree levels excluding leaves
  for (u in 1:nlevels){
    # number of nodes at this level
    nodes <- ifelse(u<2, nchild^(nlevels-u)- nchild + 1, nchild^(nlevels-u))
    # number of variables in each node
    nvNew <- nchild^u
    nvNew_last <- nchild^u - nchild
    # updated particles/normalizing constant
    xNew <- matrix(0, nrow = Nparticles, ncol = d)

    for (i in 1:nodes){
      # get children indices
      ci <- child_indices(i, nvNew)
      # lightweight mixture resampling
      if(adaptive){
        out <- marginal_msv_light()
      } else {
        out <- marginal_msv_light_fixed_theta()
      }
      # update after mixture resampling
      indices <- out$resampled_indices
      xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv-1)], x[indices[, 2], (ci[1]+nv):ci[2]])
    }
    x <- xNew
    nv <- nvNew
    nv_last <- nvNew_last
  }
  return(x)
}
