marginal_msv_nofactor_light_fixed_theta <- function(i, u, nv, ci, W, Nparticles, theta, mu, Phi, Sigma, Lambda, x, history, obs){
  nchild <- 2
  # resample on each children
  if(u == 1){
    # child 1
    indices1 <- stratified_resample(W[, nchild*(i-1)+1], theta*Nparticles)
    # child 2 (with random permutation)
    indices2 <- sample(stratified_resample(W[, nchild*i], theta*Nparticles))
  }
  else{
    # child 1
    indices1 <- sample.int(Nparticles, size = theta*Nparticles, replace = TRUE)
    # child 2
    indices2 <- sample.int(Nparticles, size = theta*Nparticles, replace = TRUE)
  }

  # mixture weights
  lWmix <- rep(0, times = theta*Nparticles)
  for (n in 1:(theta*Nparticles)) {
    # contribution of f_{t, u}
    mx <- x[indices1[n], ]
    mx[(ci[1]+nv):ci[2]] <- x[indices2[n], (ci[1]+nv):ci[2]]
    # # contribution of f_{t, u}
    # transition_node <- -0.5*mx[ci[1]:ci[2]] %*% solve(Sigma0[ci[1]:ci[2], ci[1]:ci[2]]) %*% mx[ci[1]:ci[2]]
    # transition_left <- -0.5*mx[ci[1]:(ci[1]+nv-1)] %*% solve(Sigma0[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1), drop = FALSE]) %*%
    #   mx[ci[1]:(ci[1]+nv-1)]
    # transition_right <- -0.5*mx[(ci[1]+nv):ci[2]] %*% solve(Sigma0[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2], drop = FALSE]) %*%
    #   mx[(ci[1]+nv):ci[2]]
    # contribution of g_{t,u}
    obs_covariance <- diag(sqrt(exp(mx))) %*% SigmaV %*% diag(sqrt(exp(mx)))
    obs_node <- -0.5*obs_current[ci[1]:ci[2]] %*% solve(obs_covariance[ci[1]:ci[2], ci[1]:ci[2]]) %*% obs_current[ci[1]:ci[2]]
    obs_left <- -0.5*obs_current[ci[1]:(ci[1]+nv-1)] %*% solve(obs_covariance[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1), drop = FALSE]) %*%
      obs_current[ci[1]:(ci[1]+nv-1)]
    obs_right <- -0.5*obs_current[(ci[1]+nv):ci[2]] %*% solve(obs_covariance[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2], drop = FALSE]) %*%
      obs_current[(ci[1]+nv):ci[2]]

    lWmix <- obs_node - obs_left - obs_right
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}

marginal_msv_nofactor_light_fixed_theta_first_step <- function(i, u, nv, ci, W, Nparticles, theta, SigmaV, Sigma0, x, obs){
  nchild <- 2
  # resample on each children
  if(u == 1){
    # child 1
    indices1 <- stratified_resample(W[, nchild*(i-1)+1], theta*Nparticles)
    # child 2 (with random permutation)
    indices2 <- sample(stratified_resample(W[, nchild*i], theta*Nparticles))
  }
  else{
    # child 1
    indices1 <- sample.int(Nparticles, size = theta*Nparticles, replace = TRUE)
    # child 2
    indices2 <- sample.int(Nparticles, size = theta*Nparticles, replace = TRUE)
  }

  # mixture weights
  lWmix <- rep(0, times = theta*Nparticles)
  for (n in 1:(theta*Nparticles)) {
    mx <- x[indices1[n], ]
    mx[(ci[1]+nv):ci[2]] <- x[indices2[n], (ci[1]+nv):ci[2]]
    # contribution of f_{t, u}
    transition_node <- -0.5*mx[ci[1]:ci[2]] %*% solve(Sigma0[ci[1]:ci[2], ci[1]:ci[2]]) %*% mx[ci[1]:ci[2]]
    transition_left <- -0.5*mx[ci[1]:(ci[1]+nv-1)] %*% solve(Sigma0[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1), drop = FALSE]) %*%
      mx[ci[1]:(ci[1]+nv-1)]
    transition_right <- -0.5*mx[(ci[1]+nv):ci[2]] %*% solve(Sigma0[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2], drop = FALSE]) %*%
      mx[(ci[1]+nv):ci[2]]
    # contribution of g_{t,u}
    obs_covariance <- diag(sqrt(exp(mx))) %*% SigmaV %*% diag(sqrt(exp(mx)))
    obs_node <- -0.5*obs[ci[1]:ci[2]] %*% solve(obs_covariance[ci[1]:ci[2], ci[1]:ci[2]]) %*% obs[ci[1]:ci[2]]
    obs_left <- -0.5*obs[ci[1]:(ci[1]+nv-1)] %*% solve(obs_covariance[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1), drop = FALSE]) %*%
      obs[ci[1]:(ci[1]+nv-1)]
    obs_right <- -0.5*obs[(ci[1]+nv):ci[2]] %*% solve(obs_covariance[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2], drop = FALSE]) %*%
      obs[(ci[1]+nv):ci[2]]

    lWmix[n] <- transition_node - transition_left - transition_right + obs_node - obs_left - obs_right
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}
