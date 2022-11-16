marginal_msv_light_fixed_theta_first_step <- function(i, u, nv, ci, W, Nparticles, theta, SigmaV, Sigma0, x, obs){
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
    obs_covariance <- diag(sqrt(exp(mx[ci[1]:ci[2]]))) %*% SigmaV[ci[1]:ci[2], ci[1]:ci[2]] %*% diag(sqrt(exp(mx[ci[1]:ci[2]])))
    obs_node <- -0.5*obs[ci[1]:ci[2]] %*% solve(obs_covariance) %*% obs[ci[1]:ci[2]]
    obs_left <- -0.5*obs[ci[1]:(ci[1]+nv-1)] %*% solve(obs_covariance[1:nv, 1:nv, drop = FALSE]) %*%
      obs[ci[1]:(ci[1]+nv-1)]
    obs_right <- -0.5*obs[(ci[1]+nv):ci[2]] %*% solve(obs_covariance[(nv+1):(2*nv), (nv+1):(2*nv), drop = FALSE]) %*%
      obs[(ci[1]+nv):ci[2]]
    lWmix[n] <- transition_node - transition_left - transition_right + obs_node - obs_left - obs_right
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}


marginal_msv_light_fixed_theta_first_step_vectorized <- function(i, u, nv, ci, W, Nparticles, theta, SigmaV, Sigma0, x, obs){
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
  mx <- x[indices1, ]
  mx[, (ci[1]+nv):ci[2]] <- x[indices2, (ci[1]+nv):ci[2]]
  # contribution of f_{t, u}
  transition_node <- apply(mx, 1, function(x) {-0.5*x[ci[1]:ci[2]] %*% solve(Sigma0[ci[1]:ci[2], ci[1]:ci[2]]) %*% x[ci[1]:ci[2]]})
  transition_left <- apply(mx, 1, function(x) {-0.5*x[ci[1]:(ci[1]+nv-1)] %*% solve(Sigma0[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1)]) %*% x[ci[1]:(ci[1]+nv-1)]})
  transition_right <- apply(mx, 1, function(x) {-0.5*x[(ci[1]+nv):ci[2]] %*% solve(Sigma0[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2]]) %*% x[(ci[1]+nv):ci[2]]})
  # contribution of g_{t,u}
  obs_covariance <- array(apply(mx, 1, function(x) {diag(sqrt(exp(x[ci[1]:ci[2]]))) %*% SigmaV[ci[1]:ci[2], ci[1]:ci[2]] %*% diag(sqrt(exp(x[ci[1]:ci[2]])))}), dim = c(2*nv, 2*nv, theta*Nparticles))
  obs_node <- apply(obs_covariance, 3, function(x) {-0.5*obs[ci[1]:ci[2]] %*% solve(x) %*% obs[ci[1]:ci[2]]})
  obs_left <- apply(obs_covariance, 3, function(x) {-0.5*obs[ci[1]:(ci[1]+nv-1)] %*% solve(x[1:nv, 1:nv, drop = FALSE]) %*% obs[ci[1]:(ci[1]+nv-1)]})
  obs_right <- apply(obs_covariance, 3, function(x) {-0.5*obs[(ci[1]+nv):ci[2]] %*% solve(x[(nv+1):(2*nv), (nv+1):(2*nv), drop = FALSE]) %*% obs[(ci[1]+nv):ci[2]]})
  lWmix <- transition_node - transition_left - transition_right + obs_node - obs_left - obs_right
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}

marginal_msv_light_fixed_theta <- function(i, u, nv, ci, W, Nparticles, theta, SigmaX, SigmaV, SigmaUV, phi,
                                                    x, obs_current, obs_past, history){
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
  # precompute mean vector
  mu <- t(apply(history, 1, function(x){ phi*x + SigmaUV %*% solve(SigmaV) %*% diag(1/sqrt(exp(x))) %*% obs_past}))
  for (n in 1:(theta*Nparticles)) {
    # contribution of f_{t, u}
    mx <- x[indices1[n], ]
    mx[(ci[1]+nv):ci[2]] <- x[indices2[n], (ci[1]+nv):ci[2]]
    centred_mx <- -sweep(mu, 2, mx, '-')
    # # contribution of f_{t, u}
    transition_node <- -0.5*mean(apply(centred_mx[, ci[1]:ci[2]], 1, function(x) {x %*% solve(SigmaX[ci[1]:ci[2], ci[1]:ci[2]]) %*% x }))
    transition_left <- -0.5*mean(apply(centred_mx[, ci[1]:(ci[1]+nv-1), drop = F], 1, function(x) {x %*% solve(SigmaX[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1)]) %*% x }))
    transition_right <- -0.5*mean(apply(centred_mx[, (ci[1]+nv):ci[2], drop = F], 1, function(x) {x %*% solve(SigmaX[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2]]) %*% x }))
    # contribution of g_{t,u}
    obs_covariance <- diag(sqrt(exp(mx))) %*% SigmaV %*% diag(sqrt(exp(mx)))
    obs_node <- -0.5*obs_current[ci[1]:ci[2]] %*% solve(obs_covariance[ci[1]:ci[2], ci[1]:ci[2]]) %*% obs_current[ci[1]:ci[2]]
    obs_left <- -0.5*obs_current[ci[1]:(ci[1]+nv-1)] %*% solve(obs_covariance[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1), drop = FALSE]) %*%
      obs_current[ci[1]:(ci[1]+nv-1)]
    obs_right <- -0.5*obs_current[(ci[1]+nv):ci[2]] %*% solve(obs_covariance[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2], drop = FALSE]) %*%
      obs_current[(ci[1]+nv):ci[2]]
    lWmix[n] <- transition_node - transition_left - transition_right + obs_node - obs_left - obs_right
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}

marginal_msv_light_fixed_theta_vectorized <- function(i, u, nv, ci, W, Nparticles, theta, SigmaX, SigmaV, SigmaUV, phi,
                                           x, obs_current, obs_past, history){
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
  # precompute mean vector
  mu <- t(apply(history, 1, function(x){ phi*x + SigmaUV %*% solve(SigmaV) %*% diag(1/sqrt(exp(x))) %*% obs_past}))
  for (n in 1:(theta*Nparticles)) {
    # contribution of f_{t, u}
    mx <- x[indices1[n], ]
    mx[(ci[1]+nv):ci[2]] <- x[indices2[n], (ci[1]+nv):ci[2]]
    centred_mx <- -sweep(mu, 2, mx, '-')
    # # contribution of f_{t, u}
    transition_node <- -0.5*mean(apply(centred_mx[, ci[1]:ci[2]], 1, function(x) {x %*% solve(SigmaX[ci[1]:ci[2], ci[1]:ci[2]]) %*% x }))
    transition_left <- -0.5*mean(apply(centred_mx[, ci[1]:(ci[1]+nv-1), drop = F], 1, function(x) {x %*% solve(SigmaX[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1)]) %*% x }))
    transition_right <- -0.5*mean(apply(centred_mx[, (ci[1]+nv):ci[2], drop = F], 1, function(x) {x %*% solve(SigmaX[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2]]) %*% x }))
    lWmix[n] <- transition_node - transition_left - transition_right
  }
  # contribution of g_{t,u}
  mx <- x[indices1, ]
  mx[, (ci[1]+nv):ci[2]] <- x[indices2, (ci[1]+nv):ci[2]]
  obs_covariance <- array(apply(mx, 1, function(x) {diag(sqrt(exp(x))) %*% SigmaV %*% diag(sqrt(exp(x)))}), dim = c(d, d, theta*Nparticles))
  obs_node <- apply(obs_covariance, 3, function(x) {-0.5*obs_current[ci[1]:ci[2]] %*% solve(x[ci[1]:ci[2], ci[1]:ci[2]]) %*% obs_current[ci[1]:ci[2]]})
  obs_left <- apply(obs_covariance, 3, function(x) {-0.5*obs_current[ci[1]:(ci[1]+nv-1)] %*% solve(x[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1)]) %*% obs_current[ci[1]:(ci[1]+nv-1)]})
  obs_right <- apply(obs_covariance, 3, function(x) {-0.5*obs_current[(ci[1]+nv):ci[2]] %*% solve(x[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2]]) %*% obs_current[(ci[1]+nv):ci[2]]})
  lWmix <- lWmix + obs_node - obs_left - obs_right
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}
