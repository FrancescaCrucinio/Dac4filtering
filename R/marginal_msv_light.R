marginal_msv_light_first_step <- function(ess_target, i, u, nv, ci, lW, Nparticles, theta, SigmaV, Sigma0, x, obs){
  nchild <- 2
  # mixture weights
  lWmix <- rep(0, times = Nparticles)
  for (n in 1:Nparticles) {
    mx <- x[n, ]
    # contribution of f_{t, u}
    transition_node <- -0.5*mx[ci[1]:ci[2]] %*% solve(Sigma0[ci[1]:ci[2], ci[1]:ci[2]]) %*% mx[ci[1]:ci[2]]
    transition_left <- -0.5*mx[ci[1]:(ci[1]+nv-1)] %*% solve(Sigma0[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1), drop = FALSE]) %*% mx[ci[1]:(ci[1]+nv-1)]
    transition_right <- -0.5*mx[(ci[1]+nv):ci[2]] %*% solve(Sigma0[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2], drop = FALSE]) %*% mx[(ci[1]+nv):ci[2]]
    # contribution of g_{t,u}
    obs_covariance <- diag(sqrt(exp(mx[ci[1]:ci[2]]))) %*% SigmaV[ci[1]:ci[2], ci[1]:ci[2]] %*% diag(sqrt(exp(mx[ci[1]:ci[2]])))
    obs_node <- -0.5*obs[ci[1]:ci[2]] %*% solve(obs_covariance) %*% obs[ci[1]:ci[2]]
    obs_left <- -0.5*obs[ci[1]:(ci[1]+nv-1)] %*% solve(obs_covariance[1:nv, 1:nv, drop = FALSE]) %*% obs[ci[1]:(ci[1]+nv-1)]
    obs_right <- -0.5*obs[(ci[1]+nv):ci[2]] %*% solve(obs_covariance[(nv+1):(2*nv), (nv+1):(2*nv), drop = FALSE]) %*% obs[(ci[1]+nv):ci[2]]
    lWmix[n] <- transition_node - transition_left - transition_right + obs_node - obs_left - obs_right
  }
  if(u == 1){
    lWmix <- lWmix + lW[, (nchild*(i-1)+1)] + lW[, i*nchild]
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # build ESS
  ess_s <- sum(Wmix)
  ess_ss <- sum(Wmix^2)
  ess <- ess_s^2/ess_ss
  # first permutation
  permutation <- 1:Nparticles
  theta <- 1
  while (ess < ess_target & theta <= ceiling(sqrt(Nparticles))) {
    theta <- theta+1
    new_perm <- sample.int(Nparticles)
    # mixture weights
    lWmix_perm <- rep(0, times = Nparticles)
    for (n in 1:Nparticles) {
      mx <- x[n, ]
      mx[(ci[1]+nv):ci[2]] <- x[new_perm[n], (ci[1]+nv):ci[2]]
      # contribution of f_{t, u}
      transition_node <- -0.5*mx[ci[1]:ci[2]] %*% solve(Sigma0[ci[1]:ci[2], ci[1]:ci[2]]) %*% mx[ci[1]:ci[2]]
      transition_left <- -0.5*mx[ci[1]:(ci[1]+nv-1)] %*% solve(Sigma0[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1), drop = FALSE]) %*% mx[ci[1]:(ci[1]+nv-1)]
      transition_right <- -0.5*mx[(ci[1]+nv):ci[2]] %*% solve(Sigma0[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2], drop = FALSE]) %*% mx[(ci[1]+nv):ci[2]]
      # contribution of g_{t,u}
      obs_covariance <- diag(sqrt(exp(mx[ci[1]:ci[2]]))) %*% SigmaV[ci[1]:ci[2], ci[1]:ci[2]] %*% diag(sqrt(exp(mx[ci[1]:ci[2]])))
      obs_node <- -0.5*obs[ci[1]:ci[2]] %*% solve(obs_covariance) %*% obs[ci[1]:ci[2]]
      obs_left <- -0.5*obs[ci[1]:(ci[1]+nv-1)] %*% solve(obs_covariance[1:nv, 1:nv, drop = FALSE]) %*% obs[ci[1]:(ci[1]+nv-1)]
      obs_right <- -0.5*obs[(ci[1]+nv):ci[2]] %*% solve(obs_covariance[(nv+1):(2*nv), (nv+1):(2*nv), drop = FALSE]) %*% obs[(ci[1]+nv):ci[2]]
      lWmix_perm[n] <- transition_node - transition_left - transition_right + obs_node - obs_left - obs_right
    }
    if(u == 1){
      lWmix_perm <- lWmix_perm + lW[, (nchild*(i-1)+1)] + lW[new_perm, i*nchild]
    }
    permutation <- c(permutation, new_perm)
    max.lWmix <- max(lWmix_perm)
    Wmix <- exp(lWmix_perm - max.lWmix)
    # build ESS
    ess_s <- ess_s + sum(Wmix)
    ess_ss <- ess_ss + sum(Wmix^2)
    ess <- ess_s^2/ess_ss
    lWmix <- c(lWmix, lWmix_perm)
  }
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(rep(1:Nparticles, times = theta)[indices], permutation[indices]),
              "target_reached" = (ess >= ess_target), "resampled_particles_lW" = lWmix[indices]))
}


marginal_msv_light <- function(ess_target, i, u, nv, ci, lW, Nparticles, theta, SigmaX, SigmaV, SigmaUV, phi,
                                           x, obs_current, obs_past, history){
  nchild <- 2
  # mixture weights
  lWmix <- rep(0, times = Nparticles)
  # precompute mean vector
  mu <- t(apply(history, 1, function(s){ phi*s + SigmaUV %*% solve(SigmaV) %*% diag(1/sqrt(exp(s))) %*% obs_past}))
  for (n in 1:Nparticles) {
    # contribution of f_{t, u}
    mx <- x[n, ]
    centred_mx <- -sweep(mu, 2, mx, '-')
    # # contribution of f_{t, u}
    transition_node <- -0.5*mean(apply(centred_mx[, ci[1]:ci[2]], 1, function(s) {s %*% solve(SigmaX[ci[1]:ci[2], ci[1]:ci[2]]) %*% s }))
    transition_left <- -0.5*mean(apply(centred_mx[, ci[1]:(ci[1]+nv-1), drop = F], 1, function(s) {s %*% solve(SigmaX[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1)]) %*% s }))
    transition_right <- -0.5*mean(apply(centred_mx[, (ci[1]+nv):ci[2], drop = F], 1, function(s) {s %*% solve(SigmaX[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2]]) %*% s }))
    # contribution of g_{t,u}
    obs_covariance <- diag(sqrt(exp(mx[ci[1]:ci[2]]))) %*% SigmaV[ci[1]:ci[2], ci[1]:ci[2]] %*% diag(sqrt(exp(mx[ci[1]:ci[2]])))
    obs_node <- -0.5*obs_current[ci[1]:ci[2]] %*% solve(obs_covariance) %*% obs_current[ci[1]:ci[2]]
    obs_left <- -0.5*obs_current[ci[1]:(ci[1]+nv-1)] %*% solve(obs_covariance[1:nv, 1:nv, drop = FALSE]) %*% obs_current[ci[1]:(ci[1]+nv-1)]
    obs_right <- -0.5*obs_current[(ci[1]+nv):ci[2]] %*% solve(obs_covariance[(nv+1):(2*nv), (nv+1):(2*nv), drop = FALSE]) %*% obs_current[(ci[1]+nv):ci[2]]
    lWmix[n] <- transition_node - transition_left - transition_right + obs_node - obs_left - obs_right
  }
  if(u == 1){
    lWmix <- lWmix + lW[, (nchild*(i-1)+1)] + lW[, i*nchild]
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # build ESS
  ess_s <- sum(Wmix)
  ess_ss <- sum(Wmix^2)
  ess <- ess_s^2/ess_ss
  # first permutation
  permutation <- 1:Nparticles
  theta <- 1
  while (ess < ess_target & theta <= ceiling(sqrt(Nparticles))) {
    theta <- theta+1
    new_perm <- sample.int(Nparticles)
    # mixture weights
    lWmix_perm <- rep(0, times = Nparticles)
    for (n in 1:Nparticles) {
      # contribution of f_{t, u}
      mx <- x[n, ]
      mx[(ci[1]+nv):ci[2]] <- x[new_perm[n], (ci[1]+nv):ci[2]]
      centred_mx <- -sweep(mu, 2, mx, '-')
      # # contribution of f_{t, u}
      transition_node <- -0.5*mean(apply(centred_mx[, ci[1]:ci[2]], 1, function(s) {s %*% solve(SigmaX[ci[1]:ci[2], ci[1]:ci[2]]) %*% s }))
      transition_left <- -0.5*mean(apply(centred_mx[, ci[1]:(ci[1]+nv-1), drop = F], 1, function(s) {s %*% solve(SigmaX[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1)]) %*% s }))
      transition_right <- -0.5*mean(apply(centred_mx[, (ci[1]+nv):ci[2], drop = F], 1, function(s) {s %*% solve(SigmaX[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2]]) %*% s }))
      # contribution of g_{t,u}
      obs_covariance <- diag(sqrt(exp(mx[ci[1]:ci[2]]))) %*% SigmaV[ci[1]:ci[2], ci[1]:ci[2]] %*% diag(sqrt(exp(mx[ci[1]:ci[2]])))
      obs_node <- -0.5*obs_current[ci[1]:ci[2]] %*% solve(obs_covariance) %*% obs_current[ci[1]:ci[2]]
      obs_left <- -0.5*obs_current[ci[1]:(ci[1]+nv-1)] %*% solve(obs_covariance[1:nv, 1:nv, drop = FALSE]) %*% obs_current[ci[1]:(ci[1]+nv-1)]
      obs_right <- -0.5*obs_current[(ci[1]+nv):ci[2]] %*% solve(obs_covariance[(nv+1):(2*nv), (nv+1):(2*nv), drop = FALSE]) %*% obs_current[(ci[1]+nv):ci[2]]
      lWmix[n] <- transition_node - transition_left - transition_right + obs_node - obs_left - obs_right
    }
    if(u == 1){
      lWmix_perm <- lWmix_perm + lW[, (nchild*(i-1)+1)] + lW[new_perm, i*nchild]
    }
    permutation <- c(permutation, new_perm)
    max.lWmix <- max(lWmix_perm)
    Wmix <- exp(lWmix_perm - max.lWmix)
    # build ESS
    ess_s <- ess_s + sum(Wmix)
    ess_ss <- ess_ss + sum(Wmix^2)
    ess <- ess_s^2/ess_ss
    lWmix <- c(lWmix, lWmix_perm)
  }
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(rep(1:Nparticles, times = theta)[indices], permutation[indices]),
              "target_reached" = (ess >= ess_target), "resampled_particles_lW" = lWmix[indices]))
}
