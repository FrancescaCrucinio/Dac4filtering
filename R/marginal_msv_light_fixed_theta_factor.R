# Lightweight resampling for linear Gaussian SSM
marginal_msv_light_fixed_theta_factor <- function(i, u, nv, ci, W, Nparticles, theta, mu, Phi, Sigma, Lambda, x, history, obs){
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
  all_nodes_means <- sweep(sweep(history[, ci[1]:ci[2]], 2, mu[ci[1]:ci[2]]) %*% diag(Phi[ci[1]:ci[2]]), 2,
                           mu[ci[1]:ci[2]], FUN = "+")
  all_nodes_means_left <- sweep(sweep(history[, ci[1]:(ci[1]+nv-1), drop = FALSE], 2, mu[ci[1]:(ci[1]+nv-1), drop = FALSE]) %*% diag(Phi[ci[1]:(ci[1]+nv-1)], nrow = nv),
                                2, mu[ci[1]:(ci[1]+nv-1), drop = FALSE], FUN = "+")
  all_nodes_means_right <- sweep(sweep(history[, (ci[1]+nv):ci[2], drop = FALSE], 2, mu[(ci[1]+nv):ci[2], drop = FALSE]) %*% diag(Phi[(ci[1]+nv):ci[2]], nrow = length((ci[1]+nv):ci[2])),
                                 2, mu[(ci[1]+nv):ci[2], drop = FALSE], FUN = "+")

  # mixture weights
  lWmix <- rep(0, times = theta*Nparticles)
  for (n in 1:(theta*Nparticles)) {
    # contribution of f_{t, u}
    mx <- x[indices1[n], ]
    mx[(ci[1]+nv):ci[2]] <- x[indices2[n], (ci[1]+nv):ci[2]]
    transition_node <- mean(exp(-rowSums(sweep(all_nodes_means, 2, mx[ci[1]:ci[2]])^2 %*% diag(1/Sigma[ci[1]:ci[2]]^2))))
    transition_left <- mean(exp(-rowSums(sweep(all_nodes_means_left, 2, mx[ci[1]:(ci[1]+nv-1)])^2 %*% diag(1/Sigma[ci[1]:(ci[1]+nv-1)]^2, nrow = nv))))
    transition_right <- mean(exp(-rowSums(sweep(all_nodes_means_right, 2, mx[(ci[1]+nv):ci[2]])^2 %*% diag(1/Sigma[(ci[1]+nv):ci[2]]^2, nrow = length((ci[1]+nv):ci[2])))))
    # contribution of g_{t, u}
    factor_covariance <- Lambda %*% diag(exp(mx[(p+1):d])) %*% t(Lambda)
    obs_covariance <- factor_covariance[ci[1]:p, ci[1]:p] + diag(exp(mx[ci[1]:p]))
    obs_covariance <- (obs_covariance + t(obs_covariance))/2
    if(p-ci[1]-nv > 0){ # if there are observations on the right
      obs_covariance_right <- factor_covariance[(ci[1]+nv):p, (ci[1]+nv):p] + diag(exp(mx[(ci[1]+nv):p]))
      obs_right <- -0.5*obs[(ci[1]+nv):p] %*% solve(obs_covariance_right) %*% obs[(ci[1]+nv):p]

    } else {
      obs_right <- 0
    }
    obs_node <- -0.5*obs[ci[1]:p] %*% solve(obs_covariance) %*% obs[ci[1]:p]
    obs_left <- -0.5*obs[ci[1]:(ci[1]+nv-1)] %*% diag(1/exp(mx[ci[1]:(ci[1]+nv-1)])) %*% obs[ci[1]:(ci[1]+nv-1)]
    lWmix[n] <- log(transition_node) - log(transition_right) - log(transition_left) +
      obs_node - obs_left - obs_right
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}
