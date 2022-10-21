# Lightweight resampling for linear Gaussian SSM
marginal_msv_light_fixed_theta <- function(i, u, nv, ci, W, Nparticles, theta, lambda, tau, x, history){
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
  all_nodes_means <- Phi[ci[1]:ci[2]] * sweep(history[ci[1]:ci[2], ], mu[ci[1]:ci[2]], 2) + mu[ci[1]:ci[2]]

  # mixture weights
  for (n in 1:(theta*Nparticles)) {
    mx <- x[indices1[n], ]
    mx[(ci[1]+nv):ci[2]] <- x[indices2[n], (ci[1]+nv):ci[2]]
  }

  lWmix <-
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}
