# Lightweight resampling for linear Gaussian SSM
marginal_lgssm_light <- function(i, u, nv, ci, W, Nparticles, theta, lambda, tau, x, history, memory = FALSE){
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
  if(memory){
    lWmix <- rep(0, times = theta*Nparticles)
    for (n in 1:(theta*Nparticles)) {
      lWmix[n] <- -0.5*lambda * (lambda *x[indices1[n], (ci[1]+nv-1)]^2/(tau+lambda) -
                                   2*x[indices1[n], (ci[1]+nv-1)] * x[indices2[n], (ci[1]+nv)])+
        log(mean(exp(-0.5*(x[indices1[n], (ci[1]+nv)]-history[, ci[1]+nv])^2/(tau+lambda) -
                       0.5*lambda*tau*history[, ci[1]+nv] * x[indices1[n], (ci[1]+nv-1)]/(tau+lambda)))) -
        log(mean(exp(-0.5*(x[indices1[n], (ci[1]+nv)]-history[, ci[1]+nv])^2/(tau+lambda))))
    }
  } else {
    lWmix <- -0.5*lambda * (lambda *x[indices1, (ci[1]+nv-1)]^2/(tau+lambda) -
                              2*x[indices1, (ci[1]+nv-1)] * x[indices2, (ci[1]+nv)]) +
      log(rowMeans(exp(-0.5*outer(x[indices1, (ci[1]+nv)], history[, ci[1]+nv], "-")^2/(tau+lambda) -
                         0.5*lambda*tau*outer(x[indices1, (ci[1]+nv-1)], history[, ci[1]+nv], "*")/(tau+lambda)))) -
      log(rowMeans(exp(-0.5*outer(x[indices1, (ci[1]+nv)], history[, ci[1]+nv], "-")^2/(tau+lambda))))
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}


# Adaptive Lightweight resampling for linear Gaussian SSM
marginal_lgssm_light_adaptive <- function(ess_target, i, u, nv, ci, lW, Nparticles, lambda, tau, x, history){
  # binary tree
  nchild <- 2
  integral_per_dimension <- array(0, dim = c(Nparticles, Nparticles, 2*nv))
  # marginalize out past (each row corresponds to a particle)
  for (j in 2:(2*nv)) {
    dimension <- (ci[1]:ci[2])[j]
    integral_per_dimension[, , j] <- -0.5*(tau+lambda)*outer(x[, dimension], 0.5*tau*history[, dimension]/(tau+lambda), "-")^2 -
      0.5*lambda*tau*outer(x[, dimension-1], history[, dimension], "*")/(tau+lambda)
  }
  if(ci[1] == 1){
    integral_per_dimension[, , 1] <- -0.5*tau*outer(x[, 1], 0.5*history[, 1], "-")^2
  } else {
    integral_per_dimension[, , 1] <- -0.5*(tau+lambda)*outer(x[, ci[1]], 0.5*tau*history[, ci[1]]/(tau+lambda), "-")^2
  }
  integral_merged <- rowMeans(exp(rowSums(integral_per_dimension, dims = 2)))
  integral_left <- rowMeans(exp(rowSums(integral_per_dimension[, , 1:nv, drop = FALSE], dims = 2)))
  integral_right <- rowMeans(exp(rowSums(integral_per_dimension[, , (nv+1):(2*nv), drop = FALSE], dims = 2)))
  integral_merged <- ifelse(all(integral_merged < .Machine$double.eps), rep(0, Nparticles), log(integral_merged))
  integral_left <- ifelse(all(integral_left< .Machine$double.eps), rep(0, Nparticles), log(integral_left))
  integral_right <- ifelse(all(integral_right< .Machine$double.eps), rep(0, Nparticles), log(integral_right))
  if(u == 1){
    lWmix <- lW[, (nchild*(i-1)+1)] + lW[, i*nchild] + lambda*x[, (ci[1]+nv-1)] * x[, (ci[1]+nv)]
  } else{
    lWmix <- lambda*x[, (ci[1]+nv-1)] * x[, (ci[1]+nv)]
  }
  lWmix <- lWmix + integral_merged - integral_left - integral_right
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
    integral_per_dimension[, , (nv+1):(2*nv)] <- 0
    # marginalize out past (each row corresponds to a particle)
    if(u == 1){
      integral_per_dimension[, , 2*nv] <- -0.5*(tau+lambda)*outer(x[new_perm, 2*nv], 0.5*tau*history[, 2*nv]/(tau+lambda), "-")^2 -
        0.5*lambda*tau*outer(x[, 2*nv-1], history[, 2*nv], "*")/(tau+lambda)
    } else {
      for (j in (nv+2):(2*nv)) {
        dimension <- (ci[1]:ci[2])[j]
        integral_per_dimension[, , j] <- -0.5*(tau+lambda)*outer(x[new_perm, dimension], 0.5*tau*history[, dimension]/(tau+lambda), "-")^2 -
          0.5*lambda*tau*outer(x[new_perm, dimension-1], history[, dimension], "*")/(tau+lambda)
      }
    }
    integral_merged <- rowMeans(exp(rowSums(integral_per_dimension, dims = 2)))
    integral_right <- rowMeans(exp(rowSums(integral_per_dimension[, , (nv+1):(2*nv), drop = FALSE], dims = 2)))
    integral_merged <- ifelse(all(integral_merged < .Machine$double.eps), rep(0, Nparticles), log(integral_merged))
    integral_right <- ifelse(all(integral_right< .Machine$double.eps), rep(0, Nparticles), log(integral_right))
    if(u == 1){
      lWmix_perm <- lW[, (nchild*(i-1)+1)] + lW[new_perm, i*nchild] + lambda*x[, (ci[1]+nv-1)] * x[new_perm, (ci[1]+nv)]
    } else{
      lWmix_perm <- lambda*x[, (ci[1]+nv-1)] * x[new_perm, (ci[1]+nv)]
    }
    lWmix_perm <- lWmix_perm + integral_merged - integral_left - integral_right
    permutation <- c(permutation, new_perm)
    max.lWmix <- max(lWmix_perm)
    Wmix <- exp(lWmix_perm - max.lWmix)
    # build ESS
    ess_s <- ess_s + sum(Wmix)
    ess_ss <- ess_ss + sum(Wmix^2)
    ess <- ess_s^2/ess_ss
    lWmix <- c(lWmix, lWmix_perm)
  }
  # write.table(data.frame("u" = u, "theta" = theta), file = "data/adaptive_lgssm.csv", sep = ",", append = TRUE, quote = FALSE,
  #             col.names = FALSE, row.names = FALSE)
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(rep(1:Nparticles, times = theta)[indices], permutation[indices]),
              "target_reached" = (ess >= ess_target), "resampled_particles_lW" = lWmix[indices]))
}
