# Lightweight resampling for linear Gaussian SSM
marginal_lgssm_light <- function(i, u, nv, ci, W, Nparticles, theta, lambda, tau, x, history){
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
  integral_merged <- rep(0, times = Nparticles)
  integral_left <- rep(0, times = Nparticles)
  integral_right <- rep(0, times = Nparticles)
  for (n in 1:(theta*Nparticles)) {
    mx <- x[indices1[n], ]
    mx[(ci[1]+nv):ci[2]] <- x[indices2[n], (ci[1]+nv):ci[2]]
    integral_per_dimension <- matrix(0, nrow = Nparticles, ncol = 2*nv)
    for (j in 2:(2*nv)) {
      dimension <- (ci[1]:ci[2])[j]
      integral_per_dimension[, j] <- -0.5*(tau+lambda)*(mx[dimension] - 0.5*tau*history[, dimension]/(tau+lambda))^2 -
        0.5*lambda*tau*(mx[dimension-1] * history[, dimension])/(tau+lambda)
    }
    if(ci[1] == 1){
      integral_per_dimension[, 1] <- -0.5*tau*(mx[1] - 0.5*history[, 1])^2
    } else {
      integral_per_dimension[, 1] <- -0.5*(tau+lambda)*(mx[ci[1]] - 0.5*tau*history[, ci[1]]/(tau+lambda))^2
    }
    integral_merged[n] <- mean(exp(rowSums(integral_per_dimension)))
    integral_left[n] <- mean(exp(rowSums(integral_per_dimension[, 1:nv, drop = FALSE])))
    integral_right[n] <- mean(exp(rowSums(integral_per_dimension[, (nv+1):(2*nv), drop = FALSE])))
  }
  integral_merged <- ifelse(all(integral_merged < .Machine$double.eps), rep(0, Nparticles), log(integral_merged))
  integral_right <- ifelse(all(integral_right < .Machine$double.eps), rep(0, Nparticles), log(integral_right))
  integral_left <- ifelse(all(integral_left < .Machine$double.eps), rep(0, Nparticles), log(integral_left))
  lWmix <- - 0.5*lambda * (lambda *x[indices1, (ci[1]+nv-1)]^2/(tau+lambda) -
                             2*x[indices1, (ci[1]+nv-1)] * x[indices2, (ci[1]+nv)])
  lWmix <- lWmix + integral_merged - integral_left - integral_right
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}


# Lightweight resampling for linear Gaussian SSM
marginal_lgssm_light_vectorized <- function(i, u, nv, ci, W, Nparticles, theta, lambda, tau, x, history){
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
  integral_per_dimension <- array(0, dim = c(theta*Nparticles, Nparticles, 2*nv))
  if(u == 1){
    j <- 2*nv
    dimension <- (ci[1]:ci[2])[j]
    integral_per_dimension[, , j] <- -0.5*(tau+lambda)*outer(x[indices2, dimension], 0.5*tau*history[, dimension]/(tau+lambda), "-")^2 -
      0.5*lambda*tau*outer(x[indices1, dimension-1], history[, dimension], "*")/(tau+lambda)
  } else {
    for (j in 2:nv) {
      dimension <- (ci[1]:ci[2])[j]
      integral_per_dimension[, , j] <- -0.5*(tau+lambda)*outer(x[indices1, dimension], 0.5*tau*history[, dimension]/(tau+lambda), "-")^2 -
        0.5*lambda*tau*outer(x[indices1, dimension-1], history[, dimension], "*")/(tau+lambda)
    }
    for (j in (nv+2):(2*nv)) {
      dimension <- (ci[1]:ci[2])[j]
      integral_per_dimension[, , j] <- -0.5*(tau+lambda)*outer(x[indices2, dimension], 0.5*tau*history[, dimension]/(tau+lambda), "-")^2 -
        0.5*lambda*tau*outer(x[indices2, dimension-1], history[, dimension], "*")/(tau+lambda)
    }
    j <- nv+1
    dimension <- (ci[1]:ci[2])[j]
    integral_per_dimension[, , j] <- -0.5*(tau+lambda)*outer(x[indices2, dimension], 0.5*tau*history[, dimension]/(tau+lambda), "-")^2 -
      0.5*lambda*tau*outer(x[indices1, dimension-1], history[, dimension], "*")/(tau+lambda)
  }
  if(ci[1] == 1){
    integral_per_dimension[, , 1] <- -0.5*tau*outer(x[indices1, 1], 0.5*history[, 1], "-")^2
  } else {
    integral_per_dimension[, , 1] <- -0.5*(tau+lambda)*outer(x[indices1, ci[1]], 0.5*tau*history[, ci[1]]/(tau+lambda), "-")^2
  }
  integral_merged <- rowMeans(exp(rowSums(integral_per_dimension, dims = 2)))
  integral_left <- rowMeans(exp(rowSums(integral_per_dimension[, , 1:nv, drop = FALSE], dims = 2)))
  integral_right <- rowMeans(exp(rowSums(integral_per_dimension[, , (nv+1):(2*nv), drop = FALSE], dims = 2)))
  integral_merged <- ifelse(all(integral_merged < .Machine$double.eps), rep(0, Nparticles), log(integral_merged))
  integral_right <- ifelse(all(integral_right < .Machine$double.eps), rep(0, Nparticles), log(integral_right))
  integral_left <- ifelse(all(integral_left < .Machine$double.eps), rep(0, Nparticles), log(integral_left))
  lWmix <- -0.5*lambda * (lambda *x[indices1, (ci[1]+nv-1)]^2/(tau+lambda) - 2*x[indices1, (ci[1]+nv-1)] * x[indices2, (ci[1]+nv)]) +
          integral_merged - integral_left - integral_right
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}
