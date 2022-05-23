# Lightweight resampling DaC for linear Gaussian SSM
marginal_dac_lgssm_lightweight <- function(xOld, obs, tau, lambda, sigmaY, adaptive = FALSE){
  # number of samples for lightweight mixture
  theta <- ceiling(sqrt(Nparticles))
  # dimension and number of particles
  d <- ncol(xOld)
  Nparticles <- nrow(xOld)
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # leaves
  # number of variables
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  lW <- matrix(0, nrow = Nparticles, ncol = d)
  W <- matrix(0, nrow = Nparticles, ncol = d)
  for (i in 1:nchild^nlevels){
    sample_from_past <- sample.int(Nparticles, Nparticles, replace = TRUE)
    if (i == 1) {
      x[, i] <- 0.5*xOld[sample_from_past, i] + rnorm(Nparticles)/sqrt(tau)
    } else {
      x[, i] <- 0.5*tau*xOld[sample_from_past, i]/(tau+lambda) + rnorm(Nparticles)/sqrt(tau+lambda)
    }
    lW[, i] <- -0.5*(obs[i] - x[, i])^2/sigmaY - 0.5*log(2*pi*sigmaY)
    max.lW <- max(lW[, i])
    W[, i] <- exp(lW[, i] - max.lW)
    W[, i] <- W[, i]/sum(W[, i])
  }

  # loop over tree levels excluding leaves
  for (u in 1:nlevels){
    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    # number of variables in each node
    nvNew <- nchild^u
    # updated particles/normalizing constant
    xNew <- matrix(0, nrow = Nparticles, ncol = d)

    for (i in 1:nodes){
      # get children indices
      ci <- child_indices(i, nvNew)
      # lightweight mixture resampling
      if(adaptive){
        out <- marginal_lgssm_light_adaptive(Nparticles, i, u, nv, ci, lW, Nparticles, lambda, tau, x, xOld)
      } else {
        out <- marginal_lgssm_light(i, u, nv, ci, W, Nparticles, theta, lambda, tau, x, xOld)
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

# Mixture resampling DaC for linear Gaussian SSM
dac_lgssm <- function(xOld, obs, tau, lambda, sigmaY){
  # dimension and number of particles
  d <- ncol(xOld)
  Nparticles <- nrow(xOld)
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # number of variables
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  lW <- matrix(0, nrow = Nparticles, ncol = d)
  for (i in 1:nchild^nlevels){
    if (i == 1) {
      x[, i] <- 0.5*xOld[, i] + rnorm(Nparticles)/sqrt(tau)
    } else {
      x[, i] <- 0.5*tau*xOld[, i]/(tau+lambda) + rnorm(Nparticles)/sqrt(tau+lambda)
    }
    lW[, i] <- -0.5*(obs[i] - x[, i])^2/sigmaY - 0.5*log(2*pi*sigmaY)
    max.lW <- max(lW[, i])
    W <- exp(lW[, i] - max.lW)
  }

  # loop over tree levels excluding leaves
  for (u in 1:nlevels){
    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    # number of variables in each node
    nvNew <- nchild^u

    # updated particles/normalizing constant
    xNew <- matrix(0, nrow = Nparticles, ncol = d)

    for (i in 1:nodes){

      # mixture weights
      lWmix <- matrix(0, ncol = Nparticles, nrow = Nparticles)
      ci <- child_indices(i, nvNew)

      if(u == 1){
        for (n1 in 1:Nparticles) {
          for (n2 in 1:Nparticles) {
            lWmix[n1, n2] <- lW[n1, (nchild*(i-1)+1)] + lW[n2, i*nchild] -
              -0.5*lambda * (lambda *x[n1, (ci[1]+nv-1)]^2/(tau+lambda) -
                               2*x[n1, (ci[1]+nv-1)] * x[n2, (ci[1]+nv)]) +
              log(mean(exp(-0.5*lambda*tau*xOld[n2, ci[1]+nv] * x[n1, (ci[1]+nv-1)]/(tau+lambda))))
          }
        }
        max.lWmix <- max(lWmix)
        Wmix <- exp(lWmix - max.lWmix)
      } else {
        for (n1 in 1:Nparticles) {
          for (n2 in 1:Nparticles) {
            lWmix[n1, n2] <- -0.5*lambda * (lambda *x[n1, (ci[1]+nv-1)]^2/(tau+lambda) -
                                              2*x[n1, (ci[1]+nv-1)] * x[n1, (ci[1]+nv)]) +
              log(mean(exp(-0.5*lambda*tau*xOld[n2, ci[1]+nv] * x[n1, (ci[1]+nv-1)]/(tau+lambda))))
          }
        }
        max.lWmix <- max(lWmix)
        Wmix <- exp(lWmix - max.lWmix)
      }
      Wmix <- Wmix/sum(Wmix)
      # resampling
      indices <- stratified_resample(Wmix, Nparticles)
      # get row/column indices
      res <- rc_indices(indices, Nparticles)
      for(n in 1:Nparticles) {
        # update particles
        xNew[n, ci[1]:ci[2]] <- c(x[res[n, 1], ci[1]:(ci[1]+nv-1)], x[res[n, 2], (ci[1]+nv):ci[2]])
      }
    }
    x <- xNew
    nv <- nvNew
  }
  return(x)
}

# Linear cost DaC for linear Gaussian SSM
dac_lgssm_lc <- function(xOld, obs, tau, lambda, sigmaY){
  # dimension and number of particles
  d <- ncol(xOld)
  Nparticles <- nrow(xOld)
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # number of variables
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  W <- matrix(0, nrow = Nparticles, ncol = d)
  for (i in 1:nchild^nlevels){
    if (i == 1) {
      x[, i] <- 0.5*xOld[, i] + rnorm(Nparticles)/sqrt(tau)
    } else {
      x[, i] <- 0.5*tau*xOld[, i]/(tau+lambda) + rnorm(Nparticles)/sqrt(tau+lambda)
    }
    lW <- -0.5*(obs[i] - x[, i])^2/sigmaY - 0.5*log(2*pi*sigmaY)
    max.lW <- max(lW)
    W[, i] <- exp(lW - max.lW)
    W[, i] <- W[, i]/sum(W[, i])
  }

  # loop over tree levels excluding leaves
  for (u in 1:nlevels){
    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    # number of variables in each node
    nvNew <- nchild^u

    # updated particles/normalizing constant
    xNew <- matrix(0, nrow = Nparticles, ncol = d)
    WNew <- matrix(0, nrow = Nparticles, ncol = nodes)

    for (i in 1:nodes){
      # get children indices
      ci <- child_indices(i, nvNew)
      # resample on each children
      # child 1
      indices1 <- stratified_resample(W[, nchild*(i-1)+1], Nparticles)
      # child 2 (with random permutation)
      indices2 <- sample(stratified_resample(W[, nchild*i], Nparticles))
      # weights
      lW <- -0.5*lambda * (lambda *x[indices1, (ci[1]+nv-1)]^2/(tau+lambda) -
                             2*x[indices1, (ci[1]+nv-1)] * x[indices2, (ci[1]+nv)]) +
        log(mean(exp(-0.5*lambda*tau*xOld[indices2, ci[1]+nv] * x[indices1, (ci[1]+nv-1)]/(tau+lambda))))
      max.lW <- max(lW)
      WNew[, i] <- exp(lW - max.lW)
      # update particles
      xNew[, ci[1]:ci[2]] <- cbind(x[indices1, ci[1]:(ci[1]+nv-1)], x[indices2, (ci[1]+nv):ci[2]])
    }
    x <- xNew
    nv <- nvNew
    W <- WNew
  }
  return(x)
}

# Lightweight resampling for linear Gaussian SSM
marginal_lgssm_light <- function(i, u, nv, ci, W, Nparticles, theta, lambda, tau, x, xOld){
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
  lWmix <- -0.5*lambda * (lambda *x[indices1, (ci[1]+nv-1)]^2/(tau+lambda) -
                            2*x[indices1, (ci[1]+nv-1)] * x[indices2, (ci[1]+nv)]) +
            log(mean(exp(-0.5*lambda*tau*xOld[indices2, ci[1]+nv] * x[indices1, (ci[1]+nv-1)]/(tau+lambda))))
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}


# Adaptive Lightweight resampling for linear Gaussian SSM
marginal_lgssm_light_adaptive <- function(ess_target, i, u, nv, ci, lW, Nparticles, lambda, tau, x, xOld){
  # binary tree
  nchild <- 2
  # mixture weights
  if(u == 1){
    lWmix <- lW[, (nchild*(i-1)+1)] + lW[, i*nchild] - 0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) -
                                                                      2*x[, (ci[1]+nv-1)] * x[, (ci[1]+nv)]) +
      log(mean(exp(-0.5*lambda*tau*xOld[, ci[1]+nv] * x[, (ci[1]+nv-1)]/(tau+lambda))))
  } else{
    lWmix <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) -
                              2*x[, (ci[1]+nv-1)] * x[, (ci[1]+nv)]) +
      log(mean(exp(-0.5*lambda*tau*xOld[, ci[1]+nv] * x[, (ci[1]+nv-1)]/(tau+lambda))))
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
    if(u == 1){
      lWmix_perm <- lW[, (nchild*(i-1)+1)] + lW[new_perm, i*nchild] -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) -
                                                                                     2*x[, (ci[1]+nv-1)] * x[new_perm, (ci[1]+nv)]) +
        log(mean(exp(-0.5*lambda*tau*xOld[new_perm, ci[1]+nv] * x[, (ci[1]+nv-1)]/(tau+lambda))))
    } else{
      lWmix_perm <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) -
                                     2*x[, (ci[1]+nv-1)] * x[new_perm, (ci[1]+nv)]) +
        log(mean(exp(-0.5*lambda*tau*xOld[new_perm, ci[1]+nv] * x[, (ci[1]+nv-1)]/(tau+lambda))))
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
  write.table(data.frame("u" = u, "theta" = theta), file = "data/adaptive_lgssm.csv", sep = ",", append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(rep(1:Nparticles, times = theta)[indices], permutation[indices]),
              "target_reached" = (ess >= ess_target), "resampled_particles_lW" = lWmix[indices]))
}
