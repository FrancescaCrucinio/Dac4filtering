dac_lgssm <- function(xOld, obs, tau, lambda, sigmaY, Sigma.det){
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
  lZ <- rep(0, times = d)
  for (i in 1:nchild^nlevels){
    if (i == 1 | i == d) {
      x[, i] <- 0.5*xOld[, i] + rnorm(Nparticles)/sqrt(tau+lambda)
    } else {
      x[, i] <- 0.5*xOld[, i] + rnorm(Nparticles)/sqrt(tau+2*lambda)
    }
    lW[, i] <- -0.5*(obs[i] - x[, i])^2/sigmaY - 0.5*log(2*pi*sigmaY)
    max.lW <- max(lW[, i])
    W <- exp(lW[, i] - max.lW)
    lZ[i] <- log(mean(W)) + max.lW
  }

  # loop over tree levels excluding leaves
  for (u in 1:nlevels){
    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    # number of variables in each node
    nvNew <- nchild^u

    # updated particles/normalizing constant
    xNew <- matrix(0, nrow = Nparticles, ncol = d)
    xOldNew <- matrix(0, nrow = Nparticles, ncol = d)
    lZNew <- rep(0, times = nodes)

    for (i in 1:nodes){

      # mixture weights
      lWmix <- matrix(0, ncol = Nparticles, nrow = Nparticles)
      ci <- child_indices(i, nvNew)

      if(u == 1){
        for (n1 in 1:Nparticles) {
          for (n2 in 1:Nparticles) {
            lWmix[n1, n2] <- lW[n1, (nchild*(i-1)+1)] + lW[n2, i*nchild] +
              lambda * (x[n1, (ci[1]+nv-1)] - 0.5*xOld[n1, (ci[1]+nv-1)]) * (x[n2, (ci[1]+nv)] - 0.5*xOld[n2, (ci[1]+nv)])
          }
        }
        max.lWmix <- max(lWmix)
        Wmix <- exp(lWmix - max.lWmix)
        lZNew[i] <- log(mean(Wmix)) + max.lWmix -
          0.5*Sigma.det[[u]][nchild*(i-1)+1] - 0.5*Sigma.det[[u]][i*nchild] + 0.5*Sigma.det[[u+1]][i]
      } else {
        for (n1 in 1:Nparticles) {
          for (n2 in 1:Nparticles) {
            lWmix[n1, n2] <- lambda * (x[n1, (ci[1]+nv-1)] - 0.5*xOld[n1, (ci[1]+nv-1)]) * (x[n2, (ci[1]+nv)] - 0.5*xOld[n2, (ci[1]+nv)])
          }
        }
        max.lWmix <- max(lWmix)
        Wmix <- exp(lWmix - max.lWmix)
        lZNew[i] <- log(mean(Wmix)) + max.lWmix + lZ[(nchild*(i-1)+1)] + lZ[i*nchild] -
          0.5*Sigma.det[[u]][nchild*(i-1)+1] - 0.5*Sigma.det[[u]][i*nchild] + 0.5*Sigma.det[[u+1]][i]
      }
      Wmix <- Wmix/sum(Wmix)
      # resampling
      indices <- mult_resample(Wmix, Nparticles)
      # get row/column indices
      res <- rc_indices(indices, Nparticles)
      for(n in 1:Nparticles) {
        # update particles
        xNew[n, ci[1]:ci[2]] <- c(x[res[n, 1], ci[1]:(ci[1]+nv-1)], x[res[n, 2], (ci[1]+nv):ci[2]])
        # update xOld
        xOldNew[n, ci[1]:ci[2]] <- c(xOld[res[n, 1], ci[1]:(ci[1]+nv-1)], xOld[res[n, 2], (ci[1]+nv):ci[2]])
      }
    }
    x <- xNew
    xOld <- xOldNew
    lZ <- lZNew
    nv <- nvNew
  }
  return(cbind(x, lZ))
}

dac_lgssm_lightweight <- function(xOld, obs, tau, lambda, sigmaY, Sigma.det, m){
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
  W <- matrix(0, nrow = Nparticles, ncol = d)
  lZ <- rep(0, times = d)
  for (i in 1:nchild^nlevels){
    if (i == 1 | i == d) {
      x[, i] <- 0.5*xOld[, i] + rnorm(Nparticles)/sqrt(tau+lambda)
    } else {
      x[, i] <- 0.5*xOld[, i] + rnorm(Nparticles)/sqrt(tau+2*lambda)
    }
    lW <- -0.5*(obs[i] - x[, i])^2/sigmaY - 0.5*log(2*pi*sigmaY)
    max.lW <- max(lW)
    W[, i] <- exp(lW - max.lW)
    lZ[i] <- log(mean(W[, i])) + max.lW
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
    xOldNew <- matrix(0, nrow = Nparticles, ncol = d)
    lZNew <- rep(0, times = nodes)

    for (i in 1:nodes){
      # get children indices
      ci <- child_indices(i, nvNew)
      # resample on each children
      # child 1
      indices1 <- mult_resample(W[, nchild*(i-1)+1], m*Nparticles)
      # child 2 (with random permutation)
      indices2 <- sample(mult_resample(W[, nchild*i], m*Nparticles))
      # mixture weights
      lWmix <- lambda * (x[indices1, (ci[1]+nv-1)] - 0.5*xOld[indices1, (ci[1]+nv-1)]) * 
        (x[indices2, (ci[1]+nv)] - 0.5*xOld[indices2, (ci[1]+nv)])
      max.lWmix <- max(lWmix)
      Wmix <- exp(lWmix - max.lWmix)
      lZNew[i] <- lZ[(nchild*(i-1)+1)] + lZ[i*nchild] + log(mean(Wmix)) + max.lWmix -
        0.5*Sigma.det[[u]][nchild*(i-1)+1] - 0.5*Sigma.det[[u]][i*nchild] + 0.5*Sigma.det[[u+1]][i]
      # resampling the new population
      indices <- mult_resample(Wmix/sum(Wmix), Nparticles)
      # update particles
      xNew[, ci[1]:ci[2]] <- rbind(x[indices1[indices], ci[1]:(ci[1]+nv-1)], x[indices2[indices], (ci[1]+nv):ci[2]])
      xOldNew[, ci[1]:ci[2]] <- rbind(xOld[indices1[indices], ci[1]:(ci[1]+nv-1)], xOld[indices2[indices], (ci[1]+nv):ci[2]])
      }
    x <- xNew
    xOld <- xOldNew
    lZ <- lZNew
    nv <- nvNew
    W <- matrix(1/Nparticles, nrow = Nparticles, ncol = nodes)
  }
  return(cbind(x, lZ))
}
