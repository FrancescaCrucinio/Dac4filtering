dac_lgssm <- function(xOld, obs, tau, lambda, sigmaY){
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
              0.5*lambda * (lambda*x[n1, (ci[1]+nv-1)]^2/(tau+lambda) -
                2*x[n1, (ci[1]+nv-1)] * (x[n2, (ci[1]+nv)] - 0.5*tau*xOld[n2, (ci[1]+nv)]/(tau+lambda))
                )
          }
        }
        max.lWmix <- max(lWmix)
        Wmix <- exp(lWmix - max.lWmix)
      } else {
        for (n1 in 1:Nparticles) {
          for (n2 in 1:Nparticles) {
            lWmix[n1, n2] <- -0.5*lambda * (lambda*x[n1, (ci[1]+nv-1)]^2/(tau+lambda) -
                              2*x[n1, (ci[1]+nv-1)] * (x[n2, (ci[1]+nv)] - 0.5*tau*xOld[n2, (ci[1]+nv)]/(tau+lambda))
                              )
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
        # update xOld
        xOld[n, ci[1]:ci[2]] <- xOld[res[n, 1], ci[1]:ci[2]]
      }
    }
    x <- xNew
    nv <- nvNew
  }
  return(x)
}

dac_lgssm_lightweight <- function(xOld, obs, tau, lambda, sigmaY, M = NULL){
  if(is.null(M)) {
    # number of samples for lightweight mixture (no adaptation)
    M <- ceiling(sqrt(Nparticles))
  }
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
    if (i == 1) {
      x[, i] <- 0.5*xOld[, i] + rnorm(Nparticles)/sqrt(tau)
    } else {
      x[, i] <- 0.5*tau*xOld[, i]/(tau+lambda) + rnorm(Nparticles)/sqrt(tau+lambda)
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
      if(M == "adaptive") {
        out <- lgssm_adaptive_light(Nparticles, i, u, nv, ci, lW, Nparticles, lambda, tau, x, xOld[, ci[1]+nv])
        # update after mixture resampling
        indices <- out$resampled_indices
        xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv-1)], x[indices[, 2], (ci[1]+nv):ci[2]])
        xOld[, ci[1]:ci[2]] <- xOld[indices[, 1], ci[1]:ci[2]]
        # tempering
        # if(!out$target_reached){
        #   historyIndex <- array(1:Nparticles, dim = c(Nparticles, d, d))
        #   tempering_out <- lgssm_tempering
        #   # update particles
        #   xNew <- tempering_out$x
        #   # update history
        #   xOld <- tempering_out$xOld
        # }
      }
      else{
        indices <- lgssm_light(i, u, nv, ci, W, Nparticles, M, lambda, tau, x, xOld[, ci[1]+nv])
        xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv-1)], x[indices[, 2], (ci[1]+nv):ci[2]])
        xOld[, ci[1]:ci[2]] <- xOld[indices[, 1], ci[1]:ci[2]]
      }
      }
    x <- xNew
    nv <- nvNew
  }
  return(x)
}

dac_lgssm_lc <- function(xOld, obs, tau, lambda, sigmaY){
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
    xOldNew <- matrix(0, nrow = Nparticles, ncol = d)
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
                             2*x[indices1, (ci[1]+nv-1)] * (x[indices2, (ci[1]+nv)] - 0.5*tau*xOld[indices2, (ci[1]+nv)]/(tau+lambda))
                          )
      max.lW <- max(lW)
      WNew[, i] <- exp(lW - max.lW)
      # update particles
      xNew[, ci[1]:ci[2]] <- cbind(x[indices1, ci[1]:(ci[1]+nv-1)], x[indices2, (ci[1]+nv):ci[2]])
      # update xOld
      xOld[, ci[1]:ci[2]] <- xOld[indices1, ci[1]:ci[2]]
    }
    x <- xNew
    nv <- nvNew
    W <- WNew
  }
  return(x)
}
