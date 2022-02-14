dac_lgssm_lc_crossover <- function(history, obs, tau, lambda, sigmaY){
  # dimension and number of particles
  d <- ncol(history[, , 1])
  Nparticles <- nrow(history[, , 1])
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # leaves
  # number of variables
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  W <- matrix(0, nrow = Nparticles, ncol = d)
  # history indices
  historyIndex <- array(1:Nparticles, dim = c(Nparticles, d, d))
  for (i in 1:nchild^nlevels){
    if (i == 1) {
      x[, i] <- 0.5*history[, i, 1] + rnorm(Nparticles)/sqrt(tau)
    } else {
      x[, i] <- 0.5*tau*history[, i, 1]/(tau+lambda) + rnorm(Nparticles)/sqrt(tau+lambda)
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
    # updated history
    historyIndexNew <- array(0, dim = c(Nparticles, d, nodes))

    for (i in 1:nodes){
      # get children indices
      ci <- child_indices(i, nvNew)
      # resample on each children
      # child 1
      indices1 <- stratified_resample(W[, nchild*(i-1)+1], Nparticles)
      # child 2 (with random permutation)
      indices2 <- sample(stratified_resample(W[, nchild*i], Nparticles))
      # update history
      historyIndex[, , nchild*(i-1)+1] <- historyIndex[indices1, , nchild*(i-1)+1]
      historyIndex[, , nchild*i] <- historyIndex[indices2, , nchild*i]
      # update particles
      x[, ci[1]:ci[2]] <- cbind(x[indices1, ci[1]:(ci[1]+nv-1)], x[indices2, (ci[1]+nv):ci[2]])
      if(u > 1){
        # mutation
        historyIndexNew[, , i] <- crossover(i, nodes, x, history, historyIndex, tau, lambda)
      }
      else{ # at the leaf level all histories are the same
        historyIndexNew[, , i] <- historyIndex[, , i]
      }
      # weights
      lW <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) -
                             2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*history[historyIndex[, (ci[1]+nv), nchild*i], (ci[1]+nv), 1]/(tau+lambda))
      )
      max.lW <- max(lW)
      W[, i] <- exp(lW - max.lW)
    }
    nv <- nvNew
    historyIndex <- historyIndexNew
  }
  return(x)
}

dac_lgssm_lightweight_crossover <- function(history, obs, tau, lambda, sigmaY, M = NULL){
  if(is.null(M)) {
    # number of samples for lightweight mixture (no adaptation)
    M <- ceiling(sqrt(Nparticles))
  }
  # dimension and number of particles
  d <- ncol(history[, , 1])
  Nparticles <- nrow(history[, , 1])
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # leaves
  # number of variables
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  lW <- matrix(0, nrow = Nparticles, ncol = d)
  W <- matrix(0, nrow = Nparticles, ncol = d)
  # history indices
  historyIndex <- array(1:Nparticles, dim = c(Nparticles, d, d))
  for (i in 1:nchild^nlevels){
    if (i == 1) {
      x[, i] <- 0.5*history[, i, 1] + rnorm(Nparticles)/sqrt(tau)
    } else {
      x[, i] <- 0.5*tau*history[, i, 1]/(tau+lambda) + rnorm(Nparticles)/sqrt(tau+lambda)
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

    # updated particles
    xNew <- matrix(0, nrow = Nparticles, ncol = d)
    # updated history
    historyIndexNew <- array(0, dim = c(Nparticles, d, nodes))
    for (i in 1:nodes){
      # get children indices
      ci <- child_indices(i, nvNew)
      if(u > 1){
        # mutation
        historyIndexNew[, , i] <- crossover(i, nodes, x, history, historyIndex, tau, lambda)
      }
      else{ # at the leaf level all histories are the same
        historyIndexNew[, , i] <- historyIndex[, , i]
      }
      # lightweight mixture resampling
      if(M == "adaptive") {
        # extract history for mixture weights
        xOld <- history[historyIndex[, ci[1]+nv, nchild*i], ci[1]+nv, 1]
        out <- lgssm_adaptive_light(Nparticles, i, u, nv, ci, lW, Nparticles, lambda, tau, x, xOld)
        # update after mixture resampling
        indices <- out$resampled_indices
        xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv-1)], x[indices[, 2], (ci[1]+nv):ci[2]])
        historyIndexNew[, , i] <- historyIndexNew[indices[, 1], , i]
        # tempering
        if(!out$target_reached){
          print(paste("tempering"))
          tempering_out <- lgssm_tempering(ci, nv, lambda, tau, xNew, history[, , 1], historyIndexNew,
                                           out$resampled_particles_lW, Nparticles, 1 - 1e-05, 1e-01)
          # update particles
          xNew <- tempering_out$x
          # update history
          historyIndexNew[, , i] <- historyIndexNew[tempering_out$resampled_indices, , i]
        }
      }
      else{
        xOld <- history[historyIndex[, ci[1]+nv, nchild*i], ci[1]+nv, 1]
        indices <- lgssm_light(i, u, nv, ci, W, Nparticles, M, lambda, tau, x, xOld)
      }
    }
    x <- xNew
    nv <- nvNew
    historyIndex <- historyIndexNew
  }
  return(x)
}

dac_lgssm_crossover <- function(history, obs, tau, lambda, sigmaY){
  # dimension and number of particles
  d <- ncol(history[, , 1])
  Nparticles <- nrow(history[, , 1])
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # leaves
  # number of variables
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  lW <- matrix(0, nrow = Nparticles, ncol = d)
  # history indices
  historyIndex <- array(1:Nparticles, dim = c(Nparticles, d, d))
  for (i in 1:nchild^nlevels){
    if (i == 1) {
      x[, i] <- 0.5*history[, i, 1] + rnorm(Nparticles)/sqrt(tau)
    } else {
      x[, i] <- 0.5*tau*history[, i, 1]/(tau+lambda) + rnorm(Nparticles)/sqrt(tau+lambda)
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

    # updated particles
    xNew <- matrix(0, nrow = Nparticles, ncol = d)
    # updated history
    historyIndexNew <- array(0, dim = c(Nparticles, d, nodes))
    for (i in 1:nodes){
      if(u > 1){
        # mutation
        historyIndexNew[, , i] <- crossover(i, nodes, x, history, historyIndex, tau, lambda)
      }
      else{ # at the leaf level all histories are the same
        historyIndexNew[, , i] <- historyIndex[, , i]
      }

      # mixture weights
      lWmix <- matrix(0, ncol = Nparticles, nrow = Nparticles)
      ci <- child_indices(i, nvNew)
      if(u == 1){
        for (n1 in 1:Nparticles) {
          for (n2 in 1:Nparticles) {
            lWmix[n1, n2] <- lW[n1, (nchild*(i-1)+1)] + lW[n2, i*nchild] -
              0.5*lambda * (lambda *x[n1, (ci[1]+nv-1)]^2/(tau+lambda) -
                              2*x[n1, (ci[1]+nv-1)] * (x[n2, (ci[1]+nv)] - 0.5*tau*history[historyIndex[n2, ci[1]+nv, nchild*i], ci[1]+nv, 1]/(tau+lambda))
              )
          }
        }
        max.lWmix <- max(lWmix)
        Wmix <- exp(lWmix - max.lWmix)
      } else {
        for (n1 in 1:Nparticles) {
          for (n2 in 1:Nparticles) {
            lWmix[n1, n2] <- -0.5*lambda * (lambda *x[n1, (ci[1]+nv-1)]^2 -
                                              2*x[n1, (ci[1]+nv-1)] * (x[n2, (ci[1]+nv)] - 0.5*tau*history[historyIndex[n2, ci[1]+nv, nchild*i], ci[1]+nv, 1])
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
        # update history
        historyIndexNew[n, , i] <- historyIndexNew[res[n, 1], , i]
      }
    }
    x <- xNew
    nv <- nvNew
    historyIndex <- historyIndexNew
  }
  return(x)
}
