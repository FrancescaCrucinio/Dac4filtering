dac_lgssm_lightweight_crossover_likelihood_tempering <- function(history, obs, tau, lambda, sigmaY, M = NULL){
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
  beta <- 2/(2-2^(-nlevels))
  for (i in 1:nchild^nlevels){
    if (i == 1) {
      x[, i] <- 0.5*history[, i, 1] + rnorm(Nparticles)/sqrt(tau)
    } else {
      x[, i] <- 0.5*tau*history[, i, 1]/(tau+lambda) + rnorm(Nparticles)/sqrt(tau+lambda)
    }
    lW[, i] <- -beta*0.5*(obs[i] - x[, i])^2/sigmaY
    max.lW <- max(lW[, i])
    W[, i] <- exp(lW[, i] - max.lW)
    W[, i] <- W[, i]/sum(W[, i])
  }
  # loop over tree levels excluding leaves
  for (u in 1:nlevels){
    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    nodes_dimension <- nchild^(u-1)
    betaNew <- (2-2^(-u))/(2-2^(-nlevels))

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
        out <- lgssm_adaptive_light_likelihood_tempering(Nparticles, i, u, nv, ci, lW, Nparticles, lambda, tau, x, xOld, betaNew - beta, obs)
        # update after mixture resampling
        indices <- out$resampled_indices
        xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv-1)], x[indices[, 2], (ci[1]+nv):ci[2]])
        historyIndexNew[, , i] <- historyIndexNew[indices[, 1], , i]

        # # tempering
        # tempering_out <- lgssm_tempering_crossover(ci, i, nv, lambda, tau, sigmaY, obs, xNew, history[, , 1], historyIndex,
        #                                  historyIndexNew, out$resampled_particles_lW, Nparticles, 1 - 1e-05, 1/nodes_dimension)
        # # update particles
        # xNew <- tempering_out$x
        # # update history
        # historyIndexNew <- tempering_out$history_index_updated
      }
      else{
        xOld <- history[historyIndex[, ci[1]+nv, nchild*i], ci[1]+nv, 1]
        out <- lgssm_light(i, u, nv, ci, W, Nparticles, M, lambda, tau, x, xOld)
        indices <- out$resampled_indices
        xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv-1)], x[indices[, 2], (ci[1]+nv):ci[2]])
        historyIndexNew[, , i] <- historyIndexNew[indices[, 1], , i]
      }
    }
    x <- xNew
    nv <- nvNew
    historyIndex <- historyIndexNew
    beta <- betaNew
  }
  return(x)
}

# xOldv is the ci[1]+nv component of xOld (main code)
lgssm_adaptive_light_likelihood_tempering <- function(ess_target, i, u, nv, ci, lW, Nparticles, lambda, tau, x, xOldv, beta_diff, obs){
  # binary tree
  nchild <- 2
  # mixture weights
  if(u == 1){
    lWmix <- lW[, (nchild*(i-1)+1)] + lW[, i*nchild] -
      0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOldv/(tau+lambda)) -
                      beta_diff*0.5*sweep(x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]])^2/sigmaY
      )
  } else{
    lWmix <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOldv/(tau+lambda)) -
                              beta_diff*0.5*sweep(x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]])^2/sigmaY
    )
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # build ESS
  ess_s <- sum(Wmix)
  ess_ss <- sum(Wmix^2)
  ess <- ess_s^2/ess_ss
  # first permutation
  permutation <- 1:Nparticles
  m <- 1
  while (ess < ess_target & m <= ceiling(sqrt(Nparticles))) {
    m <- m+1
    new_perm <- sample.int(Nparticles)
    if(u == 1){
      lWmix_perm <- lW[, (nchild*(i-1)+1)] + lW[new_perm, i*nchild] -
        0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[new_perm, (ci[1]+nv)] - 0.5*tau*xOldv[new_perm])
        )  -
        beta_diff*0.5*sweep(x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]])^2/sigmaY
    } else{
      lWmix_perm <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[new_perm, (ci[1]+nv)] - 0.5*tau*xOldv[new_perm])
      )  -
        beta_diff*0.5*sweep(x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]])^2/sigmaY
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
  # write.table(data.frame("u" = u, "m" = m), file = "data/adaptive_lgssm.csv", sep = ",", append = TRUE, quote = FALSE,
  #             col.names = FALSE, row.names = FALSE)
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(rep(1:Nparticles, times = m)[indices], permutation[indices]),
              "target_reached" = (ess >= ess_target), "resampled_particles_lW" = lWmix[indices]))
}
