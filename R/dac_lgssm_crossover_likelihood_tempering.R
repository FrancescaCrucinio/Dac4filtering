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
  beta <- 1/(2-2^(-nlevels))
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
  # saveRDS(list("x" = x, "W"=W), file = paste0("/Users/francescacrucinio/Documents/Dac4filtering/test/lt_node0.rds"))
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
        if(!out$target_reached){
        # tempering
        tempering_out <- lgssm_tempering_crossover_likelihood_tempering(ci, i, nv, lambda, tau, sigmaY, obs, x, history[, , 1], historyIndex, historyIndexNew,
                                                      out$resampled_particles_lW, Nparticles, 1-1e-05, 1/nodes_dimension, betaNew - beta, beta)
        # update particles
        xNew <- tempering_out$x
        # update history
        historyIndexNew <- tempering_out$history_index_updated
        } else {
          # mcmc move
          updated_particles <- lgssm_mcmc_move_likelihood_tempering(x, ci, i, nv, sigmaY, tau, lambda,
                                                  history[historyIndex[, ci[1]+nv, nchild*i], ci[1]+nv, 1], history[, , 1],
                                                  historyIndex, obs, out$resampled_particles_lW, 1/nodes_dimension,
                                                  1, betaNew - beta, beta)
          x <- updated_particles$x
        }
      }
      else{
        xOld <- history[historyIndex[, ci[1]+nv, nchild*i], ci[1]+nv, 1]
        out <- lgssm_light(i, u, nv, ci, W, Nparticles, M, lambda, tau, x, xOld)
        indices <- out$resampled_indices
        xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv-1)], x[indices[, 2], (ci[1]+nv):ci[2]])
        historyIndexNew[, , i] <- historyIndexNew[indices[, 1], , i]
      }
    }
    # saveRDS(x, file = paste0("/Users/francescacrucinio/Documents/Dac4filtering/test/lt_node", u, ".rds"))
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
      0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOldv/(tau+lambda))) -
                      rowSums(beta_diff*0.5*sweep(x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]])^2)/sigmaY
  } else{
    lWmix <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOldv/(tau+lambda))) -
                              rowSums(beta_diff*0.5*sweep(x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]])^2)/sigmaY
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
        ) - rowSums(beta_diff*0.5*sweep(x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]])^2)/sigmaY
    } else{
      lWmix_perm <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[new_perm, (ci[1]+nv)] - 0.5*tau*xOldv[new_perm])
      ) - rowSums(beta_diff*0.5*sweep(x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]])^2)/sigmaY
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
  # write.table(data.frame("u" = u, "m" = m), file = "data/adaptive_lgssm_cv_lt.csv", sep = ",", append = TRUE, quote = FALSE,
  # col.names = FALSE, row.names = FALSE)
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(rep(1:Nparticles, times = m)[indices], permutation[indices]),
              "target_reached" = (ess >= ess_target), "resampled_particles_lW" = lWmix[indices]))
}


lgssm_tempering_crossover_likelihood_tempering <- function(ci, i, nv, lambda, tau, sigmaY, obs, x, xOld, historyIndex, historyIndexNew,
                                                           after_mix_lW, ess_target, ess_decay_threshold, mcmc_sd, beta_diff, beta){
  # tree topology
  nchild <- 2
  Nparticles <- length(after_mix_lW)
  # alpha after mixture resampling
  current_alpha <- bisection_ess(after_mix_lW, ess_target)
  lOmega <- rep(0, times = Nparticles)
  same_sign <- FALSE
  # extract history for mixture weights
  xOldv <- xOld[historyIndex[, ci[1]+nv, nchild*i], ci[1]+nv]
  while (current_alpha < 1 & !same_sign) {
    res_bisection <- bisection_cess(lOmega, after_mix_lW, current_alpha, ess_decay_threshold)
    new_alpha <- res_bisection$alpha_star
    same_sign <- res_bisection$same_sign
    newlOmega <- lOmega + (new_alpha - current_alpha)*after_mix_lW
    omega <- exp(newlOmega - max(newlOmega))
    ess <- sum(omega)^2/sum(omega^2)
    if(ess < Nparticles / 2){
      resampled_indices <- stratified_resample(omega/sum(omega), Nparticles)
      newlOmega <- rep(0, times = Nparticles)
      x[, ci[1]:ci[2]] <- x[resampled_indices, ci[1]:ci[2]]
      historyIndex[, , i] <- historyIndex[resampled_indices, , i]
      historyIndexNew[, , i] <- historyIndexNew[resampled_indices, , i]
      after_mix_lW <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOldv)) -
        rowSums(beta_diff*0.5*sweep(x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]])^2)/sigmaY
    }
    updated_particles <- lgssm_mcmc_move_likelihood_tempering(x, ci, i, nv,
                                sigmaY, tau, lambda, xOldv, xOld, historyIndex, obs, after_mix_lW, mcmc_sd, new_alpha, beta_diff, beta)
    x <- updated_particles$x
    current_alpha <- new_alpha
    lOmega <- newlOmega
  }
  # alpha = 1
  new_alpha <- 1
  newlOmega <- lOmega + (1 - current_alpha)*after_mix_lW
  omega <- exp(newlOmega - max(newlOmega))
  ess <- sum(omega)^2/sum(omega^2)
  if(ess < Nparticles / 2){
    resampled_indices <- stratified_resample(omega/sum(omega), Nparticles)
    newlOmega <- rep(0, times = Nparticles)
    x[, ci[1]:ci[2]] <- x[resampled_indices, ci[1]:ci[2]]
    historyIndex[, , i] <- historyIndex[resampled_indices, , i]
    historyIndexNew[, , i] <- historyIndexNew[resampled_indices, , i]
    after_mix_lW <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOldv)) -
      rowSums(beta_diff*0.5*sweep(x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]])^2)/sigmaY
  }
  updated_particles <- lgssm_mcmc_move_likelihood_tempering(x, ci, i, nv,
                            sigmaY, tau, lambda, xOldv, xOld, historyIndex, obs, after_mix_lW, mcmc_sd, new_alpha, beta_diff, beta)
  x <- updated_particles$x
  return(list("x" = x, "history_index_updated" = historyIndexNew, "lWmix" = after_mix_lW))
}


lgssm_mcmc_move_likelihood_tempering <- function(x, ci, i, nv,
            sigmaY, tau, lambda, xOldv, xOld, historyIndex, obs, after_mix_lW, mcmc_sd, new_alpha, beta_diff, beta){
  Nparticles <- nrow(x)
  nchild <- 2
  propose_x <- x
  propose_x[, ci[1]:ci[2]] <- x[, ci[1]:ci[2]] + mcmc_sd*matrix(rnorm(Nparticles*length(ci[1]:ci[2])), nrow = Nparticles)
  # mixture weight for proposed particle
  after_mix_lW_new <- -0.5*lambda * (lambda *propose_x[, (ci[1]+nv-1)]^2 - 2*propose_x[, (ci[1]+nv-1)] * (propose_x[, (ci[1]+nv)] - 0.5*tau*xOldv)) -
    rowSums(beta_diff*0.5*sweep(propose_x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]])^2)/sigmaY
  # mixture weights ratio
  r2 <- after_mix_lW_new - after_mix_lW
  if(nv == 1){# children are leaves
    child_weights_proposal <- - beta*0.5*(obs[ci[1]] - propose_x[, ci[1]])^2/sigmaY - beta*0.5*(obs[ci[2]] - propose_x[, ci[2]])^2/sigmaY
    child_weights_current <- - beta*0.5*(obs[ci[1]] - x[, ci[1]])^2/sigmaY - beta*0.5*(obs[ci[2]] - x[, ci[2]])^2/sigmaY
    r2 <- r2 + child_weights_proposal - child_weights_current
  }
  # observation ratio
  r_obs <- rowSums((sweep(x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]]))^2-(sweep(propose_x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]]))^2)*beta/(2*sigmaY)
  # product of children ratio
  if(ci[1] == 1){
    r1 <- 0.5*tau*(x[, 1] - 0.5*xOld[historyIndex[, 1, nchild*(i-1)+1], 1])^2 -
      0.5*tau*(propose_x[, 1] - 0.5*xOld[historyIndex[, 1, nchild*(i-1)+1], 1])^2
  }
  else {
    r1 <- 0.5*(tau+lambda)*(x[, ci[1]] - 0.5*tau*xOld[historyIndex[, ci[1], nchild*(i-1)+1], ci[1]]/(tau+lambda))^2 -
      0.5*(tau+lambda)*(propose_x[, ci[1]] - 0.5*tau*xOld[historyIndex[, ci[1], nchild*(i-1)+1], ci[1]]/(tau+lambda))^2
  }
  if(nv == 1){# children are leaves
    r1 <- r1 + 0.5*(tau+lambda)*(x[, ci[2]] - 0.5*tau*xOld[historyIndex[, ci[2], nchild*i], ci[2]]/(tau+lambda))^2 -
      0.5*(tau+lambda)*(propose_x[, ci[2]] - 0.5*tau*xOld[historyIndex[, ci[2], nchild*i], ci[2]]/(tau+lambda))^2
  } else {
    right_ancestor_coordinates <- cbind(c(historyIndex[, (ci[1]+nv+1):ci[2], nchild*i]), rep((ci[1]+nv+1):ci[2], each = Nparticles))
    right_ancestor <- matrix(xOld[right_ancestor_coordinates, drop = FALSE], nrow = Nparticles)
    left_ancestor_coordinates <- cbind(c(historyIndex[, (ci[1]+1):(ci[1]+nv-1), nchild*(i-1)+1]), rep((ci[1]+1):(ci[1]+nv-1), each = Nparticles))
    left_ancestor <- matrix(xOld[left_ancestor_coordinates, drop = FALSE], nrow = Nparticles)
    r1 <- r1 + 0.5*(tau+lambda)*rowSums((x[, (ci[1]+1):(ci[1]+nv-1), drop = FALSE] - (0.5*tau*left_ancestor + lambda*x[, ci[1]:(ci[1]+nv-2), drop = FALSE])/(tau+lambda))^2 -
                                          (propose_x[, (ci[1]+1):(ci[1]+nv-1), drop = FALSE] - (0.5*tau*left_ancestor + lambda*propose_x[, ci[1]:(ci[1]+nv-2), drop = FALSE])/(tau+lambda))^2) +
      0.5*(tau+lambda)*rowSums((x[, (ci[1]+nv+1):ci[2], drop = FALSE] - (0.5*tau*right_ancestor + lambda*x[, (ci[1]+nv):(ci[2]-1), drop = FALSE])/(tau+lambda))^2 -
                                 (propose_x[, (ci[1]+nv+1):ci[2], drop = FALSE] - (0.5*tau*right_ancestor + lambda*propose_x[, (ci[1]+nv):(ci[2]-1), drop = FALSE])/(tau+lambda))^2)
  }
  mh_ratio <- r1 + r_obs + new_alpha*r2
  accepted <- runif(Nparticles) <= exp(mh_ratio)
  # print(paste(sum(accepted)/Nparticles))
  x[accepted , ] <- propose_x[accepted, ]
  after_mix_lW[accepted] <- after_mix_lW_new[accepted]
  return(list("x" = x, "after_mix_lW" = after_mix_lW))
}
