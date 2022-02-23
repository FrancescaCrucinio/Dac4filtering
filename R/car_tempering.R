car_tempering <- function(ci, i, nv, nvNew, sigmaX, sigmaY, obs, x, xOld, historyIndex, historyIndexNew, after_mix_lW, ess_target, ess_decay_threshold, mcmc_sd){
  # tree topology
  nchild <- 2
  Nparticles <- length(after_mix_lW)
  # alpha after mixture resampling
  current_alpha <- bisection_ess(after_mix_lW, ess_target)
  lOmega <- rep(0, times = Nparticles)
  same_sign <- FALSE
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
      after_mix_lW <- rep(0, times = Nparticles)
      for (n in 1:Nparticles) {
        # get ancestor
        right_ancestor_coordinates <- cbind(historyIndex[n, (ci[1]+nv):d, nchild*i], (ci[1]+nv):d)
        right_ancestor <- xOld[right_ancestor_coordinates]
        # merge the two children nodes
        mx <- x[n, ci[1]:ci[2]]
        after_mix_lW[n] <- sum(x[n, ci[1]:(ci[1]+nv-1)])*sum(x[n, (ci[1]+nv):ci[2]])/(d*sigmaX) -
          (sum(cumsum(mx[1:(nvNew-1)]/d)^2) - sum(cumsum(x[n, ci[1]:(ci[1]+nv-1)][seq(length.out = (nv-1))]/d)^2) -
             sum(cumsum(x[n, (ci[1]+nv):ci[2]][seq(length.out = (nv-1))]/d)^2))/(2*sigmaX) -
          sum(x[n, ci[1]:(ci[1]+nv-1)])*sum(cumsum(rev(right_ancestor)))/(d^2*sigmaX)
      }
    }
    updated_particles <- car_mcmc_move(x, ci, i, nv, nvNew, sigmaY, sigmaX, xOld, historyIndex, obs, after_mix_lW, mcmc_sd, new_alpha)
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
    after_mix_lW <- rep(0, times = Nparticles)
    for (n in 1:Nparticles) {
      # get ancestor
      right_ancestor_coordinates <- cbind(historyIndex[n, (ci[1]+nv):d, nchild*i], (ci[1]+nv):d)
      right_ancestor <- xOld[right_ancestor_coordinates]
      # merge the two children nodes
      mx <- x[n, ci[1]:ci[2]]
      after_mix_lW[n] <- sum(x[n, ci[1]:(ci[1]+nv-1)])*sum(x[n, (ci[1]+nv):ci[2]])/(d*sigmaX) -
        (sum(cumsum(mx[1:(nvNew-1)]/d)^2) - sum(cumsum(x[n, ci[1]:(ci[1]+nv-1)][seq(length.out = (nv-1))]/d)^2) -
           sum(cumsum(x[n, (ci[1]+nv):ci[2]][seq(length.out = (nv-1))]/d)^2))/(2*sigmaX) -
        sum(x[n, ci[1]:(ci[1]+nv-1)])*sum(cumsum(rev(right_ancestor)))/(d^2*sigmaX)
    }
  }
  updated_particles <- car_mcmc_move(x, ci, i, nv, nvNew, sigmaY, sigmaX, xOld, historyIndex, obs, after_mix_lW, mcmc_sd, new_alpha)
  x <- updated_particles$x
  return(list("x" = x, "history_index_updated" = historyIndexNew, "lWmix" = after_mix_lW))
}


car_mcmc_move <- function(x, ci, i, nv, nvNew, sigmaY, sigmaX, xOld, historyIndex, obs, after_mix_lW, mcmc_sd, new_alpha){
  Nparticles <- nrow(x)
  nchild <- 2
  propose_x <- x
  propose_x[, ci[1]:ci[2]] <- x[, ci[1]:ci[2]] + mcmc_sd*matrix(rnorm(Nparticles*length(ci[1]:ci[2])), nrow = Nparticles)
  # mixture weight for proposed particle
  after_mix_lW_new <- rep(0, times = Nparticles)
  r1 <- rep(0, times = Nparticles)
  for (n in 1:Nparticles) {
    # get ancestors
    right_ancestor_coordinates <- cbind(historyIndex[n, (ci[1]+nv):d, nchild*i], (ci[1]+nv):d)
    right_ancestor <- xOld[right_ancestor_coordinates]
    left_ancestor_coordinates <- cbind(historyIndex[n, ci[1]:d, nchild*(i-1)+1], ci[1]:d)
    left_ancestor <- xOld[left_ancestor_coordinates]
    # merge the two children nodes
    propose_mx <- propose_x[n, ci[1]:ci[2]]
    mx <- x[n, ci[1]:ci[2]]
    after_mix_lW_new[n] <- sum(propose_x[n, ci[1]:(ci[1]+nv-1)])*sum(propose_x[n, (ci[1]+nv):ci[2]])/(d*sigmaX) -
      (sum(cumsum(propose_mx[1:(nvNew-1)]/d)^2) - sum(cumsum(propose_x[n, ci[1]:(ci[1]+nv-1)][seq(length.out = (nv-1))]/d)^2) -
         sum(cumsum(propose_x[n, (ci[1]+nv):ci[2]][seq(length.out = (nv-1))]/d)^2))/(2*sigmaX) -
      sum(propose_x[n, ci[1]:(ci[1]+nv-1)])*sum(cumsum(rev(right_ancestor)))/(d^2*sigmaX)
    # product of children ratio
    r1[n] <- 0.5*sum(-(propose_x[n, ci[1]:(ci[1]+nv-1)] - (cumsum(propose_x[n, ci[1]:(ci[1]+nv-1)][seq(length.out = (nv-1))] + cumsum(rev(left_ancestor)))/d))^2 +
                               (x[n, ci[1]:(ci[1]+nv-1)] - (cumsum(x[n, ci[1]:(ci[1]+nv-1)][seq(length.out = (nv-1))] + cumsum(rev(left_ancestor)))/d))^2)/sigmaX +
      0.5*sum(-(propose_x[n, (ci[1]+nv):ci[2]] - (cumsum(x[n, (ci[1]+nv):ci[2]][seq(length.out = (nv-1))] + cumsum(rev(right_ancestor)))/d))^2 +
                       (x[n, (ci[1]+nv):ci[2]] - (cumsum(x[n, (ci[1]+nv):ci[2]][seq(length.out = (nv-1))] + cumsum(rev(right_ancestor)))/d))^2)
  }
  # mixture weights ratio
  r2 <- after_mix_lW_new - after_mix_lW
  if(nv == 1){# children are leaves
    child_weights_proposal <- - 0.5*(obs[ci[1]] - propose_x[, ci[1]])^2/sigmaY - 0.5*(obs[ci[2]] - propose_x[, ci[2]])^2/sigmaY
    child_weights_current <- - 0.5*(obs[ci[1]] - x[, ci[1]])^2/sigmaY - 0.5*(obs[ci[2]] - x[, ci[2]])^2/sigmaY
    r2 <- r2 + child_weights_proposal - child_weights_current
  }
  # observation ratio
  r_obs <- rowSums((sweep(x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]]))^2-(sweep(propose_x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]]))^2)/(2*sigmaY)
  mh_ratio <- r1 + r_obs + new_alpha*r2
  accepted <- runif(Nparticles) <= exp(mh_ratio)
  print(paste(sum(accepted)/Nparticles))
  x[accepted , ] <- propose_x[accepted, ]
  after_mix_lW[accepted] <- after_mix_lW_new[accepted]
  return(list("x" = x, "after_mix_lW" = after_mix_lW))
}


