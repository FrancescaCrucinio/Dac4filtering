lgssm_tempering_crossover <- function(ci, i, nv, lambda, tau, sigmaY, obs, x, xOld, historyIndex, historyIndexNew, after_mix_lW, ess_target, ess_decay_threshold, mcmc_sd){
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
      after_mix_lW <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOldv))
      if(nv == 1){ # children are leaves
        after_mix_lW <- after_mix_lW - 0.5*(obs[ci[1]] - x[, ci[1]])^2/sigmaY - 0.5*(obs[ci[2]] - x[, ci[2]])^2/sigmaY
      }
    }
    updated_particles <- lgssm_mcmc_move(x, ci, i, nv, sigmaY, tau, lambda, xOldv, xOld, historyIndex, obs, after_mix_lW, mcmc_sd, new_alpha)
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
    after_mix_lW <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOldv))
    if(nv == 1){ # children are leaves
      after_mix_lW <- after_mix_lW - 0.5*(obs[ci[1]] - x[, ci[1]])^2/sigmaY - 0.5*(obs[ci[2]] - x[, ci[2]])^2/sigmaY
    }
  }
  updated_particles <- lgssm_mcmc_move(x, ci, i, nv, sigmaY, tau, lambda, xOldv, xOld, historyIndex, obs, after_mix_lW, mcmc_sd, new_alpha)
  x <- updated_particles$x
  return(list("x" = x, "history_index_updated" = historyIndexNew, "lWmix" = after_mix_lW))
}


lgssm_mcmc_move <- function(x, ci, i, nv, sigmaY, tau, lambda, xOldv, xOld, historyIndex, obs, after_mix_lW, mcmc_sd, new_alpha){
  Nparticles <- nrow(x)
  nchild <- 2
  propose_x <- x
  propose_x[, ci[1]:ci[2]] <- x[, ci[1]:ci[2]] + mcmc_sd*matrix(rnorm(Nparticles*length(ci[1]:ci[2])), nrow = Nparticles)
  # mixture weight for proposed particle
  after_mix_lW_new <- -0.5*lambda * (lambda *propose_x[, (ci[1]+nv-1)]^2 - 2*propose_x[, (ci[1]+nv-1)] * (propose_x[, (ci[1]+nv)] - 0.5*tau*xOldv))
  # mixture weights ratio
  r2 <- after_mix_lW_new - after_mix_lW
  if(nv == 1){# children are leaves
    child_weights_proposal <- - 0.5*(obs[ci[1]] - propose_x[, ci[1]])^2/sigmaY - 0.5*(obs[ci[2]] - propose_x[, ci[2]])^2/sigmaY
    child_weights_current <- - 0.5*(obs[ci[1]] - x[, ci[1]])^2/sigmaY - 0.5*(obs[ci[2]] - x[, ci[2]])^2/sigmaY
    r2 <- r2 + child_weights_proposal - child_weights_current
  }
  # observation ratio
  r_obs <- rowSums((sweep(x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]]))^2-(sweep(propose_x[, ci[1]:ci[2]], 2, obs[ci[1]:ci[2]]))^2)/(2*sigmaY)
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



# for (n in 1:Nparticles) {
#   propose_x <- x[n, ]
#   propose_x[ci[1]:ci[2]] <- x[n, ci[1]:ci[2]] + mcmc_sd*rnorm(length(ci[1]:ci[2]))
#   # mixture weight for proposed particle
#   after_mix_lW_new <- -0.5*lambda * (lambda *propose_x[(ci[1]+nv-1)]^2 - 2*propose_x[(ci[1]+nv-1)] * (propose_x[(ci[1]+nv)] - 0.5*tau*xOldv[n]))
#   # mixture weights ratio
#   r2 <- after_mix_lW_new - after_mix_lW
#   if(nv == 1){# children are leaves
#     child_weights_proposal <- - 0.5*(obs[ci[1]] - propose_x[ci[1]])^2/sigmaY - 0.5*(obs[ci[2]] - propose_x[ci[2]])^2/sigmaY
#     child_weights_current <- - 0.5*(obs[ci[1]] - x[n, ci[1]])^2/sigmaY - 0.5*(obs[ci[2]] - x[n, ci[2]])^2/sigmaY
#     r2 <- r2 + child_weights_proposal - child_weights_current
#   }
#   # observation ratio
#   r_obs <- sum((x[n, ci[1]:ci[2]] - obs[ci[1]:ci[2]])^2-(propose_x[ci[1]:ci[2]] - obs[ci[1]:ci[2]])^2)/(2*sigmaY)
#   # product of children ratio
#
#   right_ancestor_coordinates <- cbind(historyIndex[n, (ci[1]+nv):ci[2], nchild*i], (ci[1]+nv):ci[2])
#   left_ancestor_coordinates <- cbind(historyIndex[n, (ci[1]+1):(ci[1]+nv-1), nchild*(i-1)+1], (ci[1]+1):ci[2])
#   # ratio for children nodes
#   if(nv == 1){# children are leaves
#     if(ci[1] == 1){
#       r1 <- 0.5*tau*(x[n, 1] - 0.5*xOld[historyIndex[n, 1, nchild*(i-1)+1], 1])^2 -
#         0.5*tau*(propose_x[1] - 0.5*xOld[historyIndex[n, 1, nchild*(i-1)+1], 1])^2
#     }
#     else {
#       r1 <- 0.5*(tau+lambda)*(x[n, ci[1]] - 0.5*tau*xOld[historyIndex[n, ci[1], nchild*(i-1)+1], ci[1]]/(tau+lambda))^2 -
#         0.5*(tau+lambda)*(propose_x[ci[1]] - 0.5*tau*xOld[historyIndex[n, ci[1], nchild*(i-1)+1], ci[1]]/(tau+lambda))^2
#     }
#     r1 <- r1 + 0.5*(tau+lambda)*(x[n, ci[2]] - 0.5*tau*xOld[historyIndex[n, ci[2], nchild*i], ci[2]]/(tau+lambda))^2 -
#       0.5*(tau+lambda)*(propose_x[ci[2]] - 0.5*tau*xOld[historyIndex[n, ci[2], nchild*i], ci[2]]/(tau+lambda))^2
#   } else {
#     if(ci[1] == 1){
#       r1 <- 0.5*tau*(x[n, 1] - 0.5*xOld[historyIndex[n, 1, nchild*(i-1)+1], 1])^2 -
#         0.5*tau*(propose_x[1] - 0.5*xOld[historyIndex[n, 1, nchild*(i-1)+1], 1])^2
#     }
#     else {
#       r1 <- 0.5*(tau+lambda)*(xNew[n, ci[1]] - 0.5*tau*xOld[historyIndex[n, ci[1], nchild*(i-1)+1], ci[1]]/(tau+lambda))^2 -
#         (propose_x[ci[1]] - 0.5*tau*xOld[historyIndex[n, ci[1], nchild*(i-1)+1], ci[1]]/(tau+lambda))^2
#     }
#     r1 <- r1 + 0.5*(tau+lambda)*sum((x[n, (ci[1]+1):(ci[1]+nv-1)] - (0.5*tau*xOld[left_ancestor_coordinates] + lambda*x[n, ci[1]:(ci[1]+nv-2)])/(tau+lambda))^2 -
#                                       (propose_x[(ci[1]+1):(ci[1]+nv-1)] - (0.5*tau*xOld[left_ancestor_coordinates] + lambda*propose_x[ci[1]:(ci[1]+nv-2)])/(tau+lambda))^2) +
#       0.5*(tau+lambda)*sum((x[n, (ci[1]+nv+1):ci[2]] - (0.5*tau*xOld[right_ancestor_coordinates] + lambda*x[n, (ci[1]+nv-1):(ci[2]-1)])/(tau+lambda))^2 -
#                              (propose_x[, (ci[1]+nv+1):ci[2]] - (0.5*tau*xOld[right_ancestor_coordinates] + lambda*propose_x[, (ci[1]+nv-1):(ci[2]-1)])/(tau+lambda))^2)
#   }
#
#   # accept/reject
#   mh_ratio <- r1 + r_obs + new_alpha*r2
#   if(runif(1) <= exp(mh_ratio)){
#     x[n, ] <- propose_x
#   }
# }
