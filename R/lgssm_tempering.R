lgssm_tempering <- function(ci, i, nv, lambda, tau, sigmaY, obs, x, xOld, historyIndex, historyIndexNew, after_mix_lW, ess_target, ess_decay_threshold, mcmc_sd){
  nchild <- 2

  Nparticles <- length(after_mix_lW)

  current_alpha <- bisection_ess(after_mix_lW, ess_target)
  lOmega <- rep(0, times = Nparticles)
  same_sign <- FALSE
  # extract history for mixture weights
  xOldv <- xOld[historyIndex[, ci[1]+nv, nchild*i], ci[1]+nv]
  while (abs(current_alpha - 1) > 1e-01 & !same_sign) {
    res_bisection <- bisection_cess(lOmega, after_mix_lW, current_alpha, ess_decay_threshold)
    new_alpha <- res_bisection$alpha_star
    same_sign <- res_bisection$same_sign
    newlOmega <- lOmega + (new_alpha - current_alpha)*after_mix_lW
    omega <- exp(newlOmega - max(newlOmega))
    ess <- sum(omega)^2/sum(omega^2)
    if(ess < Nparticles / 2){
      resampled_indices <- stratified_resample(omega/sum(omega), Nparticles)
      newlOmega <- rep(0, times = Nparticles)
      x[, ci[1]:ci[2]] <- cbind(x[resampled_indices, ci[1]:(ci[1]+nv-1)], x[resampled_indices, (ci[1]+nv):ci[2]])
      historyIndex[, , i] <- historyIndex[resampled_indices, , i]
      historyIndexNew[, , i] <- historyIndexNew[resampled_indices, , i]
      # update lW
      after_mix_lW <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOldv))
    }
    # MCMC move
    for (n in 1:Nparticles) {
      propose_x <- x[n, ]
      propose_x[ci[1]:ci[2]] <- x[n, ci[1]:ci[2]] + mcmc_sd*rnorm(length(ci[1]:ci[2]))
      right_ancestor_coordinates <- cbind(historyIndex[n, (ci[1]+nv):ci[2], nchild*i], (ci[1]+nv):ci[2])
      left_ancestor_coordinates <- cbind(historyIndex[n, (ci[1]+1):(ci[1]+nv-1), nchild*(i-1)+1], (ci[1]+1):ci[2])
      crossedover_ancestor_coordinates <- cbind(historyIndex[n, (ci[1]+1):ci[2], nchild*(i-1)+1], (ci[1]+1):ci[2])
      # observation ratio
      r_obs <- sum((x[n, ci[1]:ci[2]] - obs[ci[1]:ci[2]])^2-(propose_x[ci[1]:ci[2]] - obs[ci[1]:ci[2]])^2)/(2*sigmaY)
      # ratio for children nodes
      if(nv == 1){# children are leaves
        if(ci[1] == 1){
          r1 <- 0.5*tau*(x[n, 1] - 0.5*xOld[historyIndex[n, 1, nchild*(i-1)+1], 1])^2 -
                0.5*tau*(propose_x[1] - 0.5*xOld[historyIndex[n, 1, nchild*(i-1)+1], 1])^2
        }
        else {
          r1 <- 0.5*(tau+lambda)*(x[n, ci[1]] - 0.5*tau*xOld[historyIndex[n, ci[1], nchild*(i-1)+1], ci[1]]/(tau+lambda))^2 -
                0.5*(tau+lambda)*(propose_x[ci[1]] - 0.5*tau*xOld[historyIndex[n, ci[1], nchild*(i-1)+1], ci[1]]/(tau+lambda))^2
        }
        r1 <- r1 + 0.5*(tau+lambda)*(x[n, ci[2]] - 0.5*tau*xOld[historyIndex[n, ci[2], nchild*i], ci[2]]/(tau+lambda))^2 -
                  0.5*(tau+lambda)*(propose_x[ci[2]] - 0.5*tau*xOld[historyIndex[n, ci[2], nchild*i], ci[2]]/(tau+lambda))^2
      } else {
        if(ci[1] == 1){
          r1 <- 0.5*tau*(x[n, 1] - 0.5*xOld[historyIndex[n, 1, nchild*(i-1)+1], 1])^2 -
            0.5*tau*(propose_x[1] - 0.5*xOld[historyIndex[n, 1, nchild*(i-1)+1], 1])^2
        }
        else {
          r1 <- 0.5*(tau+lambda)*(xNew[n, ci[1]] - 0.5*tau*xOld[historyIndex[n, ci[1], nchild*(i-1)+1], ci[1]]/(tau+lambda))^2 -
                                  (propose_x[ci[1]] - 0.5*tau*xOld[historyIndex[n, ci[1], nchild*(i-1)+1], ci[1]]/(tau+lambda))^2
        }
        r1 <- r1 + 0.5*(tau+lambda)*sum((x[n, (ci[1]+1):(ci[1]+nv-1)] - (0.5*tau*xOld[left_ancestor_coordinates] + lambda*x[n, ci[1]:(ci[1]+nv-2)])/(tau+lambda))^2 -
                (propose_x[(ci[1]+1):(ci[1]+nv-1)] - (0.5*tau*xOld[left_ancestor_coordinates] + lambda*propose_x[ci[1]:(ci[1]+nv-2)])/(tau+lambda))^2) +
              0.5*(tau+lambda)*sum((x[n, (ci[1]+nv+1):ci[2]] - (0.5*tau*xOld[right_ancestor_coordinates] + lambda*x[n, (ci[1]+nv-1):(ci[2]-1)])/(tau+lambda))^2 -
                           (propose_x[, (ci[1]+nv+1):ci[2]] - (0.5*tau*xOld[right_ancestor_coordinates] + lambda*propose_x[, (ci[1]+nv-1):(ci[2]-1)])/(tau+lambda))^2)
      }
      # ratio for parent node
      r2 <- 0.5*(tau+lambda)*sum((x[n, (ci[1]+1):ci[2]] - (0.5*tau*xOld[crossedover_ancestor_coordinates] + lambda*x[n, ci[1]:(ci[2]-1)])/(tau+lambda))^2 -
                                     (propose_x[(ci[1]+1):ci[2]] - (0.5*tau*xOld[crossedover_ancestor_coordinates] + lambda*propose_x[ci[1]:(ci[2]-1)])/(tau+lambda))^2
      )
      if(ci[1] == 1){
        r2 <- r2 + 0.5*tau*((x[n, 1] - 0.5*xOld[historyIndex[n, 1, nchild*(i-1)+1], 1])^2
                            - (propose_x[1] - 0.5*xOld[historyIndex[n, 1, nchild*(i-1)+1], 1])^2)
      } else {
        r2 <- r2 + 0.5*(tau+lambda)*((x[n, ci[1]] - 0.5*tau*xOld[historyIndex[n, ci[1], nchild*(i-1)+1], ci[1]]/(tau+lambda))^2
                            - (propose_x[ci[1]] - 0.5*tau*xOld[historyIndex[n, ci[1], nchild*(i-1)+1], ci[1]]/(tau+lambda))^2)
      }

      mh_ratio <- r_obs + (1-current_alpha)*r1 + current_alpha*r2

      if(runif(1) <= exp(mh_ratio)){
        x[n, ] <- propose_x
      }
    }

    current_alpha <- new_alpha
    lOmega <- newlOmega
  }
  return(list("x" = x, "history_index_updated" = historyIndexNew))
}



