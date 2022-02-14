lgssm_tempering <- function(ci, nv, lambda, tau, sigmaY, x, xOld, historyIndex, after_mix_lW, ess_target, ess_decay_threshold, mcmc_sd){
  Nparticles <- length(after_mix_lW)
  current_alpha <- bisection_ess(after_mix_lW, ess_target)
  lOmega <- rep(0, times = Nparticles)
  resampled_indices <- 1:Nparticles
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
    xNew <- matrix(0, nrow = Nparticles, ncol = d)
    if(ess < Nparticles / 2){
      resampled_indices <- stratified_resample(omega/sum(omega), Nparticles)
      newlOmega <- rep(0, times = Nparticles)
      xNew[, ci[1]:ci[2]] <- cbind(x[resampled_indices, ci[1]:(ci[1]+nv-1)], x[resampled_indices, (ci[1]+nv):ci[2]])
      # update lW
      after_mix_lW <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOldv)
      )
    }
    # MCMC move
    for (n in 1:Nparticles) {
      propose_x <- xNew[n, ]
      propose_x[ci[1]:ci[2]] <- xNew[n, ci[1]:ci[2]] + mcmc_sd*rnorm(length(ci[1]:ci[2]))
      right_ancestor_coordinates <- cbind(historyIndex[n, (ci[1]+nv):ci[2], nchild*i], (ci[1]+nv):ci[2])
      left_ancestor_coordinates <- cbind(historyIndex[n, (ci[1]+1):ci[2], nchild*(i-1)+1], (ci[1]+1):ci[2])
      r_obs <- sum((xNew[n, ci[1]:ci[2]] - obs[ci[1]:ci[2]])^2-(propose_x[ci[1]:ci[2]] - obs[ci[1]:ci[2]])^2)/(2*sigmaY)

      # r1 <- 0.5*(tau+lambda)*sum((xNew[n, (ci[1]+1):(ci[1]+nv-1)] - (0.5*tau*xOld[, (ci[1]+1):(ci[1]+nv-1)] + lambda*xNew[, ci[1]:(ci[1]+nv-2)])/(tau+lambda))^2
                                     # - (propose_x[, (ci[1]+1):(ci[1]+nv-1)] - (0.5*tau*xOld[, (ci[1]+1):(ci[1]+nv-1)] + lambda*propose_x[, ci[1]:(ci[1]+nv-2)])/(tau+lambda))^2) +
      #   0.5*(tau+lambda)*rowSums((xNew[, (ci[1]+nv):ci[2]] - (0.5*tau*xOld[, (ci[1]+nv):ci[2]] + lambda*xNew[, (ci[1]+nv):(ci[2]-1)])/(tau+lambda))^2
      #                            - (propose_x[, (ci[1]+nv+1):ci[2]] - (0.5*tau*xOld[, (ci[1]+nv+1):ci[2]] + lambda*propose_x[, (ci[1]+nv):(ci[2]-1)])/(tau+lambda))^2
      #   )
      # r2 <- 0.5*(tau+lambda)*sum((xNew[n, (ci[1]+1):ci[2]] - (0.5*tau*xOld[left_ancestor_coordinates] + lambda*xNew[n, ci[1]:(ci[2]-1)])/(tau+lambda))^2 -
      #                                (propose_x[(ci[1]+1):ci[2]] - (0.5*tau*xOld[left_ancestor_coordinates] + lambda*propose_x[ci[1]:(ci[2]-1)])/(tau+lambda))^2
      # )
      # if(ci[1] == 1){
      #   # r2 <- r2 + 0.5*tau*((xNew[, ci[1]] - 0.5*xOld[, ci[1]])^2 - (propose_x[, ci[1]] - 0.5*xOld[, ci[1]])^2)
      #   r2 <- r2 + 0.5*tau*((xNew[n, 1] - 0.5*xOld[historyIndex[n, 1, nchild*(i-1)+1], 1])^2 - (propose_x[1] - 0.5*xOld[historyIndex[n, 1, nchild*(i-1)+1], 1])^2)
      # } else{
      #   r2 <- r2 + 0.5*tau*((xNew[n, ci[1]] - (0.5*tau*xOld[historyIndex[n, ci[1], nchild*(i-1)+1], ci[1]] + lambda*xNew[n, ci[1]+1])/(tau+lambda))^2
      #                       - (propose_x[ci[1]] - - (0.5*tau*xOld[historyIndex[n, ci[1], nchild*(i-1)+1], ci[1]] + lambda*propose_x[ci[1]+1])/(tau+lambda))^2)
      # }


      mh_ratio <- r_obs + (1-current_alpha)*r1 + current_alpha*r2
      if(runif(1) <= exp(mh_ratio)){
        xNew[n, ] <- propose_x
      }
    }

    current_alpha <- new_alpha
    lOmega <- newlOmega
    x <- xNew
  }
  return(list("x" = x, "resampled_indices" = resampled_indices))
}

