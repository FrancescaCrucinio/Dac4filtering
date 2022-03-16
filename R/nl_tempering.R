nl_tempering <- function(ess_target, ess_decay_threshold, x, history, historyIndex_left, historyIndex_right, obs, nv, nu, after_mix_lW, cir, cic, cic_right, cir_right, mcmc_sd){

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
      for (n in 1:Nparticles){
        left_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_left[n, , ]))
        right_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_right[n, , ]))
        left_ancestor <- matrix(history[left_ancestor_coordinates], nrow = d)
        right_ancestor <- matrix(history[right_ancestor_coordinates], nrow = d)

        for (col in cic_right) {
          for (row in cir_right) {
            out_neighbours <- get_neighbours_weights(row, col, d)
            valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
            valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
            after_mix_lW[n] <- after_mix_lW[n] + log(sum(valid_weights * dnorm(x[row, col, n], mean = left_ancestor[valid_current_neighbours], sd = sqrt(sigmaX)))) -
              log(sum(valid_weights * dnorm(x[row, col, n], mean = right_ancestor[valid_current_neighbours], sd = sqrt(sigmaX))))
          }
        }
      }
      if(u == 1){
        after_mix_lW <- after_mix_lW + c(lW_left) + c(lW_right)
      }
    }
    updated_particles <- nl_mcmc_move(x, history, historyIndex_left, historyIndex_right, obs,
                                      nv, nu, after_mix_lW, cir, cic, cic_right, cir_right, mcmc_sd)
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
    after_mix_lW <- rep(0, times = Nparticles)
    for (n in 1:Nparticles){
      left_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_left[n, , ]))
      right_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_right[n, , ]))
      left_ancestor <- matrix(history[left_ancestor_coordinates], nrow = d)
      right_ancestor <- matrix(history[right_ancestor_coordinates], nrow = d)

      for (col in cic_right) {
        for (row in cir_right) {
          out_neighbours <- get_neighbours_weights(row, col, d)
          valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
          valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
          after_mix_lW[n] <- after_mix_lW[n] + log(sum(valid_weights * dnorm(x[row, col, n], mean = left_ancestor[valid_current_neighbours], sd = sqrt(sigmaX)))) -
            log(sum(valid_weights * dnorm(x[row, col, n], mean = right_ancestor[valid_current_neighbours], sd = sqrt(sigmaX))))
        }
      }
    }
    if(u == 1){
      after_mix_lW <- after_mix_lW + c(lW_left) + c(lW_right)
    }
  }
  updated_particles <- nl_mcmc_move(x, history, historyIndex_left, historyIndex_right, obs,
                                    nv, nu, after_mix_lW, cir, cic, cic_right, cir_right, mcmc_sd)
  x <- updated_particles$x
  return(list("x" = x, "history_index_updated" = historyIndexNew, "lWmix" = after_mix_lW))
}

nl_mcmc_move <- function(x, history, historyIndex_left, historyIndex_right, obs, nv, nu, after_mix_lW, cir, cic, cic_right, cir_right, mcmc_sd){
  Nparticles <- dim(x)[3]
  d <- dim(x)[1]
  propose_x <- x
  propose_x[cir, cic, ] <- x[cir, cic, ] +
          mcmc_sd*array(rnorm(Nparticles*length(cir)*length(cic)), dim = c(length(cir), length(cic), Nparticles))
  # mixture weight for proposed particle
  after_mix_lW_new <- rep(0, times = Nparticles)
  for (n in 1:Nparticles){
    left_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_left[n, , ]))
    right_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_right[n, , ]))
    left_ancestor <- matrix(history[left_ancestor_coordinates], nrow = d)
    right_ancestor <- matrix(history[right_ancestor_coordinates], nrow = d)

    for (col in cic_right) {
      for (row in cir_right) {
        out_neighbours <- get_neighbours_weights(row, col, d)
        valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
        valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
        after_mix_lW_new[n] <- after_mix_lW_new[n] + log(sum(valid_weights * dnorm(propose_x[row, col, n], mean = left_ancestor[valid_current_neighbours], sd = sqrt(sigmaX)))) -
          log(sum(valid_weights * dnorm(propose_x[row, col, n], mean = right_ancestor[valid_current_neighbours], sd = sqrt(sigmaX))))
      }
    }
  }
  # mixture weights ratio
  r2 <- after_mix_lW_new - after_mix_lW
  if(nv == 1){# children are leaves
    child_weights_proposal <- -0.5*(nu+1)*log(1+(propose_x[cir_left, cic_left, ] - obs[cir_left, cic_left])^2/nu) -
      0.5*(nu+1)*log(1+(propose_x[cir_right, cic_right, ] - obs[cir_right, cic_right])^2/nu)
    child_weights_current <- -0.5*(nu+1)*log(1+(x[cir_left, cic_left, ] - obs[cir_left, cic_left])^2/nu) -
      0.5*(nu+1)*log(1+(x[cir_right, cic_right, ] - obs[cir_right, cic_right])^2/nu)
    r2 <- r2 + child_weights_proposal - child_weights_current
  }
  r_obs <- rep(0, times = Nparticles)
  r1 <- rep(0, times = Nparticles)
  for (n in 1:Nparticles){
    left_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_left[n, , ]))
    right_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_right[n, , ]))
    left_ancestor <- matrix(history[left_ancestor_coordinates], nrow = d)
    right_ancestor <- matrix(history[right_ancestor_coordinates], nrow = d)
    ft_ratio <- 0
    for (col in cic) {
      for (row in cir) {
        # observation ratio
        robs[n] <- robs[n] - 0.5*(nu+1)*log(1+(propose_x[row, col, n] - obs[row, col])^2/nu) + 0.5*(nu+1)*log(1+(x[row, col, n] - obs[row, col])^2/nu)
        # ft ratio
        out_neighbours <- get_neighbours_weights(row, col, d)
        r1[n] <- r1[n] + log(sum(out_neighbours$mixture_weights * dnorm(propose_x[row, col, n], mean = right_ancestor[row, col], sd = sqrt(sigmaX)))) +
          log(sum(out_neighbours$mixture_weights * dnorm(propose_x[row, col, n], mean = left_ancestor[row, col], sd = sqrt(sigmaX)))) -
          log(sum(out_neighbours$mixture_weights * dnorm(x[row, col, n], mean = right_ancestor[row, col], sd = sqrt(sigmaX)))) -
          log(sum(out_neighbours$mixture_weights * dnorm(x[row, col, n], mean = left_ancestor[row, col], sd = sqrt(sigmaX))))
      }
    }
  }
  mh_ratio <- r1 + r_obs + new_alpha*r2
  accepted <- runif(Nparticles) <= exp(mh_ratio)
  # print(paste(sum(accepted)/Nparticles))
  x[cir, cic, accepted] <- propose_x[cir, cic, accepted]
  after_mix_lW[accepted] <- after_mix_lW_new[accepted]
  return(list("x" = x, "after_mix_lW" = after_mix_lW))
}

