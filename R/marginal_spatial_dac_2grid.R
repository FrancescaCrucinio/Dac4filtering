marginal_dac_spatial_2grid <- function(history, obs, sigmaX, nu){
  # dimension and number of particles
  d <- nrow(history)
  Nparticles <- dim(history)[3]
  ess_target <- Nparticles
  y.error.prec <- matrix(c(1, -0.25, -0.25, 0, -0.25, 1, 0, -0.25, -0.25, 0, 1, -0.25, 0, -0.25, -0.25, 1), nrow = 4)
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # number of variables
  nv <- 1
  # leaves
  sample_from_past <- sample.int(Nparticles, Nparticles, replace = TRUE)
  x <- history[, , sample_from_past] + sqrt(sigmaX)*array(rnorm(Nparticles*d^2), dim = c(d, d, Nparticles))
  lW <- -0.5*(nu+1)*log(1+y.error.prec[1, 1]*sweep(x, 1:2, obs)^2/nu)
  ### tree has only one level
  # number of variables in each node
  nvNew <- 2
  # row indices of children (column 1 = left child, column 2 = right child)
  cir <- matrix(c(1:nv, (nv+1):nvNew), ncol = 2)
  # column indices of children (column 1 = left child, column 2 = right child)
  cic <- matrix(c(1:nv, (nv+1):nvNew), ncol = 2)
  #### HORIZONTAL MERGE ###
  ### Step 1
  lW_left <- lW[1, 1, ]
  lW_right <- lW[1, 2, ]
  # mixture weights
  lWmix <- rep(0, times = Nparticles)
  for (n in 1:Nparticles){
    # merged x
    mx <- x[, , n]
    # contribution of g_{t, u}
    sum_over_neighbours_obs_left <- (obs[1, 1] - mx[1, 1]) %*% y.error.prec[1, 1] %*% (obs[1, 1] - mx[1, 1])
    sum_over_neighbours_obs_right <- (obs[1, 2] - mx[1, 2]) %*% y.error.prec[3, 3] %*% (obs[1, 2] - mx[1, 2])
    sum_over_neighbours_obs_merged <- (obs[1, ] - mx[1, ]) %*% y.error.prec[c(1, 3), c(1, 3)] %*% (obs[1, ] - mx[1, ])
    # contribution of f_{t, u}
    transition_integral <- sweep(history[1, , ], 1, mx[1, ])^2/(2*sigmaX)
    transition_node <- log(sum(exp(-colSums(transition_integral))))
    transition_left <- log(sum(exp(-transition_integral[1, ])))
    transition_right <- log(sum(exp(-transition_integral[2, ])))
    lWmix[n] <- - 0.5*(nu+2)*log(1+abs(sum_over_neighbours_obs_merged)/nu) +
      0.5*(nu+1)*(log(1+abs(sum_over_neighbours_obs_left)/nu) + log(1+abs(sum_over_neighbours_obs_right)/nu)) +
      transition_node - transition_right - transition_left
  }
  lWmix <- lWmix + c(lW_left) + c(lW_right)
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # build ESS
  ess_s <- sum(Wmix)
  ess_ss <- sum(Wmix^2)
  ess <- ess_s^2/ess_ss
  # first permutation
  permutation <- 1:Nparticles
  theta <- 1
  while (ess < ess_target & theta<=ceiling(sqrt(Nparticles))) {
    theta <- theta+1
    new_perm <- sample.int(Nparticles)
    # mixture weights
    lWmix_perm <- rep(0, times = Nparticles)
    for (n in 1:Nparticles){
      # merged x
      mx <- x[, , n]
      mx[1, 2] <- x[1, 2, new_perm[n]]
      # contribution of g_{t, u}
      sum_over_neighbours_obs_left <- (obs[1, 1] - mx[1, 1]) %*% y.error.prec[1, 1] %*% (obs[1, 1] - mx[1, 1])
      sum_over_neighbours_obs_right <- (obs[1, 2] - mx[1, 2]) %*% y.error.prec[3, 3] %*% (obs[1, 2] - mx[1, 2])
      sum_over_neighbours_obs_merged <- (obs[1, ] - mx[1, ]) %*% y.error.prec[c(1, 3), c(1, 3)] %*% (obs[1, ] - mx[1, ])
      # contribution of f_{t, u}
      transition_integral <- sweep(history[1, , ], 1, mx[1, ])^2/(2*sigmaX)
      transition_node <- log(sum(exp(-colSums(transition_integral))))
      transition_left <- log(sum(exp(-transition_integral[1, ])))
      transition_right <- log(sum(exp(-transition_integral[2, ])))
      lWmix_perm[n] <- - 0.5*(nu+2)*log(1+abs(sum_over_neighbours_obs_merged)/nu) +
        0.5*(nu+1)*(log(1+abs(sum_over_neighbours_obs_left)/nu) + log(1+abs(sum_over_neighbours_obs_right)/nu)) +
        transition_node - transition_right - transition_left
    }
    lWmix_perm <- lWmix_perm + c(lW_left) + c(lW_right[new_perm])
    permutation <- c(permutation, new_perm)
    max.lWmix <- max(lWmix_perm)
    Wmix <- exp(lWmix_perm - max.lWmix)
    # build ESS
    ess_s <- ess_s + sum(Wmix)
    ess_ss <- ess_ss + sum(Wmix^2)
    ess <- ess_s^2/ess_ss
    lWmix <- c(lWmix, lWmix_perm)
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  indices <- cbind(rep(1:Nparticles, times = theta)[indices], permutation[indices])
  x[1, 1, ] <- x[1, 1, indices[, 1]]
  x[1, 2, ] <- x[1, 2, indices[, 2]]

  ### Step 2 ###
  lW_left <- lW[2, 1, ]
  lW_right <- lW[2, 2, ]
  # mixture weights
  lWmix <- rep(0, times = Nparticles)
  for (n in 1:Nparticles){
    # merged x
    mx <- x[, , n]
    # contribution of g_{t, u}
    sum_over_neighbours_obs_left <- (obs[2, 1] - mx[2, 1]) %*% y.error.prec[2, 2] %*% (obs[2, 1] - mx[2, 1])
    sum_over_neighbours_obs_right <- (obs[2, 2] - mx[2, 2]) %*% y.error.prec[4, 4] %*% (obs[2, 2] - mx[2, 2])
    sum_over_neighbours_obs_merged <- (obs[2, ] - mx[2, ]) %*% y.error.prec[c(2, 4), c(2, 4)] %*% (obs[2, ] - mx[2, ])
    # contribution of f_{t, u}
    transition_integral <- sweep(history[2, , ], 1, mx[2, ])^2/(2*sigmaX)
    transition_node <- log(sum(exp(-colSums(transition_integral))))
    transition_left <- log(sum(exp(-transition_integral[1, ])))
    transition_right <- log(sum(exp(-transition_integral[2, ])))
    lWmix[n] <- - 0.5*(nu+2)*log(1+abs(sum_over_neighbours_obs_merged)/nu) +
      0.5*(nu+1)*(log(1+abs(sum_over_neighbours_obs_left)/nu) + log(1+abs(sum_over_neighbours_obs_right)/nu)) +
      transition_node - transition_right - transition_left
  }
  lWmix <- lWmix + c(lW_left) + c(lW_right)
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # build ESS
  ess_s <- sum(Wmix)
  ess_ss <- sum(Wmix^2)
  ess <- ess_s^2/ess_ss
  # first permutation
  permutation <- 1:Nparticles
  theta <- 1
  while (ess < ess_target & theta<=ceiling(sqrt(Nparticles))) {
    theta <- theta+1
    new_perm <- sample.int(Nparticles)
    # mixture weights
    lWmix_perm <- rep(0, times = Nparticles)
    for (n in 1:Nparticles){
      # merged x
      mx <- x[, , n]
      mx[2, 2] <- x[2, 2, new_perm[n]]
      # contribution of g_{t, u}
      sum_over_neighbours_obs_left <- (obs[2, 1] - mx[2, 1]) %*% y.error.prec[2, 2] %*% (obs[2, 1] - mx[2, 1])
      sum_over_neighbours_obs_right <- (obs[2, 2] - mx[2, 2]) %*% y.error.prec[4, 4] %*% (obs[2, 2] - mx[2, 2])
      sum_over_neighbours_obs_merged <- (obs[2, ] - mx[2, ]) %*% y.error.prec[c(2, 4), c(2, 4)] %*% (obs[2, ] - mx[2, ])
      # contribution of f_{t, u}
      transition_integral <- sweep(history[2, , ], 1, mx[2, ])^2/(2*sigmaX)
      transition_node <- log(sum(exp(-colSums(transition_integral))))
      transition_left <- log(sum(exp(-transition_integral[1, ])))
      transition_right <- log(sum(exp(-transition_integral[2, ])))
      lWmix_perm[n] <- - 0.5*(nu+2)*log(1+abs(sum_over_neighbours_obs_merged)/nu) +
        0.5*(nu+1)*(log(1+abs(sum_over_neighbours_obs_left)/nu) + log(1+abs(sum_over_neighbours_obs_right)/nu)) +
        transition_node - transition_right - transition_left
    }
    lWmix_perm <- lWmix_perm + c(lW_left) + c(lW_right[new_perm])
    permutation <- c(permutation, new_perm)
    max.lWmix <- max(lWmix_perm)
    Wmix <- exp(lWmix_perm - max.lWmix)
    # build ESS
    ess_s <- ess_s + sum(Wmix)
    ess_ss <- ess_ss + sum(Wmix^2)
    ess <- ess_s^2/ess_ss
    lWmix <- c(lWmix, lWmix_perm)
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  indices <- cbind(rep(1:Nparticles, times = theta)[indices], permutation[indices])
  x[2, 1, ] <- x[2, 1, indices[, 1]]
  x[2, 2, ] <- x[2, 2, indices[, 2]]

  #### VERTICAL MERGE ###
  # mixture weights
  lWmix <- rep(0, times = Nparticles)
  for (n in 1:Nparticles){
    # merged x
    mx <- x[, , n]
    # contribution of g_{t, u}
    sum_over_neighbours_obs_left <- (obs[1, ] - mx[1, ]) %*% y.error.prec[c(1, 3), c(1, 3)] %*% (obs[1, ] - mx[1, ])
    sum_over_neighbours_obs_right <- (obs[2, ] - mx[2, ]) %*% y.error.prec[c(2, 4), c(2, 4)] %*% (obs[2, ] - mx[2, ])
    sum_over_neighbours_obs_merged <- c(obs - mx) %*% y.error.prec %*% c(obs - mx)
    # contribution of f_{t, u}
    transition_integral <- sweep(history, 1:2, mx)^2/(2*sigmaX)
    transition_node <- log(sum(exp(-colSums(transition_integral, dims = 2))))
    transition_left <- log(sum(exp(-colSums(transition_integral[1, , ]))))
    transition_right <- log(sum(exp(-colSums(transition_integral[2, , ]))))
    lWmix[n] <- - 0.5*(nu+4)*log(1+abs(sum_over_neighbours_obs_merged)/nu) +
      0.5*(nu+2)*(log(1+abs(sum_over_neighbours_obs_left)/nu) + log(1+abs(sum_over_neighbours_obs_right)/nu)) +
      transition_node - transition_right - transition_left
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # build ESS
  ess_s <- sum(Wmix)
  ess_ss <- sum(Wmix^2)
  ess <- ess_s^2/ess_ss
  # first permutation
  permutation <- 1:Nparticles
  theta <- 1
  while (ess < ess_target & theta<=ceiling(sqrt(Nparticles))) {
    theta <- theta+1
    new_perm <- sample.int(Nparticles)
    # mixture weights
    lWmix_perm <- rep(0, times = Nparticles)
    for (n in 1:Nparticles){
      # merged x
      mx <- x[, , n]
      mx[2, ] <- x[2, , new_perm[n]]
      # contribution of g_{t, u}
      sum_over_neighbours_obs_left <- (obs[1, ] - mx[1, ]) %*% y.error.prec[c(1, 3), c(1, 3)] %*% (obs[1, ] - mx[1, ])
      sum_over_neighbours_obs_right <- (obs[2, ] - mx[2, ]) %*% y.error.prec[c(2, 4), c(2, 4)] %*% (obs[2, ] - mx[2, ])
      sum_over_neighbours_obs_merged <- c(obs - mx) %*% y.error.prec %*% c(obs - mx)
      # contribution of f_{t, u}
      transition_integral <- sweep(history, 1:2, mx)^2/(2*sigmaX)
      transition_node <- log(sum(exp(-colSums(transition_integral, dims = 2))))
      transition_left <- log(sum(exp(-colSums(transition_integral[1, , ]))))
      transition_right <- log(sum(exp(-colSums(transition_integral[2, , ]))))
      lWmix_perm[n] <- - 0.5*(nu+4)*log(1+abs(sum_over_neighbours_obs_merged)/nu) +
        0.5*(nu+2)*(log(1+abs(sum_over_neighbours_obs_left)/nu) + log(1+abs(sum_over_neighbours_obs_right)/nu)) +
        transition_node - transition_right - transition_left
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
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  indices <- cbind(rep(1:Nparticles, times = theta)[indices], permutation[indices])
  x[1, , ] <- x[1, , indices[, 1]]
  x[2, , ] <- x[2, , indices[, 2]]
  return(x)
}
