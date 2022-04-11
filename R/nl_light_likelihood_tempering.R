nl_light_likelihood_tempering <- function(u, x, obs, history, historyIndex_left, historyIndex_right,
                                          cir_left, cir_right, cic_left, cic_right, lW_left, lW_right, sigmaX, nu, m, beta_diff){
  d <- dim(history)[1]
  Nparticles <- dim(history)[3]
  # binary tree
  nchild <- 2
  # resample on each children
  lW_left <- c(lW_left)
  max.lW_left <- max(lW_left)
  W_left <- exp(lW_left - max.lW_left)
  lW_right<- c(lW_right)
  max.lW_right <- max(lW_right)
  W_right <- exp(lW_right - max.lW_right)
  if(u == 1){
    # child 1
    indices1 <- stratified_resample(W_left/sum(W_left), m*Nparticles)
    # child 2 (with random permutation)
    indices2 <- sample(stratified_resample(W_right/sum(W_right), m*Nparticles))
  }
  else{
    # child 1
    indices1 <- sample.int(Nparticles, size = m*Nparticles, replace = TRUE)
    # child 2
    indices2 <- sample.int(Nparticles, size = m*Nparticles, replace = TRUE)
  }
  # mixture weights
  lWmix <- rep(0, times = m*Nparticles)
  for (n in 1:(m*Nparticles)){
    left_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_left[indices1[n], , ]))
    right_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_right[indices2[n], , ]))
    left_ancestor <- matrix(history[left_ancestor_coordinates], nrow = d)
    right_ancestor <- matrix(history[right_ancestor_coordinates], nrow = d)
    for (col in cic_right) {
      for (row in cir_right) {
        out_neighbours <- get_neighbours_weights(row, col, d)
        valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
        valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
        lWmix[n] <- lWmix[n] + log(sum(valid_weights * dnorm(x[row, col, indices2[n]], mean = left_ancestor[valid_current_neighbours], sd = sqrt(sigmaX)))) -
          log(sum(valid_weights * dnorm(x[row, col, indices2[n]], mean = right_ancestor[valid_current_neighbours], sd = sqrt(sigmaX)))) +
          beta_diff*0.5*(nu+1)*log(1+(x[row, col, ] - obs[row, col])^2/nu)
      }
    }
    for (col in cic_left) {
      for (row in cir_left) {
        lWmix[n] <- lWmix[n] + beta_diff*0.5*(nu+1)*log(1+(x[row, col, ] - obs[row, col])^2/nu)
      }
    }
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}
