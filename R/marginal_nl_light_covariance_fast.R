marginal_nl_light_covariance_fast <- function(u_info, x, history, cir_left, cir_right, cic_left, cic_right,
                                   lW_left, lW_right, sigmaX, theta, tau, nu){
  d <- dim(history)[1]
  Nparticles <- dim(history)[3]
  nlevels <- log2(d)
  cic <- unique(c(cic_left, cic_right))
  cir <- unique(c(cir_left, cir_right))
  # resample on each children
  if(u_info$u == 1 & u_info$direction == "h"){
    max.lW_left <- max(lW_left)
    W_left <- exp(lW_left - max.lW_left)
    max.lW_right <- max(lW_right)
    W_right <- exp(lW_right - max.lW_right)
    # child 1
    indices1 <- stratified_resample(W_left/sum(W_left), theta*Nparticles)
    # child 2 (with random permutation)
    indices2 <- sample(stratified_resample(W_right/sum(W_right), theta*Nparticles))
  }
  else{
    # child 1
    indices1 <- sample.int(Nparticles, size = theta*Nparticles, replace = TRUE)
    # child 2
    indices2 <- sample.int(Nparticles, size = theta*Nparticles, replace = TRUE)
  }
  ##### get nodes dimenssion
  # mixture weights
  lWmix <- rep(0, times = theta*Nparticles)
  all_nodes <- as.matrix(expand.grid(cir, cic))
  how_many_nodes <- nrow(all_nodes)
  nodes_left <- rep((all_nodes[, 1] %in% cir_left) & (all_nodes[, 2] %in% cic_left), each = 5)
  nodes_right <- rep((all_nodes[, 1] %in% cir_right) & (all_nodes[, 2] %in% cic_right), each = 5)
  all_nodes_neighbours <- do.call(rbind, apply(all_nodes, 1, get_all_neighbours, d, tau, simplify = FALSE))
  all_nodes_replicates <- all_nodes[rep(1:how_many_nodes, each = 5), ]
  all_nodes_neighbours[all_nodes_neighbours[, 3] == 0, 1:2] <- 1
  for (n in 1:(theta*Nparticles)){
    sum_merged <- 0
    sum_left <- 0
    sum_right <- 0
    # merged x
    mx <- x[, , indices1[n]]
    mx[as.matrix(expand.grid(cir_right, cic_right))] <- x[as.matrix(expand.grid(cir_right, cic_right, indices2[n]))]
    # contribution of f_{t, u}
    for (m in 1:Nparticles){
      history_m <- history[, , m]
      # current node
      tmp <- all_nodes_neighbours[, 3] * exp(-(mx[all_nodes_replicates] - history_m[all_nodes_neighbours[, 1:2]])^2/(2*sigmaX))
      tmp_left <- tmp * nodes_left
      tmp_right <- tmp * nodes_right
      sum_over_neighbours_left <- rowSums(matrix(tmp_left, nrow=how_many_nodes, byrow = TRUE))
      sum_over_neighbours_right <- rowSums(matrix(tmp_right, nrow=how_many_nodes, byrow = TRUE))
      sum_over_neighbours_merged <- rowSums(matrix(tmp, nrow=how_many_nodes, byrow = TRUE))
      sum_left <- sum_left + exp(sum(log(sum_over_neighbours_left[sum_over_neighbours_left>0])))
      sum_right <- sum_right + exp(sum(log(sum_over_neighbours_right[sum_over_neighbours_right>0])))
      sum_merged <- sum_merged + exp(sum(log(sum_over_neighbours_merged[sum_over_neighbours_merged>0])))
    }
    # contribution of g_{t, u}
    tmp_obs <- all_nodes_neighbours[, 4]*(mx[all_nodes_replicates] - obs[all_nodes_replicates])*(mx[all_nodes_neighbours[, 1:2]] - obs[all_nodes_neighbours[, 1:2]])
    tmp_obs_left <- tmp_obs * nodes_left
    tmp_obs_right <- tmp_obs * nodes_right
    sum_over_neighbours_obs_left <- sum(tmp_obs_left)
    sum_over_neighbours_obs_right <-sum(tmp_obs_right)
    sum_over_neighbours_obs_merged <- sum(tmp_obs)
    lWmix[n] <- log(sum_merged) - log(sum_left) - log(sum_right) - 0.5*(nu+nodes_dimension)*log(1+abs(sum_over_neighbours_obs_merged)/nu)
      + 0.5*(nu+nodes_dimension_child)*(log(1+abs(sum_over_neighbours_obs_left)/nu) + log(1+abs(sum_over_neighbours_obs_right)/nu))
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}
