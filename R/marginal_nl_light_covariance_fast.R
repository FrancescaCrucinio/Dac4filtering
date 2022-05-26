marginal_nl_light_covariance_fast <- function(u_info, x, history, cir_left, cir_right, cic_left, cic_right,
                                   lW_left, lW_right, sigmaX, theta, tau, nu){
  d <- dim(history)[1]
  Nparticles <- dim(history)[3]
  nlevels <- log2(d)
  cic <- unique(c(cic_left, cic_right))
  cir <- unique(c(cir_left, cir_right))
  if(u_info$direction == "h"){
    u <- 2*u_info$u-1
    nodes_dimension <- nchild^u
    nodes_dimension_child <- nchild^(u-1)
  } else {
    u <- 2*u_info$u
    nodes_dimension <- nchild^u
    nodes_dimension_child <- nchild^(u-1)
  }
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


marginal_nl_light_covariance_fast_adaptive <- function(ess_target, u_info, x, history, cir_left, cir_right, cic_left, cic_right,
                                            lW_left, lW_right, sigmaX, tau, nu){
  d <- dim(history)[1]
  Nparticles <- dim(history)[3]
  nlevels <- log2(d)
  cic <- unique(c(cic_left, cic_right))
  cir <- unique(c(cir_left, cir_right))
  all_nodes <- as.matrix(expand.grid(cir, cic))
  how_many_nodes <- nrow(all_nodes)
  nodes_left <- rep((all_nodes[, 1] %in% cir_left) & (all_nodes[, 2] %in% cic_left), each = 5)
  nodes_right <- rep((all_nodes[, 1] %in% cir_right) & (all_nodes[, 2] %in% cic_right), each = 5)
  all_nodes_neighbours <- do.call(rbind, apply(all_nodes, 1, get_all_neighbours, d, tau, simplify = FALSE))
  all_nodes_replicates <- all_nodes[rep(1:how_many_nodes, each = 5), ]
  all_nodes_neighbours[all_nodes_neighbours[, 3] == 0, 1:2] <- 1
  # mixture weights
  lWmix <- rep(0, times = Nparticles)
  for (n in 1:Nparticles){
    sum_merged <- 0
    sum_left <- 0
    sum_right <- 0
    # merged x
    mx <- x[, , n]
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
  if(u_info$u == 1 & u_info$direction == "h"){
    lWmix <- lWmix + c(lW_left) + c(lW_right)
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
      sum_merged <- 0
      sum_left <- 0
      sum_right <- 0
      # merged x
      mx <- x[, , n]
      mx[as.matrix(expand.grid(cir_right, cic_right))] <- x[as.matrix(expand.grid(cir_right, cic_right, new_perm[n]))]
      for (m in 1:Nparticles){
        history_m <- history[, , m]
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
      lWmix_perm[n] <- log(sum_merged) - log(sum_left) - log(sum_right) - 0.5*(nu+nodes_dimension)*log(1+abs(sum_over_neighbours_obs_merged)/nu)
      + 0.5*(nu+nodes_dimension_child)*(log(1+abs(sum_over_neighbours_obs_left)/nu) + log(1+abs(sum_over_neighbours_obs_right)/nu))
    }
    if(u_info$u == 1 & u_info$direction == "h"){
      lWmix_perm <- lWmix_perm + c(lW_left) + c(lW_right)
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
  write.table(data.frame("u" = u_info$u, "direction" = u_info$direction, "theta" = theta), file = "data/adaptive_nl.csv", sep = ",", append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(rep(1:Nparticles, times = theta)[indices], permutation[indices]),
              "target_reached" = (ess >= ess_target), "resampled_particles_lW" = lWmix[indices]))
}

