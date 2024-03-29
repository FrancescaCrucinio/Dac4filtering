marginal_spatial_light <- function(u_info, x, obs, cir_left, cir_right, cic_left, cic_right, lW_left, lW_right, sigmaX, tau, tau_diag, nu, history){
  Nparticles <- dim(history)[3]
  theta <- ceiling(sqrt(Nparticles))
  d <- dim(x)[1]
  Nparticles <- dim(x)[3]
  cic <- unique(c(cic_left, cic_right))
  cir <- unique(c(cir_left, cir_right))
  nchild <- 2
  if(u_info$direction == "h"){
    u <- 2*u_info$u-1
    nodes_dimension <- nchild^u
    nodes_dimension_child <- nchild^(u-1)
  } else {
    u <- 2*u_info$u
    nodes_dimension <- nchild^u
    nodes_dimension_child <- nchild^(u-1)
  }
  all_nodes <- as.matrix(expand.grid(cir, cic))
  all_nodes_left <- cbind(all_nodes[(all_nodes[, 1] %in% cir_left) & (all_nodes[, 2] %in% cic_left), 1],
                          all_nodes[(all_nodes[, 1] %in% cir_left) & (all_nodes[, 2] %in% cic_left), 2])
  all_nodes_right <- cbind(all_nodes[(all_nodes[, 1] %in% cir_right) & (all_nodes[, 2] %in% cic_right), 1],
                           all_nodes[(all_nodes[, 1] %in% cir_right) & (all_nodes[, 2] %in% cic_right), 2])
  how_many_nodes <- nrow(all_nodes)
  how_many_nodes_left <- nrow(all_nodes_left)
  how_many_nodes_right <- nrow(all_nodes_right)
  all_nodes_history <- cbind(do.call(rbind, replicate(Nparticles, all_nodes, simplify=FALSE)), rep(1:Nparticles, each = how_many_nodes))
  all_nodes_history_left <- cbind(do.call(rbind, replicate(Nparticles, all_nodes_left, simplify=FALSE)), rep(1:Nparticles, each =how_many_nodes_left))
  all_nodes_history_right <- cbind(do.call(rbind, replicate(Nparticles, all_nodes_right, simplify=FALSE)), rep(1:Nparticles, each = how_many_nodes_right))
  all_nodes_neighbours <- do.call(rbind, apply(all_nodes, 1, get_all_neighbours, d, tau, tau_diag, simplify = FALSE))
  neighbours_in_node <- (all_nodes_neighbours[, 1] %in% cir) & (all_nodes_neighbours[, 2] %in% cic)
  neighbours_in_left <- ((all_nodes_neighbours[, 1] %in% cir_left) & (all_nodes_neighbours[, 2] %in% cic_left)) *
    rep((all_nodes[, 1] %in% cir_left) & (all_nodes[, 2] %in% cic_left), each = 5)
  neighbours_in_right <- ((all_nodes_neighbours[, 1] %in% cir_right) & (all_nodes_neighbours[, 2] %in% cic_right)) *
    rep((all_nodes[, 1] %in% cir_right) & (all_nodes[, 2] %in% cic_right), each = 5)
  all_nodes_replicates <- all_nodes[rep(1:how_many_nodes, each = 5), ]
  all_nodes_neighbours[all_nodes_neighbours[, 3] == 0, 1:2] <- 1
  # resample on each children
  if(u == 1){
    lW_left <- c(lW_left)
    max.lW_left <- max(lW_left)
    W_left <- exp(lW_left - max.lW_left)
    lW_right<- c(lW_right)
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
  for (n in 1:(theta*Nparticles)){
    # merged x
    mx <- x[, , indices1[n]]
    mx[as.matrix(expand.grid(cir_right, cic_right))] <- x[as.matrix(expand.grid(cir_right, cic_right, indices2[n]))]
    # contribution of g_{t, u}
    tmp_obs <- neighbours_in_node*all_nodes_neighbours[, 4]*(obs[all_nodes_replicates] - mx[all_nodes_replicates])*(obs[all_nodes_neighbours[, 1:2]] - mx[all_nodes_neighbours[, 1:2]])
    tmp_obs_left <- tmp_obs * neighbours_in_left
    tmp_obs_right <- tmp_obs * neighbours_in_right
    sum_over_neighbours_obs_left <- sum(tmp_obs_left)
    sum_over_neighbours_obs_right <-sum(tmp_obs_right)
    sum_over_neighbours_obs_merged <- sum(tmp_obs)
    # contribution of f_{t, u}
    transition_node <- mean(exp(-rowSums(sweep(matrix(history[all_nodes_history], ncol = how_many_nodes, byrow = TRUE), 2, mx[all_nodes])^2/(2*sigmaX))))
    transition_left <- mean(exp(-rowSums(sweep(matrix(history[all_nodes_history_left], ncol = how_many_nodes_left, byrow = TRUE), 2, mx[all_nodes_left])^2/(2*sigmaX))))
    transition_right <- mean(exp(-rowSums(sweep(matrix(history[all_nodes_history_right], ncol = how_many_nodes_right, byrow = TRUE), 2, mx[all_nodes_right])^2/(2*sigmaX))))
    # transition_node <- ifelse(all(transition_node < .Machine$double.eps), 0, log(transition_node))
    # transition_left[n] <- ifelse(all(transition_left< .Machine$double.eps), 0, log(transition_left))
    # transition_right <- ifelse(all(transition_right< .Machine$double.eps), 0, log(transition_right))
    lWmix[n] <- - 0.5*(nu+nodes_dimension)*log(1+abs(sum_over_neighbours_obs_merged)/nu) +
      0.5*(nu+nodes_dimension_child)*(log(1+abs(sum_over_neighbours_obs_left)/nu) + log(1+abs(sum_over_neighbours_obs_right)/nu)) +
      log(transition_node) - log(transition_right) - log(transition_left)
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}
