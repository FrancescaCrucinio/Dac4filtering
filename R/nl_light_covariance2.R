nl_light_covariance2 <- function(u, obs, x, history, historyIndex_left, historyIndex_right, cir_left, cir_right, cic_left, cic_right,
                                lW_left, lW_right, sigmaX, m, tau, nu){
  d <- dim(history)[1]
  Nparticles <- dim(history)[3]
  # binary tree
  nchild <- 2
  nodes_dimension <- nchild^u
  nodes_dimension_child <- nchild^(u-1)
  cic <- unique(c(cic_left, cic_right))
  cir <- unique(c(cir_left, cir_right))
  # resample on each children
  if(u == 1){
    lW_left <- c(lW_left)
    max.lW_left <- max(lW_left)
    W_left <- exp(lW_left - max.lW_left)
    lW_right<- c(lW_right)
    max.lW_right <- max(lW_right)
    W_right <- exp(lW_right - max.lW_right)
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

    # latent state contribution + right child observation contribution
    obs_weight_right <- 0
    # left child observation contribution
    obs_weight_left <- 0
    # merged node observation contribution
    obs_weight_merged <- 0
    # merged x
    mx <- x[, , indices1[n]]
    mx[as.matrix(expand.grid(cir_right, cic_right))] <- x[as.matrix(expand.grid(cir_right, cic_right, indices2[n]))]
    for (col in cic) {
      for (row in cir) {
        out_neighbours <- get_neighbours_weights(row, col, d)
        valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
        valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
        main_node <- which(valid_weights == max(valid_weights))
        obs_precision <- rep(tau, length(valid_weights))
        obs_precision[main_node] <- 1
        vcn_node_coordinates <- (valid_current_neighbours[, 1] %in% cir) & (valid_current_neighbours[, 2] %in% cic)
        vcn_left_coordinates <- (valid_current_neighbours[, 1] %in% cir_left) & (valid_current_neighbours[, 2] %in% cic_left)
        vcn_right_coordinates <- (valid_current_neighbours[, 1] %in% cir_right) & (valid_current_neighbours[, 2] %in% cic_right)
        vcn_node <- valid_current_neighbours[vcn_node_coordinates, , drop = FALSE]
        vcn_left <- valid_current_neighbours[vcn_left_coordinates, , drop = FALSE]
        vcn_right <- valid_current_neighbours[vcn_right_coordinates, , drop = FALSE]
        # current node
        current_node_fit <- obs[valid_current_neighbours[main_node, , drop = FALSE]] -
          mx[valid_current_neighbours[main_node, , drop = FALSE]]
        obs_fit <- obs[vcn_node] - mx[vcn_node]
        obs_weight_merged <- obs_weight_merged + current_node_fit*sum(obs_fit*obs_precision[vcn_node_coordinates])
        # left node
        if((valid_current_neighbours[main_node, 2] %in% cic_left) & (valid_current_neighbours[main_node, 1] %in% cir_left)){
          current_node_fit <- obs[valid_current_neighbours[main_node, , drop = FALSE]] -
            x[valid_current_neighbours[main_node, 1], valid_current_neighbours[main_node, 2], indices1[n]]
          obs_fit <- obs[vcn_left] - x[cbind(vcn_left, rep(indices1[n], nrow(vcn_left)))]
          obs_weight_left <- obs_weight_left + current_node_fit*sum(obs_fit*obs_precision[vcn_left_coordinates])
        }
        # right node
        if((valid_current_neighbours[main_node, 2] %in% cic_right) & (valid_current_neighbours[main_node, 1] %in% cir_right)){
          current_node_fit <- obs[valid_current_neighbours[main_node, , drop = FALSE]] -
            x[valid_current_neighbours[main_node, 1], valid_current_neighbours[main_node, 2], indices2[n]]
          obs_fit <- obs[vcn_right] - x[cbind(vcn_right, rep(indices2[n], nrow(vcn_right)))]
          obs_weight_right <- obs_weight_right + current_node_fit*sum(obs_fit*obs_precision[vcn_right_coordinates])
        }
        # contribute of f
        if((col %in% cic_right) & (row %in% cir_right)){ # right child
          lWmix[n] <- lWmix[n] + log(sum(valid_weights * exp(-(x[row, col, indices2[n]] - left_ancestor[valid_current_neighbours])^2/(2*sigmaX)))) -
            log(sum(valid_weights * exp(-(x[row, col, indices2[n]] - right_ancestor[valid_current_neighbours])^2/(2*sigmaX))))
        }
      }
    }
    lWmix[n] <- lWmix[n] - 0.5*(nu+nodes_dimension)*log(1+obs_weight_merged/nu) + 0.5*(nu+nodes_dimension_child)*(log(1+obs_weight_left/nu) + log(1+obs_weight_right/nu))
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}

