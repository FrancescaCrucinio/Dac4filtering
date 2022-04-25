nl_light_covariance <- function(u, obs, x, history, historyIndex_left, historyIndex_right, cir_left, cir_right, cic_left, cic_right,
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
      current_node_fit <- obs[col, col] - mx[col, col]
      for (row in cir) {
        node_distance <- abs(row - col)
        if(node_distance <= 1 & col %in% cir){
          obs_fit <- obs[row, col] - mx[row, col]
          obs_weight_merged <- obs_weight_merged + current_node_fit*obs_fit*tau^(node_distance)
        }
        if((col %in% cic_left) & (row %in% cir_left)){ # left child
          current_node_fit <- obs[col, col] - x[col, col, indices1[n]]
          node_distance <- abs(row - col)
          if(node_distance <= 1 & col %in% cir_left){
            obs_fit <- obs[row, col] - x[row, col, indices1[n]]
            obs_weight_left <- obs_weight_left + current_node_fit*obs_fit*tau^(node_distance)
          }
        }
        if((col %in% cic_right) & (row %in% cir_right)){ # right child
          current_node_fit <- obs[col, col] - x[col, col, indices2[n]]
          out_neighbours <- get_neighbours_weights(row, col, d)
          valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
          valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
          lWmix[n] <- lWmix[n] + log(sum(valid_weights * exp(-(x[row, col, indices2[n]] - left_ancestor[valid_current_neighbours])^2/(2*sigmaX)))) -
            log(sum(valid_weights * exp(-(x[row, col, indices2[n]] - right_ancestor[valid_current_neighbours])^2/(2*sigmaX))))
          node_distance <- abs(row - col)
          if(node_distance <= 1 & col %in% cir_right){
            obs_fit <- obs[row, col] - x[row, col, indices2[n]]
            obs_weight_right <- obs_weight_right + current_node_fit*obs_fit*tau^(node_distance)
          }
        }
      }
    }
    lWmix[n] <- lWmix[n] - 0.5*(nu+nodes_dimension)*log(1+obs_weight_merged/nu)
          + 0.5*(nu+nodes_dimension_child)*(log(1+obs_weight_left/nu) + log(1+obs_weight_right/nu))
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}


nl_adaptive_light_covariance <- function(ess_target, u, obs, x, history, historyIndex_left, historyIndex_right, cir_left, cir_right, cic_left, cic_right,
                                         lW_left, lW_right, sigmaX, tau, nu, u_info){
  d <- dim(history)[1]
  Nparticles <- dim(history)[3]
  # binary tree
  nchild <- 2
  nodes_dimension <- nchild^u
  nodes_dimension_child <- nchild^(u-1)
  cic <- unique(c(cic_left, cic_right))
  cir <- unique(c(cir_left, cir_right))
  # mixture weights
  lWmix <- rep(0, times = Nparticles)
  for (n in 1:Nparticles){
    left_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_left[n, , ]))
    right_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_right[n, , ]))
    left_ancestor <- matrix(history[left_ancestor_coordinates], nrow = d)
    right_ancestor <- matrix(history[right_ancestor_coordinates], nrow = d)

    # latent state contribution + right child observation contribution
    obs_weight_right <- 0
    # left child observation contribution
    obs_weight_left <- 0
    # merged node observation contribution
    obs_weight_merged <- 0
    # merged x
    mx <- x[, , n]
    mx[as.matrix(expand.grid(cir_right, cic_right))] <- x[as.matrix(expand.grid(cir_right, cic_right, n))]
    for (col in cic) {
      current_node_fit <- obs[col, col] - mx[col, col]
      for (row in cir) {
        node_distance <- abs(row - col)
        if(node_distance <= 1 & col %in% cir){
          obs_fit <- obs[row, col] - mx[row, col]
          obs_weight_merged <- obs_weight_merged + current_node_fit*obs_fit*tau^(node_distance)
        }
        if((col %in% cic_left) & (row %in% cir_left)){ # left child
          current_node_fit <- obs[col, col] - x[col, col, n]
          node_distance <- abs(row - col)
          if(node_distance <= 1 & col %in% cir_left){
            obs_fit <- obs[row, col] - x[row, col, n]
            obs_weight_left <- obs_weight_left + current_node_fit*obs_fit*tau^(node_distance)
          }
        }
        if((col %in% cic_right) & (row %in% cir_right)){ # right child
          current_node_fit <- obs[col, col] - x[col, col, n]
          out_neighbours <- get_neighbours_weights(row, col, d)
          valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
          valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
          lWmix[n] <- lWmix[n] + log(sum(valid_weights * exp(-(x[row, col, n] - left_ancestor[valid_current_neighbours])^2/(2*sigmaX)))) -
            log(sum(valid_weights * exp(-(x[row, col, n] - right_ancestor[valid_current_neighbours])^2/(2*sigmaX))))
          node_distance <- abs(row - col)
          if(node_distance <= 1 & col %in% cir_right){
            obs_fit <- obs[row, col] - x[row, col, n]
            obs_weight_right <- obs_weight_right + current_node_fit*obs_fit*tau^(node_distance)
          }
        }
      }
    }
    lWmix[n] <- lWmix[n] - 0.5*(nu+nodes_dimension)*log(1+obs_weight_merged/nu)
        + 0.5*(nu+nodes_dimension_child)*(log(1+obs_weight_left/nu) + log(1+obs_weight_right/nu))
  }
  if(u == 1){
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
  m <- 1
  while (ess < ess_target & m<=ceiling(sqrt(Nparticles))) {
    m <- m+1
    new_perm <- sample.int(Nparticles)
    # mixture weights
    lWmix_perm <- rep(0, times = Nparticles)
    for (n in 1:Nparticles) {
      left_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_left[n, , ]))
      right_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_right[new_perm[n], , ]))
      left_ancestor <- matrix(history[left_ancestor_coordinates], nrow = d)
      right_ancestor <- matrix(history[right_ancestor_coordinates], nrow = d)

      # latent state contribution + right child observation contribution
      obs_weight_right <- 0
      # left child observation contribution
      obs_weight_left <- 0
      # merged node observation contribution
      obs_weight_merged <- 0
      # merged x
      mx <- x[, , n]
      mx[as.matrix(expand.grid(cir_right, cic_right))] <- x[as.matrix(expand.grid(cir_right, cic_right, new_perm[n]))]
      print(paste(cic, cir))
      for (col in cic) {

        current_node_fit <- obs[col, col] - mx[col, col]
        for (row in cir) {
          print(paste(col, row))
          node_distance <- abs(row - col)
          if(node_distance <= 1 & col %in% cir){
            print(paste("squared", current_node_fit))
            obs_fit <- obs[row, col] - mx[row, col]
            print(paste("product", obs_fit))
            obs_weight_merged <- obs_weight_merged + current_node_fit*obs_fit*tau^(node_distance)
          }
          if((col %in% cic_left) & (row %in% cir_left)){ # left child
            current_node_fit <- obs[col, col] - x[col, col, n]
            node_distance <- abs(row - col)
            if(node_distance <= 1 & col %in% cir_left){
              obs_fit <- obs[row, col] - x[row, col, n]
              obs_weight_left <- obs_weight_left + current_node_fit*obs_fit*tau^(node_distance)
            }
          }
          if((col %in% cic_right) & (row %in% cir_right)){ # right child
            current_node_fit <- obs[col, col] - x[col, col, new_perm[n]]
            out_neighbours <- get_neighbours_weights(row, col, d)
            valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
            valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
            print(paste(lWmix_perm[n]))
            lWmix_perm[n] <- lWmix_perm[n] + log(sum(valid_weights * exp(-(x[row, col, new_perm[n]] - left_ancestor[valid_current_neighbours])^2/(2*sigmaX)))) -
              log(sum(valid_weights * exp(-(x[row, col, new_perm[n]] - right_ancestor[valid_current_neighbours])^2/(2*sigmaX))))
            node_distance <- abs(row - col)
            if(node_distance <= 1 & col %in% cir_right){
              obs_fit <- obs[row, col] - x[row, col, new_perm[n]]
              obs_weight_right <- obs_weight_right + current_node_fit*obs_fit*tau^(node_distance)
            }
          }
        }
      }
      print(paste("weight", lWmix_perm[n]))
      lWmix_perm[n] <- lWmix_perm[n] - 0.5*(nu+nodes_dimension)*log(1+obs_weight_merged/nu)
          + 0.5*(nu+nodes_dimension_child)*(log(1+obs_weight_left/nu) + log(1+obs_weight_right/nu))
      print(paste("weight", lWmix_perm[n]))
      print(paste("merge", obs_weight_merged))
      if(is.nan(lWmix_perm[n])){ break}
    }

    if(u == 1){
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

    print(paste(ess))
  }
  write.table(data.frame("u" = u, "m" = m), file = "data/adaptive_nl_cov.csv", sep = ",", append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(rep(1:Nparticles, times = m)[indices], permutation[indices]),
              "target_reached" = (ess >= ess_target), "resampled_particles_lW" = lWmix[indices]))
}

