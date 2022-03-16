nl_light_covariance <- function(u, x, history, historyIndex_left, historyIndex_right, cir_left, cir_right, cic_left, cic_right,
                     lW_left, lW_right, sigmaX, Nparticles, m, d, nodes_dimension){
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
  for (n in 1:m*Nparticles){
    left_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_left[indices1[n], , ]))
    right_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_right[indices2[n], , ]))
    left_ancestor <- matrix(history[left_ancestor_coordinates], nrow = d)
    right_ancestor <- matrix(history[right_ancestor_coordinates], nrow = d)

    # latent state contribution + right child observation contribution
    obs_weight_right <- 0
    for (col in cic_right) {
      current_node_fit <- obs[col, col] - x[col, col, indices2[n]]
      for (row in cir_right) {
        out_neighbours <- get_neighbours_weights(row, col, d)
        valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
        valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
        lWmix[n] <- lWmix[n] + log(sum(valid_weights * dnorm(x[row, col, indices2[n]], mean = left_ancestor[valid_current_neighbours], sd = sqrt(sigmaX)))) -
          log(sum(valid_weights * dnorm(x[row, col, indices2[n]], mean = right_ancestor[valid_current_neighbours], sd = sqrt(sigmaX))))
        node_distance <- abs(row - col)
        if(node_distance <= 1 & col %in% cir_right){
          obs_fit <- obs[row, col] - x[row, col, indices2[n]]
          obs_weight_right <- obs_weight_right + current_node_fit*obs_fit/4^(node_distance)
        }
      }
    }
    # left child observation contribution
    obs_weight_left <- 0
    for (col in cic_left) {
      current_node_fit <- obs[col, col] - x[col, col, indices1[n]]
      for (row in cir_left) {
        node_distance <- abs(row - col)
        if(node_distance <= 1 & col %in% cir_left){
          obs_fit <- obs[row, col] - x[row, col, indices1[n]]
          obs_weight_left <- obs_weight_left + current_node_fit*obs_fit/4^(node_distance)
        }
      }
    }
    # merged node observation contribution
    obs_weight_merged <- 0
    # merged x
    mx <- x[, , indices1[n]]
    mx[as.matrix(expand.grid(cir_right, cic_right))] <- x[as.matrix(expand.grid(cir_right, cic_right, indices2[n]))]
    for (col in unique(c(cic_left, cic_right))) {
      current_node_fit <- obs[col, col] - mx[col, col]
      for (row in unique(c(cir_left, cir_right))) {
        node_distance <- abs(row - col)
        if(node_distance <= 1 & col %in% c(cir_left, cir_right)){
          obs_fit <- obs[row, col] - mx[row, col]
          obs_weight_merged <- obs_weight_merged + current_node_fit*obs_fit/4^(node_distance)
        }
      }
    }
    lWmix[n] <- lWmix[n] + 0.5*(nu+nodes_dimension)*(-log(1+obs_weight_merged/nu) + log(1+obs_weight_left/nu) + log(1+obs_weight_right/nu))
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}


nl_adaptive_light_covariance <- function(ess_target, u, x, history, historyIndex_left, historyIndex_right, cir_left, cir_right, cic_left, cic_right,
                                         lW_left, lW_right, sigmaX, Nparticles, m, d, nodes_dimension){
  # binary tree
  nchild <- 2
  # mixture weights
  lWmix <- rep(0, times = Nparticles)
  for (n in 1:Nparticles){
    left_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_left[n, , ]))
    right_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_right[n, , ]))
    left_ancestor <- matrix(history[left_ancestor_coordinates], nrow = d)
    right_ancestor <- matrix(history[right_ancestor_coordinates], nrow = d)

    # latent state contribution + right child observation contribution
    obs_weight_right <- 0
    for (col in cic_right) {
      current_node_fit <- obs[col, col] - x[col, col, n]
      for (row in cir_right) {
        out_neighbours <- get_neighbours_weights(row, col, d)
        valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
        valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
        lWmix[n] <- lWmix[n] + log(sum(valid_weights * dnorm(x[row, col, n], mean = left_ancestor[valid_current_neighbours], sd = sqrt(sigmaX)))) -
          log(sum(valid_weights * dnorm(x[row, col, n], mean = right_ancestor[valid_current_neighbours], sd = sqrt(sigmaX))))
        node_distance <- abs(row - col)
        if(node_distance <= 1 & col %in% cir_right){
          obs_fit <- obs[row, col] - x[row, col, n]
          obs_weight_right <- obs_weight_right + current_node_fit*obs_fit/4^(node_distance)
        }
      }
    }
    # left child observation contribution
    obs_weight_left <- 0
    for (col in cic_left) {
      current_node_fit <- obs[col, col] - x[col, col, n]
      for (row in cir_left) {
        node_distance <- abs(row - col)
        if(node_distance <= 1 & col %in% cir_left){
          obs_fit <- obs[row, col] - x[row, col, n]
          obs_weight_left <- obs_weight_left + current_node_fit*obs_fit/4^(node_distance)
        }
      }
    }
    # merged node observation contribution
    obs_weight_merged <- 0
    # merged x
    mx <- x[, , n]
    mx[as.matrix(expand.grid(cir_right, cic_right))] <- x[as.matrix(expand.grid(cir_right, cic_right, n))]
    for (col in unique(c(cic_left, cic_right))) {
      current_node_fit <- obs[col, col] - mx[col, col]
      for (row in unique(c(cir_left, cir_right))) {
        node_distance <- abs(row - col)
        if(node_distance <= 1 & col %in% c(cir_left, cir_right)){
          obs_fit <- obs[row, col] - mx[row, col]
          obs_weight_merged <- obs_weight_merged + current_node_fit*obs_fit/4^(node_distance)
        }
      }
    }
    lWmix[n] <- lWmix[n] + 0.5*(nu+nodes_dimension)*(log(1+obs_weight_merged/nu) - log(1+obs_weight_left/nu) - log(1+obs_weight_right/nu))
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
      for (col in cic_right) {
        current_node_fit <- obs[col, col] - x[col, col, new_perm[n]]
        for (row in cir_right) {
          out_neighbours <- get_neighbours_weights(row, col, d)
          valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
          valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
          lWmix_perm[n] <- lWmix_perm[n] + log(sum(valid_weights * dnorm(x[row, col, n], mean = left_ancestor[valid_current_neighbours], sd = sqrt(sigmaX)))) -
            log(sum(valid_weights * dnorm(x[row, col, new_perm[n]], mean = right_ancestor[valid_current_neighbours], sd = sqrt(sigmaX))))
          node_distance <- abs(row - col)
          if(node_distance <= 1 & col %in% cir_right){
            obs_fit <- obs[row, col] - x[row, col, new_perm[n]]
            obs_weight_right <- obs_weight_right + current_node_fit*obs_fit/4^(node_distance)
          }
        }
      }
      # left child observation contribution
      obs_weight_left <- 0
      for (col in cic_left) {
        current_node_fit <- obs[col, col] - x[col, col, n]
        for (row in cir_left) {
          node_distance <- abs(row - col)
          if(node_distance <= 1 & col %in% cir_left){
            obs_fit <- obs[row, col] - x[row, col, n]
            obs_weight_left <- obs_weight_left + current_node_fit*obs_fit/4^(node_distance)
          }
        }
      }
      # merged node observation contribution
      obs_weight_merged <- 0
      # merged x
      mx <- x[, , n]
      mx[as.matrix(expand.grid(cir_right, cic_right))] <- x[as.matrix(expand.grid(cir_right, cic_right, new_perm[n]))]
      for (col in c(cic_left, cic_right)) {
        current_node_fit <- obs[col, col] - mx[col, col]
        for (row in c(cir_left, cic_right)) {
          node_distance <- abs(row - col)
          print(paste(row, col))
          if(node_distance <= 1 & col %in% c(cir_left, cir_right)){
            obs_fit <- obs[row, col] - mx[row, col]
            obs_weight_merged <- obs_weight_merged + current_node_fit*obs_fit/4^(node_distance)
            print(paste(obs_weight_merged))
          }
        }
      }
      lWmix_perm[n] <- lWmix_perm[n] + 0.5*(nu+nodes_dimension)*(log(1+obs_weight_merged/nu) - log(1+obs_weight_left/nu) - log(1+obs_weight_right/nu))

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
  }
  write.table(data.frame("u" = u, "m" = m), file = "data/adaptive_nl.csv", sep = ",", append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(rep(1:Nparticles, times = m)[indices], permutation[indices]),
              "target_reached" = (ess >= ess_target), "resampled_particles_lW" = lWmix[indices]))
}

