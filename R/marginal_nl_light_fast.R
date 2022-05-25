marginal_nl_light_fast <- function(u_info, x, history, cir_left, cir_right, cic_left, cic_right,
                              lW_left, lW_right, sigmaX, theta){
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
  # mixture weights
  lWmix <- rep(0, times = theta*Nparticles)
  all_nodes <- as.matrix(expand.grid(cir, cic))
  how_many_nodes <- nrow(all_nodes)
  nodes_left <- rep((all_nodes[, 1] %in% cir_left) & (all_nodes[, 2] %in% cic_left), each = 5)
  nodes_right <- rep((all_nodes[, 1] %in% cir_right) & (all_nodes[, 2] %in% cic_right), each = 5)
  all_nodes_neighbours <- do.call(rbind, apply(all_nodes, 1, get_all_neighbours, d, simplify = FALSE))
  all_nodes_replicates <- all_nodes[rep(1:how_many_nodes, each = 5), ]
  all_nodes_neighbours[all_nodes_neighbours[, 3] == 0, 1:2] <- 1


  # how_many_nodes <- length(all_nodes_lexicographic)
  # all_nodes_neighbours <- matrix(c(all_nodes_lexicographic - d, all_nodes_lexicographic - 1,
  #                           all_nodes_lexicographic, all_nodes_lexicographic + 1,
  #                           all_nodes_lexicographic + d), nrow = how_many_nodes)
  # all_nodes_neighbours_distances <- rep(1/(1+delta), times = 5*how_many_nodes)
  # all_nodes_neighbours_distances[(2*how_many_nodes+1):(3*how_many_nodes)] <- 1/delta
  # all_nodes_neighbours_valid <- (all_nodes_neighbours > 0) & (all_nodes_neighbours <= d^2)
  # mixture_weights <- matrix(all_nodes_neighbours_distances*all_nodes_neighbours_valid, nrow = how_many_nodes, byrow = FALSE)
  # mixture_weights <- mixture_weights/rowSums(mixture_weights)
  # mixture_weights <- c(t(mixture_weights))
  # all_nodes_valid <- rep(all_nodes_lexicographic, each = 5)
  # all_nodes_neighbours_valid_coordinates <- c(t(all_nodes_neighbours_valid))*c(t(all_nodes_neighbours))
  # all_nodes_neighbours_valid_coordinates[all_nodes_neighbours_valid_coordinates == 0] <- 1
  for (n in 1:(theta*Nparticles)){
    sum_merged <- 0
    sum_left <- 0
    sum_right <- 0
    # merged x
    mx <- x[, , indices1[n]]
    mx[as.matrix(expand.grid(cir_right, cic_right))] <- x[as.matrix(expand.grid(cir_right, cic_right, indices2[n]))]
    # ancestors_indices <- matrix(all_nodes_neighbours_valid, nrow = Nparticles, ncol = length(all_nodes_neighbours_valid), byrow = TRUE) + d^2*(0:(Nparticles-1))
    # history_all_particles <- matrix(history[ancestors_indices], ncol = length(all_nodes_neighbours_valid))
    # tmp <-  sweep(exp(-(sweep(history_all_particles, 2, mx[all_nodes_valid]))^2/(2*sigmaX)), 2, mixture_weights, "*")
    # tmp_left <- sweep(tmp, 2, rep(nodes_left, each = 5), "*")
    # tmp_right <- sweep(tmp, 2, rep(nodes_right, each = 5), "*")
    # sum_over_neighbours_all <- log(apply(tmp, 1, sum_over_nodes, how_many_nodes))
    # sum_over_neighbours_left <- log(matrix(apply(tmp_left, 1, sum_over_nodes, how_many_nodes), ncol = Nparticles))
    # sum_over_neighbours_right <- log(matrix(apply(tmp_right, 1, sum_over_nodes, how_many_nodes), ncol = Nparticles))
    # lWmix[n] <- log(sum(exp(colSums(sum_over_neighbours_all)))) - log(sum(exp(colSums(sum_over_neighbours_left)))) -
    #   log(sum(exp(colSums(sum_over_neighbours_right))))
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
    lWmix[n] <- log(sum_merged) - log(sum_left) - log(sum_right)
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}

marginal_nl_light_fast_adaptive <- function(ess_target, u_info, x, history, cir_left, cir_right, cic_left, cic_right,
                                   lW_left, lW_right, sigmaX){
  d <- dim(history)[1]
  Nparticles <- dim(history)[3]
  nlevels <- log2(d)
  cic <- unique(c(cic_left, cic_right))
  cir <- unique(c(cir_left, cir_right))
  all_nodes <- as.matrix(expand.grid(cir, cic))
  how_many_nodes <- nrow(all_nodes)
  nodes_left <- rep((all_nodes[, 1] %in% cir_left) & (all_nodes[, 2] %in% cic_left), each = 5)
  nodes_right <- rep((all_nodes[, 1] %in% cir_right) & (all_nodes[, 2] %in% cic_right), each = 5)
  all_nodes_neighbours <- do.call(rbind, apply(all_nodes, 1, get_all_neighbours, d, simplify = FALSE))
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
    lWmix[n] <- log(sum_merged) - log(sum_left) - log(sum_right)
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
      lWmix_perm[n] <- log(sum_merged) - log(sum_left) - log(sum_right)
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
