nl_light <- function(u, x, history, historyIndex_left, historyIndex_right, cir_right, cic_right,
                     lW_left, lW_right, sigmaX, m){
  d <- dim(history)[1]
  Nparticles <- dim(history)[3]
  # binary tree
  nchild <- 2
  # resample on each children
  if(u == 1){
    max.lW_left <- max(lW_left)
    W_left <- exp(lW_left - max.lW_left)
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
    tmp <- as.matrix(expand.grid(cir_right, cic_right))
    neighbours_out <- sapply(1:nrow(tmp), wrap_get_neighbours_weights, tmp, d, simplify = FALSE)
    res_tmp <- do.call(rbind, sapply(neighbours_out, `[[`, 1, simplify = FALSE))
    tmp <- tmp[rep(seq_len(nrow(tmp)), times = sapply(neighbours_out, `[[`, 2)), ]
    lWmix[n] <- log(sum(res_tmp[, 3] * exp(-(x[cbind(tmp, indices2[n])] - left_ancestor[res_tmp[, 1:2]])^2/(2*sigmaX)))) -
            log(sum(res_tmp[, 3] * exp(-(x[cbind(tmp, indices2[n])] - right_ancestor[res_tmp[, 1:2]])^2/(2*sigmaX))))
    # for (col in cic_right) {
    #   for (row in cir_right) {
    #     out_neighbours <- get_neighbours_weights(row, col, d)
    #     valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
    #     valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
    #     lWmix[n] <- lWmix[n] + log(sum(valid_weights * exp(-(x[row, col, indices2[n]] - left_ancestor[valid_current_neighbours])^2/(2*sigmaX)))) -
    #       log(sum(valid_weights * exp(-(x[row, col, indices2[n]] - right_ancestor[valid_current_neighbours])^2/(2*sigmaX))))
    #   }
    # }
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}

nl_adaptive_light <- function(ess_target, u, x, history, historyIndex_left, historyIndex_right, cir_right, cic_right,
                              lW_left, lW_right, sigmaX, u_info){
  d <- dim(history)[1]
  Nparticles <- dim(history)[3]
  # binary tree
  nchild <- 2
  # mixture weights
  lWmix <- rep(0, times = Nparticles)
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
        lWmix[n] <- lWmix[n] + log(sum(valid_weights * exp(-(x[row, col, n] - left_ancestor[valid_current_neighbours])^2/(2*sigmaX)))) -
          log(sum(valid_weights * exp(-(x[row, col, n] - right_ancestor[valid_current_neighbours])^2/(2*sigmaX))))
      }
    }
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

      for (col in cic_right) {
        for (row in cir_right) {
          out_neighbours <- get_neighbours_weights(row, col, d)
          valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
          valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
          lWmix_perm[n] <- lWmix_perm[n] + log(sum(valid_weights * exp(-(x[row, col, n] - left_ancestor[valid_current_neighbours])^2/(2*sigmaX)))) -
            log(sum(valid_weights * exp(-(x[row, col, new_perm[n]] - right_ancestor[valid_current_neighbours])^2/(2*sigmaX))))
        }
      }
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
  write.table(data.frame("u" = u_info$u, "direction" = u_info$direction, "m" = m), file = "data/adaptive_nl.csv", sep = ",", append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(rep(1:Nparticles, times = m)[indices], permutation[indices]),
              "target_reached" = (ess >= ess_target), "resampled_particles_lW" = lWmix[indices]))
}

wrap_get_neighbours_weights <- function(n, mat, d){
  out_neighbours <- get_neighbours_weights(mat[n, 1], mat[n, 2], d)
  valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
  valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
  return(list("neighbours" = matrix(c(valid_current_neighbours, valid_weights), ncol = 3),
              "how_many_neighbours" = nrow(valid_current_neighbours)))
}

