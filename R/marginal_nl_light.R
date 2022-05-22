marginal_nl_light <- function(u_info, x, history, cir_left, cir_right, cic_left, cic_right,
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
  for (n in 1:(theta*Nparticles)){
    sum_left_m <- 0
    sum_right_m <- 0
    sum_merged_m <- 0
    for (m in 1:Nparticles){
      sum_left <- 0
      sum_right <- 0
      sum_merged <- 0
      for (col in cic) {
        for (row in cir) {
          out_neighbours <- get_neighbours_weights(row, col, d)
          valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
          valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
          ancestor_coordinates <- cbind(valid_current_neighbours, rep(m, times = nrow(valid_current_neighbours)))
          # left child
          if((col %in% cic_left) & (row %in% cir_left)){
            sum_left <- sum_left +
              log(sum(valid_weights * exp(-(x[row, col, indices1[n]] - history[ancestor_coordinates])^2/(2*sigmaX*(nlevels-u_info$u+2)))/(sigmaX*sqrt(nlevels-u_info$u+2))))
          }
          # right child
          if((col %in% cic_right) & (row %in% cir_right)){
            sum_right <- sum_right +
              log(sum(valid_weights * exp(-(x[row, col, indices2[n]] - history[ancestor_coordinates])^2/(2*sigmaX*(nlevels-u_info$u+2)))/(sigmaX*sqrt(nlevels-u_info$u+2))))
          }
          # current node
          # merged x
          mx <- x[, , indices1[n]]
          mx[as.matrix(expand.grid(cir_right, cic_right))] <- x[as.matrix(expand.grid(cir_right, cic_right, indices2[n]))]
          sum_merged <- sum_merged +
            log(sum(valid_weights * exp(-(mx[row, col] - history[ancestor_coordinates])^2/(2*sigmaX*(nlevels-u_info$u+1)))/(sigmaX*sqrt(nlevels-u_info$u+1))))
        }
      }
      sum_merged_m <- sum_merged_m + exp(sum_merged)
      sum_left_m <- sum_left_m + exp(sum_left)
      sum_right_m <- sum_right_m + exp(sum_right)
    }
    lWmix[n] <- log(sum_merged_m) - log(sum_left_m) - log(sum_right_m)
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}
