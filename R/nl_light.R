nl_light <- function(i, u, nv, ci, W, Nparticles, m, eta, x){
  # binary tree
  nchild <- 2
  valid_neighbours <- sapply(c(ci), neighbours_lattice, d = d)
  u_valid_neighbours <- lapply(valid_neighbours, in_node, ci = c(ci))
  left_valid_neighbours <- lapply(valid_neighbours[1:nv], in_node, ci = ci[, 1])
  right_valid_neighbours <- lapply(valid_neighbours[(nv+1):nvNew], in_node, ci = ci[, 2])
  # resample on each children
  if(u == 1){
    # child 1
    indices1 <- stratified_resample(W[, nchild*(i-1)+1], m*Nparticles)
    # child 2 (with random permutation)
    indices2 <- sample(stratified_resample(W[, nchild*i], m*Nparticles))
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
    # merged x
    mx <- x[indices1[n], ]
    mx[ci[, 2]] <-x[indices2[n], ci[, 2]]
    current_node_fit <- obs[c(ci)] - mx
    neighbours_fit <- lapply(u_valid_neighbours, neighbours_fit_fun, mx, obs)
    left_neighbours_fit <- lapply(left_valid_neighbours, neighbours_fit_fun, x[indices1[n], ], obs)
    right_neighbours_fit <- lapply(right_valid_neighbours, neighbours_fit_fun, x[indices2[n], ], obs)
    product_neighbours_fit <- unlist(sapply(1:nvNew, neighbours_fit_product, current_node_fit, neighbours_fit))
    left_product_neighbours_fit <- unlist(sapply(1:nv, neighbours_fit_product, current_node_fit, left_neighbours_fit))
    right_product_neighbours_fit <- unlist(sapply((nv+1):nvNew, neighbours_fit_product, current_node_fit, right_neighbours_fit))
    lWmix[n] <- 0.5*(nu+d)*log(1 + sum(sum(current_node_fit^2, 0.25*product_neighbours_fit))/nu) -
      0.5*(nu+d)*log(1 + sum(sum(current_node_fit[ci[, 1]]^2, 0.25*left_product_neighbours_fit))/nu) -
      0.5*(nu+d)*log(1 + sum(sum(current_node_fit[ci[, 2]]^2, 0.25*right_product_neighbours_fit))/nu)
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}


neighbours_fit_fun <- function(coordinate, mx, obs){
  obs[coordinate] - mx[coordinate]
}
neighbours_fit_product <- function(coordinate, current_node_fit, neighbours_fit){
  current_node_fit[coordinate]*unlist(neighbours_fit[coordinate])
}
