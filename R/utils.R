bisection_ess <- function(lW, ess_target){
  W <- exp(lW - max(lW))
  ess_alpha <- function(alpha){ sum(W^(alpha))^2/sum(W^(2*alpha)) - ess_target}
  alpha_star <- bisect(ess_alpha, 0, 1)$root
  return(alpha_star)
}
bisection_cess <- function(lOmega, lW, current_alpha, ess_decay_threshold){
  W <- exp(lW - max(lW))
  omega <- exp(lOmega - max(lOmega))
  omega <- omega/sum(omega)
  cess_alpha <- function(alpha){ sum(omega*W^(alpha-current_alpha))^2/sum(omega*W^(2*(alpha-current_alpha))) - ess_decay_threshold}
  if(sign(cess_alpha(current_alpha)) == sign(cess_alpha(1))){
    return(list("alpha_star" = current_alpha, "same_sign" = TRUE))
  } else{
    alpha_star <- bisect(cess_alpha, current_alpha, 1)$root
    return(list("alpha_star" = alpha_star, "same_sign" = FALSE))
  }
}
# function returning index of first and last child of node i in level u
child_indices <- function(i, nv){
  rbind(nv*(i-1)+1, i*nv)
}
# function returning index of first and last child of node i in level u for 2D lattice
# (each column corresponds to one child)
child_indices_lattice <- function(d, u, i, nv, nodes){
  if((u %% 2) == 0){ # connect vertically
    matrix(c((nv*(i-1)+1):(i*nv), (nv*(i-1)+1):(i*nv)+d), ncol = 2)
  } else { # connect horizontally
    if(u == 1){
      matrix((2*(i-1)+1):(i*2), ncol = 2)
    } else {
      nvOld <- nchild^(u-2)
      matrix(c((nvOld*(i-1)+1):(i*nvOld), (nvOld*(i-1)+1):(i*nvOld)+d,
               (nvOld*i+1):((i+1)*nvOld), (nvOld*i+1):((i+1)*nvOld)+d), ncol = 2)
    }
  }
}
# function returning the indices of the neighbours of i on a 2D lattice
neighbours_lattice <- function(d, i){
  # check for boundary
  if((i %% d) == 1){
    x_neighbours <- c(i-d, i+1, i+d)
  } else if((i %% d) == 0){
    x_neighbours <- c(i-d, i-1, i+d)
  } else {
    x_neighbours <- c(i-d, i-1, i+1, i+d)
  }
  valid <- (x_neighbours <= d^2) & (x_neighbours > 0)
  return(x_neighbours[valid])
}
# function checking if given indices are in current node
in_node <- function(x, ci){
  x[x %in% ci]
}
# get row and column indices for mixture resampling
rc_indices <- function(indices, dims){
  r <- indices -  trunc(indices/dims)*dims
  r[r == 0] <- dims
  c <- pmin(trunc((indices -1 )/dims)+1, dims)
  return(cbind(r, c))
}
#' Kolmogorov-Smirnov distance
#'
#' @param x a vector of samples
#' @param N mid point of vector x
#' @return The Kolmogorov-Smirnov distance between x[1:N] and x[N:end]
ks_dist <- function(x, N){
  ks.test(x[1:N], x[N:length(x)])$statistic
}
#' Wasserstein-1 distance
#'
#' @param x a vector of samples
#' @param N mid point of vector x
#' @return The Wasserstein-1 distance between x[1:N] and x[N:end]
w1_dist <- function(x, N){
  wasserstein1d(x[1:N], x[N:length(x)])
}
### RESAMPLING
mult_resample <- function(W, N){
  W <- c(W)
  # vector to store number of offsprings
  indices <- rep(0, times = N)

  # inverse transform
  un <- sort(runif(N))
  s <- W[1]
  m <- 1
  for(i in 1:N){
    while(s < un[i]){
      m <- m+1
      s <- s+W[m]
    }
    indices[i] <- m
  }
  return(indices)
}

#' Stratified resampling
#'
#' @param W vector of weights
#' @param N number of samples to output
#' @return A vector of N indices corresponding to equally weighted samples from the weighted distribution W
stratified_resample <- function(W, N){
  W <- c(W)
  # vector to store number of offsprings
  indices <- rep(0, times = N)

  # inverse transform
  # partition (0, 1] into N disjoint intervals
  # and draw one uniform for every interval
  un <- (0:(N-1)+runif(N))/N
  s <- W[1]
  m <- 1
  for(i in 1:N){
    while(s < un[i]){
      m <- m+1
      s <- s+W[m]
    }
    indices[i] <- m
  }
  return(indices)
}
