# Index of first and last child of node i in level u
child_indices <- function(i, nv){
  rbind(nv*(i-1)+1, i*nv)
}

# Get row and column indices for mixture resampling
rc_indices <- function(indices, dims){
  r <- indices -  trunc(indices/dims)*dims
  r[r == 0] <- dims
  c <- pmin(trunc((indices -1 )/dims)+1, dims)
  return(cbind(r, c))
}

# Kolmogorov-Smirnov distance
ks_dist <- function(x, N){
  ks.test(x[1:N], x[N:length(x)])$statistic
}

# Wasserstein-1 distance
w1_dist <- function(x, N){
  wasserstein1d(x[1:N], x[N:length(x)])
}

# Stratified resampling
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

# Get neighbours weights for the nonlinear spatial model
get_all_neighbours <- function(rowcol, d, tau = NULL, tau_diag = NULL){
  delta <- 1
  neighbours_distance <- c(1, 1, 0, 1, 1)
  x_neighbours <- matrix(c(-1, 0, 0, 0, 1, 0, -1, 0, 1, 0), ncol = 2)
  current_x_neighbours <- x_neighbours + matrix(rep(c(rowcol[1], rowcol[2]), each = 5), ncol = 2)
  current_valid <- as.logical(rowSums((current_x_neighbours > 0) * (current_x_neighbours < d+1)) - 1)
  mixture_weights <- 1/(neighbours_distance+delta) * current_valid
  mixture_weights <- mixture_weights/sum(mixture_weights)
  out <- cbind(current_x_neighbours, mixture_weights)
  if(!is.null(tau)) {
    obs_weights <- tau * neighbours_distance + tau_diag *(1-neighbours_distance)
    obs_weights <- obs_weights * current_valid
    out <- cbind(out, obs_weights)
  }
  return(out)
}
