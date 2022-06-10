# Get neighbours weights for the nonlinear model
get_neighbours_weights <- function(row, col, d){
  delta <- 1
  neighbours_distance <- c(1, 1, 0, 1, 1)
  out <- neighbours_lattice(row, col, d)
  mixture_weights <- 1/(neighbours_distance+delta) * out$current_valid
  mixture_weights <- mixture_weights/sum(mixture_weights)
  return(list("mixture_weights" = mixture_weights, "current_x_neighbours" = out$current_x_neighbours))
}

# Get indices of the neighbours of i on a 2D lattice
neighbours_lattice <- function(row, col, d){
  x_neighbours <- matrix(c(-1, 0, 0, 0, 1, 0, -1, 0, 1, 0), ncol = 2)
  current_x_neighbours <- x_neighbours + matrix(rep(c(row, col), each = 5), ncol = 2)
  current_valid <- as.logical(rowSums((current_x_neighbours > 0) * (current_x_neighbours < d+1)) - 1)
  return(list("current_x_neighbours" = current_x_neighbours, "current_valid" = current_valid))
}

# Sample from mixture for nonlinear model
sample_mixture <- function(n, mixture_weights, current_x_neighbours, xOld){
  # sample component of mixture
  mixture_component <- sample.int(5, size = 1, prob = mixture_weights)
  component_coordinates <- current_x_neighbours[mixture_component, ]
  xMean <- xOld[component_coordinates[1], component_coordinates[2], n]
  return(xMean)
}

# Sample from mixture for nonlinear model
marginal_sample_mixture <- function(n, mixture_weights, current_x_neighbours, xOld, sigmaX){
  Nparticles <- dim(xOld)[3]
  # sample component of mixture
  mixture_component <- sample.int(5, size = 1, prob = mixture_weights)
  component_coordinates <- current_x_neighbours[mixture_component, ]
  # sample history component
  xMean <- xOld[component_coordinates[1], component_coordinates[2], sample.int(Nparticles, 1)]
  return(xMean)
}

# Get neighbours weights for the nonlinear model
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
    # obs_weights <- tau^(neighbours_distance)* current_valid
    obs_weights <- tau * neighbours_distance + tau_diag *(1-neighbours_distance)
    obs_weights <- obs_weights * current_valid
    out <- cbind(out, obs_weights)
  }
  # current_x_neighbours_lexicographic <- current_x_neighbours[, 1] + (current_x_neighbours[, 2] - 1)*d
  # current_x_neighbours_covariance <- full_covariance[sort(current_x_neighbours_lexicographic), sort(current_x_neighbours_lexicographic)]
  # current_x_neighbours_precision <- solve(current_x_neighbours_covariance)
  return(out)
}

# Get indices of the neighbours of i on a 2D lattice
get_previous_neighbours <- function(row, col, d){
  x_neighbours <- matrix(c(-1, 0, 0, -1), ncol = 2)
  current_x_neighbours <- x_neighbours + matrix(rep(c(row, col), each = 2), ncol = 2)
  current_valid <- as.logical(rowSums((current_x_neighbours > 0) * (current_x_neighbours < d+1)) - 1)
  return(list("current_x_neighbours" = current_x_neighbours, "current_valid" = current_valid))
}
