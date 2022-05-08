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

# Sample frm mixture for nonlinear model
sample_mixture <- function(n, mixture_weights, current_x_neighbours, xOld){
  # sample component of mixture
  mixture_component <- sample.int(5, size = 1, prob = mixture_weights)
  component_coordinates <- current_x_neighbours[mixture_component, ]
  xMean <- xOld[component_coordinates[1], component_coordinates[2], n]
  return(xMean)
}
