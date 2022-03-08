get_neighbours_weights <- function(row, col, d){
  neighbours_distance <- c(1, 1, 0, 1, 1)
  out <- neighbours_lattice(row, col, d)
  mixture_weights <- 1/(neighbours_distance+delta) * out$current_valid
  mixture_weights <- mixture_weights/sum(mixture_weights)
  return(list("mixture_weights" = mixture_weights, "current_x_neighbours" = out$current_x_neighbours))
}
# function returning the indices of the neighbours of i on a 2D lattice
neighbours_lattice <- function(row, col, d){
  x_neighbours <- matrix(c(-1, 0, 0, 0, 1, 0, -1, 0, 1, 0), ncol = 2)
  current_x_neighbours <- x_neighbours + matrix(rep(c(row, col), each = 5), ncol = 2)
  current_valid <- as.logical(rowSums((current_x_neighbours > 0) * (current_x_neighbours < d+1)) - 1)
  return(list("current_x_neighbours" = current_x_neighbours, "current_valid" = current_valid))
}
sample_mixture <- function(n, mixture_weights, current_x_neighbours, xOld){
  # sample component of mixture
  mixture_component <- sample.int(5, size = 1, prob = mixture_weights)
  component_coordinates <- current_x_neighbours[mixture_component, ]
  xMean <- xOld[component_coordinates[1], component_coordinates[2], n]
  return(xMean)
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

# function checking if given indices are in current node
in_node <- function(x, ci){
  x[x %in% ci]
}



neighbours_fit_fun <- function(coordinate, mx, obs){
  obs[coordinate] - mx[coordinate]
}
neighbours_fit_product <- function(coordinate, current_node_fit, neighbours_fit){
  current_node_fit[coordinate]*unlist(neighbours_fit[coordinate])
}
