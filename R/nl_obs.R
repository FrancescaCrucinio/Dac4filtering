nl_obs <- function(d, sigmaX, nu, delta, Time.Step){
  y.error.prec <- matrix(0, nrow = d^2, ncol = d^2)
  diag(y.error.prec) <- 1
  diag(y.error.prec[-1, ]) <- 1/4
  vertical_neighbours <- ((0:d^2) * (d^2 + 1) + d+1)
  y.error.prec[vertical_neighbours[(vertical_neighbours <= d^4)]] <- 1/4
  y.error.prec[upper.tri(y.error.prec)] = t(y.error.prec)[upper.tri(y.error.prec)]
  x <- array(0, dim = c(d, d, Time.step+1))
  y <- array(0, dim = c(d, d, Time.step))
  # initial state
  x[, , 1] <- sqrt(sigmaX)*matrix(rnorm(d^2), nrow = 2)

  x_neighbours <- matrix(c(-1, 0, 0, 0, 1, 0, -1, 0, 1, 0), ncol = 2)
  neighbours_distance <- c(1, 1, 0, 1, 1)
  for (t in 2:(Time.step+1)) {
    for (row in 1:d) {
      for(col in 1:d){
        print(paste(row, col))
          current_x_neighbours <- x_neighbours + matrix(rep(c(row, col), each = 5), ncol = 2)
          current_valid <- as.logical(rowSums((current_x_neighbours > 0) * (current_x_neighbours < d+1)) - 1)
          mixture_weights <- 1/(neighbours_distance+delta) * current_valid
          mixture_weights <- mixture_weights/sum(mixture_weights)
          mixture_component <- sample.int(5, size = 1, prob = mixture_weights)
          component_coordinates <- current_x_neighbours[mixture_component, ]
          x[row, col, t] <- x[component_coordinates[1], component_coordinates[2], t-1] + sqrt(sigmaX)*rnorm(1)
      }
    }
    y[, , t-1] <- matrix(rmvt(1, mu = c(x[, , t]), df = nu, sigma = inv(y.error.prec)), ncol = d)

  }
  return(x, y)
}
