# Generate observations from a nonlinear model
nl_obs <- function(d, sigmaX, nu, delta, y.error.var, Time.step){
  x <- array(0, dim = c(d, d, Time.step+1))
  y <- array(0, dim = c(d, d, Time.step))
  y_iid <- array(0, dim = c(d, d, Time.step))
  # initial state
  x[, , 1] <- sqrt(sigmaX)*matrix(rnorm(d^2), nrow = d)

  for (t in 2:(Time.step+1)) {
    for (row in 1:d) {
      for(col in 1:d){
        out_neighbours <- get_neighbours_weights(row, col, d)
        mixture_component <- sample.int(5, size = 1, prob = out_neighbours$mixture_weights)
        component_coordinates <- out_neighbours$current_x_neighbours[mixture_component, ]
        x[row, col, t] <- x[component_coordinates[1], component_coordinates[2], t-1] + sqrt(sigmaX)*rnorm(1)
        y_iid[row, col, t-1] <- x[row, col, t] + rt(1, df = nu)
      }
    }
    y[, , t-1] <- matrix(rmvt(1, mu = c(x[, , t]), S = y.error.var, df = 5), ncol = d)

  }
  return(list("x" = x, "y" = y, "yiid" = y_iid))
}
