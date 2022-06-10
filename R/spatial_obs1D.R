# Generate observations from a spatial model
spatial_obs1D <- function(sigmaX, nu, tau_diag, Time.step){
  y.error.prec <- tau_diag
  y.error.var <- 1/y.error.prec
  x <- rep(0, times = Time.step+1)
  y <- rep(0, times = Time.step)
  # initial state
  x[1] <- sqrt(sigmaX)*rnorm(1)

  for (t in 2:(Time.step+1)) {
    x[t] <- x[t-1] + sqrt(sigmaX)*rnorm(1)
    y[t-1] <- x[t] + rt.scaled(1, df = nu, mean = 0, sd = tau_diag)
  }
  return(list("x" = x, "y" = y))
}
