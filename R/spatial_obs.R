# Generate observations from a spatial model
spatial_obs <- function(d, sigmaX, nu, tau, tau_diag, Time.step){
  y.error.prec <- matrix(0, nrow = d^2, ncol = d^2)
  diag(y.error.prec) <- tau_diag
  diag(y.error.prec[-1, ]) <- tau
  vertical_neighbours <- ((0:d^2) * (d^2 + 1) + d+1)
  y.error.prec[vertical_neighbours[(vertical_neighbours <= d^4)]] <- tau
  y.error.prec[upper.tri(y.error.prec)] <- t(y.error.prec)[upper.tri(y.error.prec)]
  y.error.var <- inv(y.error.prec)
  y.error.var[upper.tri(y.error.var)] <- t(y.error.var)[upper.tri(y.error.var)]
  x <- array(0, dim = c(d, d, Time.step+1))
  y <- array(0, dim = c(d, d, Time.step))
  # initial state
  x[, , 1] <- sqrt(sigmaX)*matrix(rnorm(d^2), nrow = d)

  for (t in 2:(Time.step+1)) {
    x[, , t] <- x[, , t-1] + sqrt(sigmaX)*matrix(rnorm(d^2), nrow = d, ncol = d)
    tmp = rmvt(1, mu = c(x[, , t]), S = y.error.var, df = nu)
    y[, , t-1] <- matrix(rmvt(1, mu = c(x[, , t]), S = y.error.var, df = nu), ncol = d)
  }

  return(list("x" = x, "y" = y, "variance" = y.error.var*nu/(nu-2)))
}
