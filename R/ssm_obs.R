# Generate observations from a linear gaussian SSM
ssm_obs <- function(mu0, Sigma0, y.coeff, x.coeff, x.error.prec, y.error.var, Time.step){
  p <- nrow(y.error.var)
  y <- matrix(0, nrow = Time.step, ncol = p)
  # initial state
  x <- mvrnorm(n = 1, mu0, Sigma0)
  # loop over time
  for(i in 1:Time.step){
    x <- t(rmvnp(1, c(x.coeff%*%x), x.error.prec))
    y[i, ] <- mvrnorm(n = 1, y.coeff%*%x, y.error.var)
  }
  return(y)
}
