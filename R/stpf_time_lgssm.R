stpf_time_lgssm <- function(tau, lambda, sigmaY, Nparticles, x0, y, M = NULL){
  # dimension
  d <- ncol(y)
  # time interval
  Time.step <- nrow(y)
  lZ <- rep(0, times = Time.step)
  # m <- matrix(0, nrow = Time.step, ncol = d)
  # v <- matrix(0, nrow = Time.step, ncol = d)
  # initial value
  x <- x0

  if(is.null(M)) {
    # number of particles per island
    M <- d
  }

  for (t in 1:Time.step) {
    res_stpf <- stpf_lgssm(x, y[t, ], tau, lambda, sigmaY)
    # x <- res_stpf[1][[1]]
    lZ[t] <- res_stpf[2][[1]]
    # m[t, ] <- colMeans(x, dims = 2)
    # v[t, ] <- colMeans(x^2, dims = 2) - m[t, ]^2
  }
  return(res_stpf)
}
