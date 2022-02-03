nsmc_time_lgssm <- function(tau, lambda, sigmaY, Nparticles, x0, y, M = NULL, marginals = NULL){
  # dimension
  d <- ncol(y)
  # time interval
  Time.step <- nrow(y)
  m <- matrix(0, nrow = Time.step, ncol = d)
  v <- matrix(0, nrow = Time.step, ncol = d)
  # initial value
  x <- x0

  if(is.null(M)) {
    # number of particles per island
    M <- 2*d
  }

  for (t in 1:Time.step) {
    x <- nsmc_lgssm(x, y[t, ], tau, lambda, sigmaY, M)
    m[t, ] <- colMeans(x)
    v[t, ] <- colVars(x)
  }
  if(!is.null(marginals)){
    # compare marginals at last time step
    ks_nsmc <- apply(rbind(x, marginals), ks_dist, N = Nparticles, MARGIN = 2)
    w1_nsmc <- apply(rbind(x, marginals), ks_dist, N = Nparticles, MARGIN = 2)
  }
  out <- list("m" = m, "v" = v, "ks" = ks_nsmc, "w1" = w1_nsmc)
  return(out)
}

nsmc_time_car <- function(sigmaX, sigmaY, Nparticles, x0, y, M = NULL, marginals = NULL){
  # dimension
  d <- ncol(y)
  # time interval
  Time.step <- nrow(y)
  m <- matrix(0, nrow = Time.step, ncol = d)
  v <- matrix(0, nrow = Time.step, ncol = d)
  # initial value
  x <- x0

  if(is.null(M)) {
    # number of particles per island
    M <- d
  }

  for (t in 1:Time.step) {
    x <- nsmc_car(x, y[t, ], sigmaX, sigmaY, M)
    m[t, ] <- colMeans(x)
    v[t, ] <- colVars(x)
  }
  if(!is.null(marginals)){
    # compare marginals at last time step
    ks_nsmc <- apply(rbind(x, marginals), ks_dist, N = Nparticles, MARGIN = 2)
    w1_nsmc <- apply(rbind(x, marginals), ks_dist, N = Nparticles, MARGIN = 2)
  }
  out <- list("m" = m, "v" = v, "ks" = ks_nsmc, "w1" = w1_nsmc, "x" = x)
  return(out)
}
