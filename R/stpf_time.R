stpf_time_lgssm <- function(tau, lambda, sigmaY, Nparticles, x0, y, M = NULL, marginals = NULL){
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
    x <- stpf_lgssm(x, y[t, ], tau, lambda, sigmaY)
    m[t, ] <- colMeans(x, dims = 2)
    v[t, ] <- colMeans(x^2, dims = 2) - m[t, ]^2
  }
  if(!is.null(marginals)){
    # compare marginals at last time step
    ks_stpf <- apply(rbind(matrix(x, ncol = d, nrow = Nparticles*M), marginals), ks_dist, N = M*Nparticles, MARGIN = 2)
    w1_stpf <- apply(rbind(matrix(x, ncol = d, nrow = Nparticles*M), marginals), w1_dist, N = M*Nparticles, MARGIN = 2)
  }
  out <- list("m" = m, "v" = v, "ks" = ks_stpf, "w1" = w1_stpf)
  return(out)
}


stpf_time_car <- function(sigmaX, sigmaY, Nparticles, x0, y, M = NULL, marginals = NULL){
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
    x <- stpf_car(x, y[t, ], sigmaX, sigmaY)
    m[t, ] <- colMeans(x, dims = 2)
    v[t, ] <- colMeans(x^2, dims = 2) - m[t, ]^2
  }
  if(!is.null(marginals)){
    # compare marginals at last time step
    ks_dac <- apply(rbind(matrix(x, ncol = d, nrow = Nparticles*M), marginals), ks_dist, N = Nparticles, MARGIN = 2)
    w1_dac <- apply(rbind(matrix(x, ncol = d, nrow = Nparticles*M), marginals), w1_dist, N = Nparticles, MARGIN = 2)
  }
  out <- list("m" = m, "v" = v, "ks" = ks_dac, "w1" = w1_dac)
  return(out)
}
