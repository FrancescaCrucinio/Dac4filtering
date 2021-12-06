dac_time_lgssm <- function(tau, lambda, sigmaY, Nparticles, x0, y, method = "light", M = NULL, marginals = NULL){
  # dimension
  d <- ncol(y)
  # time interval
  Time.step <- nrow(y)
  m <- matrix(0, nrow = Time.step, ncol = d)
  v <- matrix(0, nrow = Time.step, ncol = d)
  # initial value
  x <- x0

  for (t in 1:Time.step) {
    switch(method,
           lc={
             # linear cost dac
             res_dac <- dac_lgssm_lc(x, y[t, ], tau, lambda, sigmaY, Sigma.det)
           },
           mix={
             # dac with mixture weights
             res_dac <- dac_lgssm(x, y[t, ], tau, lambda, sigmaY, Sigma.det)
           }, # dac with lightweight mixture weights (no adaptation)
           light={
             res_dac <- dac_lgssm_lightweight(x, y[t, ], tau, lambda, sigmaY, Sigma.det)
           }, # dac with adaptive lightweight mixture weights
           {
             res_dac <- dac_lgssm_lightweight(x, y[t, ], tau, lambda, sigmaY, Sigma.det, "adaptive")
           }
    )
    x <- res_dac[, 1:d]
    m[t, ] <- colMeans(x)
    v[t, ] <- colVars(x)
  }
  if(!is.null(marginals)){
    # compare marginals at last time step
    ks_dac <- apply(rbind(x, marginals), ks_dist, N = Nparticles, MARGIN = 2)
    w1_dac <- apply(rbind(x, marginals), w1_dist, N = Nparticles, MARGIN = 2)
  }
  else{
    ks_dac <- NULL
    w1_dac <- NULL
  }
  out <- list("m" = m, "v" = v, "ks" = ks_dac, "w1" = w1_dac)
  return(out)
}
