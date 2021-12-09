dac_time_lgssm_crossover <- function(tau, lambda, sigmaY, Nparticles, x0, y, method = "adaptive", M = NULL, marginals = NULL){
  # dimension
  d <- ncol(y)
  # time interval
  Time.step <- nrow(y)
  m <- matrix(0, nrow = Time.step, ncol = d)
  v <- matrix(0, nrow = Time.step, ncol = d)
  # initial value
  history <- array(0, dim = c(Nparticles, d, 2))
  history[, , 1] <- x0

  for (t in 1:Time.step) {
    switch(method,
           lc={
             # linear cost dac
             res_dac <- dac_lgssm_lc_crossover(history, y[t, ], tau, lambda, sigmaY)
           },
           mix={
             # dac with mixture weights
             res_dac <- dac_lgssm_crossover(history, y[t, ], tau, lambda, sigmaY)
           }, # dac with lightweight mixture weights (no adaptation)
           light={
             res_dac <- dac_lgssm_lightweight_crossover(history, y[t, ], tau, lambda, sigmaY)
           }, # dac with adaptive lightweight mixture weights
           {
             res_dac <- dac_lgssm_lightweight_crossover(history, y[t, ], tau, lambda, sigmaY, "adaptive")
           }
    )
    x <- res_dac[, 1:d]
    history[, , 2] <- history[, , 1]
    history[, , 1] <- x
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
