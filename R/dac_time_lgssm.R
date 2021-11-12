dac_time_lgssm <- function(tau, lambda, sigmaY, Nparticles, x0, y, method = "light", M = NULL, marginals = NULL){
  # dimension
  d <- ncol(y)
  # time interval
  Time.step <- nrow(y)
  lZ <- rep(0, times = Time.step)
  m <- matrix(0, nrow = Time.step, ncol = d)
  v <- matrix(0, nrow = Time.step, ncol = d)
  # initial value
  x <- x0

  # precompute determinants of precision matrices
  Sigma.det <- vector(mode = "list", length = log2(d)+1)
  Sigma.det[[1]] <- log(diag(x.error.prec))
  # loop over tree levels excluding leaves
  nlevels <- log2(d)
  nchild <- 2
  for (u in 1:nlevels){
    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    # number of variables in each node
    nvNew <- nchild^u
    tmp <- rep(0, times = nodes)
    for (i in 1:nodes){
      ci <- child_indices(i, nvNew)
      # determinant of precision of variables in node i at level u
      tmp[i] <- det(x.error.prec[ci[1]:ci[2], ci[1]:ci[2]])
    }
    Sigma.det[[u+1]] <- log(tmp)
  }
  if(is.null(M)) {
    # parameter for lightweight mixture
    M <- ceiling(sqrt(Nparticles))
  }
  for (t in 1:Time.step) {
    switch(method,
           lc={
             # linear cost dac
             res_dac <- dac_lgssm_lc(x, y[t, ], tau, lambda, sigmaY, Sigma.det)
           },
           mix={
             # dac with mixture weights
             res_dac <- dac_lgssm(x, y[t, ], tau, lambda, sigmaY, Sigma.det)
           }, # dac with lightweight mixture weights
           {
             res_dac <- dac_lgssm_lightweight(x, y[t, ], tau, lambda, sigmaY, Sigma.det, M)
           }
    )
    x <- res_dac[, 1:d]
    lZ[t] <- res_dac[1, d+1]
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
  out <- list("m" = m, "v" = v, "lZ" = lZ, "ks" = ks_dac, "w1" = w1_dac)
  return(out)
}
