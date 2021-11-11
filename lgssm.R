set.seed(1234)
# dimension
d <- 8
# initial state
mu0 <- rep(0, times = d)
Sigma0 <- diag(x = 1, d, d)
# parameters
tau <- 1
lambda <- 1
sigmaY <- 0.5^2
# coefficient and precision of the state equation
x.coeff <- diag(x = 0.5, d, d)
x.error.prec <- diag(x = (tau + 2*lambda), d, d)
x.error.prec[row(x.error.prec) - col(x.error.prec) == 1] <- -lambda
x.error.prec[row(x.error.prec) - col(x.error.prec) == -1] <- -lambda
x.error.prec[1, 1] <- x.error.prec[1, 1] - lambda
x.error.prec[d, d] <- x.error.prec[d, d] - lambda
x.error.var <- solve(x.error.prec)

# coefficient and covariance of the observation equation
y.error.var <- diag(x = sigmaY, d, d)
y.coeff <- diag(x = 1, d, d)

# number of time steps
Time.step <- 10

# get observations
y <- lgssm_obs(mu0, Sigma0, y.coeff, x.coeff, x.error.prec, y.error.var, Time.step)

# Kalman filter
res_KF <- fkf(a0 = mu0, P0 = 0.5^2*Sigma0 + x.error.var, dt = as.matrix(rep(0, times = d)),
              ct = as.matrix(rep(0, times = d)), Tt = x.coeff, Zt = y.coeff, HHt = x.error.var, GGt = y.error.var, yt = t(y))
true_ll <- res_KF$logLik
true_means <- t(res_KF$att)
true_variances <- matrix(0, ncol = d, nrow = Time.step)
for (t in 1:Time.step){
  true_variances[t, ] <- diag(res_KF$Ptt[, , t])
}
# samples from marginals at last time step
marginals <- matrix(0, nrow = 10^5, ncol = d)
for(i in 1:d){
  marginals[, i] <- rnorm(10^5, mean = true_means[Time.step, i], sd = sqrt(true_variances[Time.step, i]))
}

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


Nparticles <- 100*d
M <- 2*d
Nrep <- 1
# dac
se_dac <- array(0, dim = c(Time.step, d, Nrep))
vse_dac <- array(0, dim = c(Time.step, d, Nrep))
Zrep_dac <- matrix(0, nrow = Time.step, ncol = Nrep)
trep_dac <- rep(0, times = Nrep)
ks_dac <- matrix(0, nrow = d, ncol = Nrep)
w1_dac <- matrix(0, nrow = d, ncol = Nrep)
# stpf
se_stpf <- array(0, dim = c(Time.step, d, Nrep))
vse_stpf <- array(0, dim = c(Time.step, d, Nrep))
Zrep_stpf <- matrix(0, nrow = Time.step, ncol = Nrep)
trep_stpf <- rep(0, times = Nrep)
ks_stpf <- matrix(0, nrow = d, ncol = Nrep)
w1_stpf <- matrix(0, nrow = d, ncol = Nrep)
for (j in 1:Nrep){
  x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
  # dac (lightweight)
  tic()
  res_dac_light <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "light")
  runtime <- toc()
  trep_dac <- runtime$toc - runtime$tic
  Zrep_dac[, j] <- res_dac_light[1, d+1]
  #se_dac[, , j] <- (colMeans(res_dac_light[, 1:d]) - true_means)^2
  #vse_dac[, , j] <- (colVars(res_dac_light[, 1:d]) - true_variances)^2
  ks_dac[, j] <- apply(rbind(res_dac_light[, 1:d], marginals), ks_dist, N = Nparticles, MARGIN = 2)
  w1_dac[, j] <- apply(rbind(res_dac_light[, 1:d], marginals), w1_dist, N = Nparticles, MARGIN = 2)

  # stpf
  x0 <- array(mvrnorm(n = Nparticles*M, mu0, Sigma0), dim = c(Nparticles, M, d))
  tic()
  res_stpf <- stpf_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y)
  runtime <- toc()
  trep_stpf <- runtime$toc - runtime$tic
  Zrep_stpf[, j] <- res_stpf[2][[1]]
  #se_stpf[, , j] <- (res_stpf[, 1:d] - true_means)^2
  #vse_stpf[, , j] <- (res_stpf[, (d+1):(2*d)] - true_variances)^2
}
mse_dac <- apply(se_dac, c(1,2), mean)
vmse_dac <- apply(vse_dac, c(1,2), mean)
mse_stpf <- apply(se_stpf, c(1,2), mean)
vmse_stpf <- apply(vse_stpf, c(1,2), mean)

