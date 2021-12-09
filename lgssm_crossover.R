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
m1 <- diag(x = tau, d, d)
m1[1, 1] <- tau + lambda
m2 <- diag(x = tau+lambda, d, d)
m2[row(m2) - col(m2) == 1] <- -lambda
m3 <- diag(x = tau+lambda, d, d)
m3[1, 1] <- tau
x.coeff <- 0.5*m1%*%solve(m2)
x.error.prec <- (t(m2)%*%m3%*%m2)/(tau+lambda)^2
x.error.var <- solve(x.error.prec)

# coefficient and covariance of the observation equation
y.error.var <- diag(x = sigmaY, d, d)
y.coeff <- diag(x = 1, d, d)

# number of time steps
Time.step <- 10

# get observations
y <- ssm_obs(mu0, Sigma0, y.coeff, x.coeff, x.error.prec, y.error.var, Time.step)

# Kalman filter
res_KF <- fkf(a0 = mu0, P0 = 0.5^2*t(x.coeff)%*%Sigma0%*%x.coeff + x.error.var, dt = as.matrix(rep(0, times = d)),
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

Nparticles <- 1000
M <- 2*d
Nrep <- 10
# dac
se_dac <- array(0, dim = c(Time.step, d, Nrep))
vse_dac <- array(0, dim = c(Time.step, d, Nrep))
trep_dac <- rep(0, times = Nrep)
ks_dac <- matrix(0, nrow = d, ncol = Nrep)
w1_dac <- matrix(0, nrow = d, ncol = Nrep)
# stpf
se_stpf <- array(0, dim = c(Time.step, d, Nrep))
vse_stpf <- array(0, dim = c(Time.step, d, Nrep))
trep_stpf <- rep(0, times = Nrep)
ks_stpf <- matrix(0, nrow = d, ncol = Nrep)
w1_stpf <- matrix(0, nrow = d, ncol = Nrep)
# nsmc
se_nsmc <- array(0, dim = c(Time.step, d, Nrep))
vse_nsmc <- array(0, dim = c(Time.step, d, Nrep))
trep_nsmc <- rep(0, times = Nrep)
ks_nsmc <- matrix(0, nrow = d, ncol = Nrep)
w1_nsmc <- matrix(0, nrow = d, ncol = Nrep)
for (j in 1:Nrep){
  x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
  # dac (lightweight)
  tic()
  res_dac_light <- dac_time_lgssm_crossover(tau, lambda, sigmaY, Nparticles, x0, y, method = "adaptive", marginals = marginals)
  runtime <- toc()
  trep_dac <- runtime$toc - runtime$tic
  se_dac[, , j] <- (res_dac_light$m - true_means)^2
  vse_dac[, , j] <- (res_dac_light$v - true_variances)^2
  ks_dac[, j] <- res_dac_light$ks
  w1_dac[, j] <- res_dac_light$w1

  # nsmc
  tic()
  res_nsmc <- nsmc_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, M = M, marginals = marginals)
  runtime <- toc()
  trep_nsmc <- runtime$toc - runtime$tic
  se_nsmc[, , j] <- (res_nsmc$m - true_means)^2
  vse_nsmc[, , j] <- (res_nsmc$v - true_variances)^2
  ks_nsmc[, j] <- res_nsmc$ks
  w1_nsmc[, j] <- res_nsmc$w1

  # stpf
  x0 <- array(mvrnorm(n = Nparticles*M, mu0, Sigma0), dim = c(Nparticles, M, d))
  tic()
  res_stpf <- stpf_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, marginals = marginals)
  runtime <- toc()
  trep_stpf <- runtime$toc - runtime$tic
  se_stpf[, , j] <- (res_stpf$m - true_means)^2
  vse_stpf[, , j] <- (res_stpf$v - true_variances)^2
  ks_stpf[, j] <- res_stpf$ks
  w1_stpf[, j] <- res_stpf$w1
}
mse_dac <- apply(se_dac, c(1,2), mean)
vmse_dac <- apply(vse_dac, c(1,2), mean)
mse_stpf <- apply(se_stpf, c(1,2), mean)
vmse_stpf <- apply(vse_stpf, c(1,2), mean)
mse_nsmc <- apply(se_nsmc, c(1,2), mean)
vmse_nsmc <- apply(vse_nsmc, c(1,2), mean)

plot(1:Time.step, type = "l", rowMeans(mse_dac/true_variances), col = "blue", xlab=" ", ylab=" ", cex = 1.5, ylim = c(0, max(mse_dac/true_variances)))
lines(1:Time.step, rowMeans(mse_nsmc/true_variances), col = "red")
lines(1:Time.step, rowMeans(mse_stpf/true_variances), col = "green")

plot(1:Time.step, type = "l", colMeans(ks_dac), col = "blue", xlab=" ", ylab=" ", cex = 1.5, ylim = c(0, max(colMeans(ks_stpf))))
lines(1:Time.step, colMeans(ks_nsmc), col = "red")
lines(1:Time.step, colMeans(ks_stpf), col = "green")

plot(1:Time.step, type = "l", colMeans(w1_dac), col = "blue", xlab=" ", ylab=" ", cex = 1.5, ylim = c(0,  max(colMeans(w1_stpf))))
lines(1:Time.step, colMeans(w1_nsmc), col = "red")
lines(1:Time.step, colMeans(w1_stpf), col = "green")

