# Linear Gaussian SSM -- comparison of dac and dac with mixture reweighting (both lightweight and full cost)
set.seed(1234)

# dimension
d <- 4
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
true_means <- res_KF$att
true_variances <- matrix(0, ncol = d, nrow = Time.step)
for (t in 1:Time.step){
  true_variances[t, ] <- diag(res_KF$Ptt[, , t])
}



### DAC
Nparticles <- 100*d
Nrep <- 10
se <- rep(list(array(0, dim = c(Time.step, d, Nrep))), times = 3)
vse <- rep(list(array(0, dim = c(Time.step, d, Nrep))), times = 3)
Zrep <- matrix(0, nrow = Nrep, ncol = 3)
trep <- matrix(0, nrow = Nrep, ncol = 3)

for (j in 1:Nrep){
  x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
  # dac linear cost
  tic()
  res_dac <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "lc")
  runtime <- toc()
  trep[j, 1] <- runtime$toc - runtime$tic
  Zrep[j, 1] <- sum(res_dac[, 2*d+1])
  se[1][[1]][, , j] <- (res_dac[, 1:d] - t(true_means))^2
  vse[1][[1]][, , j] <- (res_dac[, (d+1):(2*d)] - true_variances)^2

  # dac with mixture weights
  tic()
  res_dac_mix <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "mix")
  runtime <- toc()
  trep[j, 2] <- runtime$toc - runtime$tic
  Zrep[j, 2] <- sum(res_dac_mix[, 2*d+1])
  se[2][[1]][, , j] <- (res_dac_mix[, 1:d] - t(true_means))^2
  vse[2][[1]][, , j] <- (res_dac_mix[, (d+1):(2*d)] - true_variances)^2

  # dac with lightweight mixture weighting
  tic()
  res_dac_light <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "light")
  runtime <- toc()
  trep[j, 3] <- runtime$toc - runtime$tic
  Zrep[j, 3] <- sum(res_dac_light[, 2*d+1])
  se[3][[1]][, , j] <- (res_dac_light[, 1:d] - t(true_means))^2
  vse[3][[1]][, , j] <- (res_dac_light[, (d+1):(2*d)] - true_variances)^2
}
mse <- array(0, dim = c(Time.step, d, 3))
vmse <- array(0, dim = c(Time.step, d, 3))
for (j in 1:3) {
  mse[, , j] <- apply(se[j][[1]], c(1,2), mean)
  vmse[, , j] <- apply(vse[j][[1]], c(1,2), mean)
}

plot(1:Time.step, type = "l", rowMeans(mse[, , 1]), col = "blue", xlab=" ", ylab=" ", cex = 1.5,
     ylim = c(0, max(mse[, , 1])))
lines(1:Time.step, rowMeans(mse[, , 2]), col = "red")
lines(1:Time.step, rowMeans(mse[, , 3]), col = "gray")
legend(1, 0.001, legend = c("dac", "dac-lw"), col=c("red", "blue"), lty=1, cex=0.8)

plot(1:Time.step, type = "l", rowMeans(mse[, , 1]/true_variances), col = "blue", xlab=" ", ylab=" ", cex = 1.5,
     ylim = c(0, max(mse[, , 1]/true_variances)))
lines(1:Time.step, rowMeans(mse[, , 2]/true_variances), col = "red")
lines(1:Time.step, rowMeans(mse[, , 3]/true_variances), col = "gray")
