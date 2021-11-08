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
Nparticles <- 1000
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

g <- as.factor(rep(c("dac", "mix", "light"), each = Time.step))
# MSE
df <- data.frame(x = rep(1:Time.step, times = 3), y = c(rowMeans(mse[, , 1]), rowMeans(mse[, , 2]), rowMeans(mse[, , 3])), g)
ggplot(data = df, aes(x = x, y= y, group = g)) +
  geom_line(aes(x = x, y= y, color = g), size = 2) +
  scale_y_continuous(trans='log10') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=15))
# RMSE
df <- data.frame(x = rep(1:Time.step, times = 3), y = c(rowMeans(mse[, , 1]/true_variances),
                                                        rowMeans(mse[, , 2]/true_variances), rowMeans(mse[, , 3]/true_variances)), g)
ggplot(data = df, aes(x = x, y= y, group = g)) +
  geom_line(aes(x = x, y= y, color = g), size = 2) +
  scale_y_continuous(trans='log10') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=15))

df_dac <- data.frame()
df_dac_mix <- data.frame()
df_dac_light <- data.frame()

res <- bench::mark("dac" = {
  x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
  res_dac <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "lc")
  lZ <- res_dac[, 2*d+1]
  se <- (res_dac[, 1:d] - t(true_means))^2
  vse <- (res_dac[, (d+1):(2*d)] - true_variances)^2
  df_dac <- data.frame(rbind(df_dac, cbind(se, vse, lZ)))
  df_dac
},
"mix" = {
  x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
  res_dac_mix <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "mix")
  lZ <- res_dac_mix[, 2*d+1]
  se <- (res_dac_mix[, 1:d] - t(true_means))^2
  vse <- (res_dac_mix[, (d+1):(2*d)] - true_variances)^2
  df_dac_mix <- data.frame(rbind(df_dac_mix, cbind(se, vse, lZ)))
},
"light" = {
  x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
  res_dac_light <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "light")
  lZ <- res_dac_light[, 2*d+1]
  se <- (res_dac_light[, 1:d] - t(true_means))^2
  vse <- (res_dac_light[, (d+1):(2*d)] - true_variances)^2
  df_dac_light <- data.frame(rbind(df_dac_light, cbind(se, vse, lZ)))
},
memory = capabilities("profmem"),
check = FALSE,
iterations = 9,
filter_gc = TRUE)
