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
iterations = Nrep,
filter_gc = TRUE)

mse <- rep(0, times = 3*Nrep)
rmse <- rep(0, times = 3*Nrep)
for (j in 1:Nrep) {
  mse[j] <- rowMeans(df_dac[j*Time.step, 1:d])
  rmse[j] <- rowMeans(df_dac[j*Time.step, 1:d]/true_variances)
  mse[Nrep+j] <- rowMeans(df_dac_mix[j*Time.step, 1:d])
  rmse[Nrep+j] <- rowMeans(df_dac_mix[j*Time.step, 1:d]/true_variances)
  mse[2*Nrep+j] <- rowMeans(df_dac_light[j*Time.step, 1:d])
  rmse[2*Nrep+j] <- rowMeans(df_dac_light[j*Time.step, 1:d]/true_variances)
}

g <- as.factor(rep(c("dac", "mix", "light"), each = Nrep))
times <- rep(res$total_time/res$n_itr, each = Nrep)
memory <- rep(res$mem_alloc/1e6, each = Nrep)
# time
df <- data.frame(x1 = times, x2 = memory, y = mse, g)
ggplot(data = df, aes(x = x1, y = y, group = g, fill = g)) +
  geom_boxplot(aes(x = x1, y= y), coef = 6) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=15))
# memory
ggplot(data = df, aes(x = x2, y = y, group = g, fill = g)) +
  geom_boxplot(aes(x = x2, y= y), coef = 6) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=15))
