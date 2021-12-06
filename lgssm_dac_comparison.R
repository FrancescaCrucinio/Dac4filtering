### Linear Gaussian SSM -- comparison of dac and dac with mixture reweighting (both lightweight and full cost)

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
res_KF <- fkf(a0 = mu0, P0 = 0.5^2*Sigma0 + x.error.var, dt = as.matrix(rep(0, times = d)),
              ct = as.matrix(rep(0, times = d)), Tt = x.coeff, Zt = y.coeff, HHt = x.error.var, GGt = y.error.var, yt = t(y))
true_ll <- res_KF$logLik
true_means <- t(res_KF$att)
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
df_dac_light_ada <- data.frame()

res <- bench::mark("dac" = {
  x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
  res_dac <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "lc")
  # lZ <- res_dac$lZ
  se <- (res_dac$m[Time.step, ] - true_means[Time.step, ])^2
 # vse <- (res_dac$v - true_variances)^2
  df_dac <- data.frame(rbind(df_dac, t(se)))
  df_dac
},
"mix" = {
  x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
  res_dac_mix <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "mix")
  # lZ <- res_dac_mix$lZ
  se <- (res_dac_mix$m[Time.step, ] - true_means[Time.step, ])^2
  # vse <- (res_dac_mix$v - true_variances)^2
  df_dac_mix <- data.frame(rbind(df_dac_mix, t(se)))
},
"light" = {
  x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
  res_dac_light <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "light")
  # lZ <- res_dac_light$lZ
  se <- (res_dac_light$m[Time.step, ] - true_means[Time.step, ])^2
  # vse <- (res_dac_light$v - true_variances)^2
  df_dac_light <- data.frame(rbind(df_dac_light, t(se)))
},
"ada" = {
  x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
  res_dac_light_ada <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "ada")
  # lZ <- res_dac_light_ada$lZ
  se <- (res_dac_light_ada$m[Time.step, ] - true_means[Time.step, ])^2
  # vse <- (res_dac_light_ada$v - true_variances)^2
  df_dac_light_ada <- data.frame(rbind(df_dac_light_ada, t(se)))
},
memory = capabilities("profmem"),
check = FALSE,
iterations = Nrep,
filter_gc = TRUE)

df <- rbind(df_dac, df_dac_mix, df_dac_light, df_dac_light_ada)
df$algo <- as.factor(rep(c("dac", "mix", "light", "light_ada"), each = nrow(df)/4))
df$runtime <- rep(res$total_time/res$n_itr, each = nrow(df)/4)
df$memory <-  rep(res$mem_alloc/1e6, each = nrow(df)/4)

ID <- 1
ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
print(ID)
filename <- paste0("resampling_comparison_d", d, "N", Nparticles, "ID", ID)write.csv(x=df, file=filename)
