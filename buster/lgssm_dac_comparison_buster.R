devtools::load_all("/storage/u1693998/Dac4filtering")
### Linear Gaussian SSM -- comparison of dac and dac with mixture reweighting (both lightweight and full cost)
ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
set.seed(1234*ID)
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
Nparticles <- 100
df <- data.frame()

# NO CROSSOVER
# dac
x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
tic()
res_dac <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "lc")
runtime <- toc()
se <- (res_dac$m[Time.step, ] - true_means[Time.step, ])^2
df <- data.frame(rbind(df, t(c(se, runtime$toc - runtime$tic))))
# mix
tic()
res_dac_mix <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "mix")
runtime <- toc()
se <- (res_dac_mix$m[Time.step, ] - true_means[Time.step, ])^2
df <- data.frame(rbind(df, t(c(se, runtime$toc - runtime$tic))))
# light
tic()
res_dac_light <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "light")
runtime <- toc()
se <- (res_dac_light$m[Time.step, ] - true_means[Time.step, ])^2
df <- data.frame(rbind(df, t(c(se, runtime$toc - runtime$tic))))
# adaptive light
tic()
res_dac_light_ada <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "ada")
runtime <- toc()
se <- (res_dac_light_ada$m[Time.step, ] - true_means[Time.step, ])^2
df <- data.frame(rbind(df, t(c(se, runtime$toc - runtime$tic))))
# CROSSOVER
# dac
x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
# tic()
# res_dac <- dac_time_lgssm_crossover(tau, lambda, sigmaY, Nparticles, x0, y, method = "lc")
# runtime <- toc()
# se <- (res_dac$m[Time.step, ] - true_means[Time.step, ])^2
# df <- data.frame(rbind(df, t(c(se, runtime$toc - runtime$tic))))
# # mix
# tic()
# res_dac_mix <- dac_time_lgssm_crossover(tau, lambda, sigmaY, Nparticles, x0, y, method = "mix")
# runtime <- toc()
# se <- (res_dac_mix$m[Time.step, ] - true_means[Time.step, ])^2
# df <- data.frame(rbind(df, t(c(se, runtime$toc - runtime$tic))))
# # light
# tic()
# res_dac_light <- dac_time_lgssm_crossover(tau, lambda, sigmaY, Nparticles, x0, y, method = "light")
# runtime <- toc()
# se <- (res_dac_light$m[Time.step, ] - true_means[Time.step, ])^2
# df <- data.frame(rbind(df, t(c(se, runtime$toc - runtime$tic))))
# adaptive light
tic()
res_dac_light_ada <- dac_time_lgssm_crossover(tau, lambda, sigmaY, Nparticles, x0, y, method = "ada")
runtime <- toc()
se <- (res_dac_light_ada$m[Time.step, ] - true_means[Time.step, ])^2
df <- data.frame(rbind(df, t(c(se, runtime$toc - runtime$tic))))

# df$algo <- as.factor(c("dac", "mix", "light", "light_ada", "dac", "mix", "light", "light_ada"))
# df$algo <- as.factor(c("dac", "light", "light_ada", "dac", "light", "light_ada"))
# df$mutation <- as.factor(c(0, 0, 0, 0, 1, 1, 1, 1))
# df$mutation <- as.factor(c(0, 0, 0, 1, 1, 1))

df$algo <- as.factor(c("light_ada_nothreshold", "light_ada_nothreshold"))
df$mutation <- as.factor(c(0, 1))
filename <- paste0("results/additional_resampling_comparison_d", d, "N", Nparticles, "ID", ID)
write.csv(x=df, file=filename)
