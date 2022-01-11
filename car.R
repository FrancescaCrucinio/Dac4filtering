set.seed(1234)
# dimension
d <- 8
# initial state
mu0 <- rep(0, times = d)
Sigma0 <- diag(x = 1, d, d)
# parameters
sigmaX <- 1
sigmaY <- 1

# coefficient and covariance of the state equation
m1 <- matrix(0, d, d)
m1[lower.tri(m1, diag = FALSE)] <- -1/d
diag(m1) <- 1
m2 <- matrix(0, d, d)
m2[upper.tri(m1, diag = TRUE)] <- 1/d
x.coeff <- solve(m1)%*%m2
x.error.prec <- sigmaX^(-2) * t(m1)%*%m1
x.error.var <- solve(x.error.prec)

# coefficient and covariance of the observation equation
y.error.var <- diag(sigmaY, d, d)
y.coeff <- diag(1, d, d)


# number of time steps
Time.step <- 1

# get observations
y <- ssm_obs(mu0, Sigma0, y.coeff, x.coeff, x.error.prec, y.error.var, Time.step)

# Kalman filter
res_KF <- fkf(a0 = mu0, P0 = x.coeff%*%Sigma0 + x.error.var, dt = as.matrix(rep(0, times = d)),
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
df <- data.frame()

x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
# dac (lightweight)
tic()
res_dac_light <- dac_time_car(sigmaX, sigmaY, Nparticles, x0, y, marginals = marginals)
runtime <- toc()
rse <- (res_dac_light$m - true_means)^2/true_variances
df <- data.frame(rbind(df, cbind(t(rse), res_dac_light$w1, res_dac_light$ks, rep(runtime$toc[[1]] - runtime$tic[[1]], times = d))))
# nsmc
tic()
res_nsmc <- nsmc_time_car(sigmaX, sigmaY, Nparticles, x0, y, M = M, marginals = marginals)
runtime <- toc()
rse <- (res_nsmc$m - true_means)^2/true_variances
df <- data.frame(rbind(df, cbind(t(rse), res_nsmc$w1, res_nsmc$ks, rep(runtime$toc[[1]] - runtime$tic[[1]], times = d))))
# stpf
x0 <- array(mvrnorm(n = Nparticles*M, mu0, Sigma0), dim = c(Nparticles, M, d))
tic()
res_stpf <- stpf_time_car(sigmaX, sigmaY, Nparticles, x0, y, marginals = marginals)
runtime <- toc()
rse <- (res_stpf$m - true_means)^2/true_variances
df <- data.frame(rbind(df, cbind(t(rse), res_stpf$w1, res_stpf$ks, rep(runtime$toc[[1]] - runtime$tic[[1]], times = d))))

df$algo <- as.factor(rep(c("dac", "nsmc", "stpf"), each = d))
