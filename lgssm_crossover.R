# devtools::load_all("/storage/u1693998/Dac4filtering")
# ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
ID <- 1
set.seed(1234*ID)
# dimension
d <- 32
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
Time.step <- 1

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
heatmap(cov2cor(res_KF$Ptt[, , Time.step]), Colv = NA, Rowv = NA)
# samples from marginals at last time step
marginals <- matrix(0, nrow = 10^5, ncol = d)
for(i in 1:d){
  marginals[, i] <- rnorm(10^5, mean = true_means[Time.step, i], sd = sqrt(true_variances[Time.step, i]))
}

Nparticles <- 100
M <- 100

x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
# dac (lightweight adaptive)
tic()
res_dac <- dac_lgssm_lightweight(x0, y[t, ], tau, lambda, sigmaY, "adaptive")
runtime <- toc()
cov_res_dac <- cov(res_dac)
heatmap(cov_res_dac, Colv = NA, Rowv = NA, scale="column")

# dac (lightweight)
tic()
res_dac_light <- dac_lgssm_lightweight(x0, y[t, ], tau, lambda, sigmaY)
runtime <- toc()
cov_res_dac_light <- cov(res_dac_light)
heatmap(cov_res_dac_light, Colv = NA, Rowv = NA, scale="column")
# nsmc
tic()
res_nsmc <- nsmc_lgssm(x0, y[t, ], tau, lambda, sigmaY, M)
runtime <- toc()
cov_res_nsmc <- cov(res_nsmc)
heatmap(cov_res_nsmc, Colv = NA, Rowv = NA, scale="column")
# stpf
x0 <- array(mvrnorm(n = Nparticles*M, mu0, Sigma0), dim = c(Nparticles, M, d))
tic()
res_stpf <- stpf_lgssm(x0, y[t, ], tau, lambda, sigmaY)
runtime <- toc()
cov_res_stpf <- cov(matrix(res_stpf, ncol = d, nrow = Nparticles*M))
heatmap(cov_res_stpf, Colv = NA, Rowv = NA, scale="column")

plot(1:Time.step, rowMeans(true_means), type = "l")
lines(1:Time.step, rowMeans(res_dac_light$m), type = "l", col = "red")
lines(1:Time.step, rowMeans(res_nsmc$m), type = "l", col = "green")
lines(1:Time.step, rowMeans(res_stpf$m), type = "l", col = "blue")

