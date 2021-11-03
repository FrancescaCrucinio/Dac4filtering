# Linear Gaussian SSM -- comparison of dac and dac with mixture reweighting (both lightweight and full cost)
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
true_means <- res_KF$att
true_variances <- matrix(0, ncol = d, nrow = Time.step)
for (t in 1:Time.step){
  true_variances[t, ] <- diag(res_KF$Ptt[, , t])
}



### DAC
Nparticles <- 100*d
Nrep <- 5
se <- rep(list(array(0, dim = c(Time.step, d, Nrep))), times = 3)
vse <- rep(list(array(0, dim = c(Time.step, d, Nrep))), times = 3)
Zrep <- matrix(0, nrow = Nrep, ncol = 3)
tRep <- matrix(0, nrow = Nrep, ncol = 3)
