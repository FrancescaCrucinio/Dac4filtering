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
true_means <- res_KF$att
true_variances <- matrix(0, ncol = d, nrow = Time.step)
for (t in 1:Time.step){
  true_variances[t, ] <- diag(res_KF$Ptt[, , t])
}
