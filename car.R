# libraries
library(MASS)
library(FKF)
library(resample)
library(LaplacesDemon)

set.seed(1234)
# dimension
d <- 4
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
Time.step <- 10

# get observations
y <- lgssm_obs(mu0, Sigma0, y.coeff, x.coeff, x.error.prec, y.error.var, Time.step)

# Kalman filter
res_KF <- fkf(a0 = mu0, P0 = x.coeff%*%Sigma0 + x.error.var, dt = as.matrix(rep(0, times = d)),
              ct = as.matrix(rep(0, times = d)), Tt = x.coeff, Zt = y.coeff, HHt = x.error.var, GGt = y.error.var, yt = t(y))
true_ll <- res_KF$logLik
true_means <- res_KF$att
true_variances <- matrix(0, ncol = d, nrow = Time.step)
for (t in 1:Time.step){
  true_variances[t, ] <- diag(res_KF$Ptt[, , t])
}

# precompute determinants of precision matrices
Sigma.det <- vector(mode = "list", length = log2(d)+1)
Sigma.det[[1]] <- log(rep(1/sigmaX, times = d))
# loop over tree levels excluding leaves
nlevels <- log2(d)
nchild <- 2
for (u in 1:nlevels){
  # number of nodes at this level
  nodes <- nchild^(nlevels-u)
  # number of variables in each node
  nv <- nchild^u
  tmp <- rep(0, times = nodes)
  for (i in 1:nodes){
    ci <- child_indices(i, nv)
    # determinant of precision of variables in node i at level u
    tmp[i] <- det(x.error.prec[ci[1]:ci[2], ci[1]:ci[2]])
  }
  Sigma.det[[u+1]] <- log(tmp)
}

# DAC
Nparticles <- 1000
Nrep <- 10
se <- array(0, dim = c(Time.step, d, Nrep))
vse <- array(0, dim = c(Time.step, d, Nrep))
Zrep <- rep(0, times = Nrep)

for (j in 1:Nrep){
  x <- array(0, dim = c(Nparticles, d, Time.step+1))
  x[, , 1] <- mvrnorm(n = Nparticles, mu0, Sigma0)
  lZ <- rep(0, times = Time.step)
  m <- matrix(0, nrow = Time.step, ncol = d)
  v <- matrix(0, nrow = Time.step, ncol = d)
  for (t in 1:Time.step) {
    res_dac <- dac_car(x[, , t], y[t, ], sigmaX, sigmaY, Sigma.det)
    x[, , t+1] <- res_dac[, 1:d]
    lZ[t] <- res_dac[1, d+1]
    m[t, ] <- colMeans(x[, , t+1])
    v[t, ] <- colVars(x[, , t+1])
  }
  Zrep[j] <- sum(lZ)
  se[, , j] <- (m - t(true_means))^2
  vse[, , j] <- (v - true_variances)^2
}
mse <- apply(se, c(1,2), mean)
vmse <- apply(vse, c(1,2), mean)
boxplot(Zrep)
abline(h = true_ll, col = "red")

plot(1:Time.step, type = "l", rowMeans(mse))
plot(1:Time.step, type = "l", rowMeans(vmse))
