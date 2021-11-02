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
Time.step <- 1

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
# precompute determinants of precision matrices
Sigma.det <- vector(mode = "list", length = log2(d)+1)
Sigma.det[[1]] <- log(diag(x.error.prec))

# loop over tree levels excluding leaves
nlevels <- log2(d)
nchild <- 2
for (u in 1:nlevels){
  # number of nodes at this level
  nodes <- nchild^(nlevels-u)
  # number of variables in each node
  nvNew <- nchild^u
  tmp <- rep(0, times = nodes)
  for (i in 1:nodes){
    ci <- child_indices(i, nvNew)
    # determinant of precision of variables in node i at level u
    tmp[i] <- det(x.error.prec[ci[1]:ci[2], ci[1]:ci[2]])
  }
  Sigma.det[[u+1]] <- log(tmp)
}
# DAC
Nparticles <- 100*d
Nrep <- 10
se <- array(0, dim = c(Time.step, d, Nrep))
vse <- array(0, dim = c(Time.step, d, Nrep))
Zrep <- rep(0, times = Nrep)
tRep <- rep(0, times = Nrep)
selw <- array(0, dim = c(Time.step, d, Nrep))
vselw <- array(0, dim = c(Time.step, d, Nrep))
Zreplw <- rep(0, times = Nrep)
tReplw <- rep(0, times = Nrep)
for (j in 1:Nrep){
  x <- array(0, dim = c(Nparticles, d, Time.step+1))
  x[, , 1] <- mvrnorm(n = Nparticles, mu0, Sigma0)
  lZ <- rep(0, times = Time.step)
  m <- matrix(0, nrow = Time.step, ncol = d)
  v <- matrix(0, nrow = Time.step, ncol = d)
  xlw <- array(0, dim = c(Nparticles, d, Time.step+1))
  xlw[, , 1] <- mvrnorm(n = Nparticles, mu0, Sigma0)
  lZlw <- rep(0, times = Time.step)
  mlw <- matrix(0, nrow = Time.step, ncol = d)
  vlw <- matrix(0, nrow = Time.step, ncol = d)
  for (t in 1:Time.step) {
    # dac
    tic()
    res_dac <- dac_lgssm(x[, , t], y[t, ], tau, lambda, sigmaY, Sigma.det)
    runtime <- toc()
    tRep[j] <- tRep[j] + runtime$toc - runtime$tic
    x[, , t+1] <- res_dac[, 1:d]
    lZ[t] <- res_dac[1, d+1]
    m[t, ] <- colMeans(x[, , t+1])
    v[t, ] <- colVars(x[, , t+1])

    # ligthweight dac
    tic()
    res_dac_lw <- dac_lgssm_lightweight(x[, , t], y[t, ], tau, lambda, sigmaY, Sigma.det, ceiling(sqrt(Nparticles)))
    runtime <- toc()
    tReplw[j] <- tReplw[j] + runtime$toc - runtime$tic
    xlw[, , t+1] <- res_dac_lw[, 1:d]
    lZlw[t] <- res_dac_lw[1, d+1]
    mlw[t, ] <- colMeans(xlw[, , t+1])
    vlw[t, ] <- colVars(xlw[, , t+1])
  }
  Zrep[j] <- sum(lZ)
  se[, , j] <- (m - t(true_means))^2
  vse[, , j] <- (v - true_variances)^2
  Zreplw[j] <- sum(lZlw)
  selw[, , j] <- (mlw - t(true_means))^2
  vselw[, , j] <- (vlw - true_variances)^2
}
mse <- apply(se, c(1,2), mean)
vmse <- apply(vse, c(1,2), mean)
mselw <- apply(selw, c(1,2), mean)
vmselw <- apply(vselw, c(1,2), mean)

# MSE
plot(1:Time.step, type = "l", rowMeans(mse), col = "blue", xlab=" ", ylab=" ",
     ylim = c(0, max(rowMeans(mse), rowMeans(mselw))), cex = 1.5)
lines(1:Time.step, rowMeans(mselw), col = "red")
legend(1, 0.001, legend = c("dac", "dac-lw"), col=c("red", "blue"), lty=1, cex=0.8)

Nparticles <- 1000
M <- ceiling(sqrt(Nparticles))
x_stpf <- array(mvrnorm(n = Nparticles*M, mu0, Sigma0), dim = c(Nparticles, M, d))
res_stpf <- stpf_lgssm(x_stpf, obs, tau, lambda, sigmaY)
