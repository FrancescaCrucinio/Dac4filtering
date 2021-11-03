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
Time.step <- 20

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

Nparticles <- 1000
M <- ceiling(sqrt(Nparticles))
Nrep <- 1
se_dac <- array(0, dim = c(Time.step, d, Nrep))
vse_dac <- array(0, dim = c(Time.step, d, Nrep))
Zrep_dac <- rep(0, times = Nrep)
tRep_dac <- rep(0, times = Nrep)
se_stpf <- array(0, dim = c(Time.step, d, Nrep))
vse_stpf <- array(0, dim = c(Time.step, d, Nrep))
Zrep_stpf <- rep(0, times = Nrep)
tRep_stpf <- rep(0, times = Nrep)
for (j in 1:Nrep){
  # dac
  x_dac <- array(0, dim = c(Nparticles, d, Time.step+1))
  x_dac[, , 1] <- mvrnorm(n = Nparticles, mu0, Sigma0)
  lZ_dac <- rep(0, times = Time.step)
  m_dac <- matrix(0, nrow = Time.step, ncol = d)
  v_dac <- matrix(0, nrow = Time.step, ncol = d)
  # stpf
  x_stpf <- array(mvrnorm(n = Nparticles*M, mu0, Sigma0), dim = c(Nparticles, M, d))
  lZ_stpf <- rep(0, times = Time.step)
  m_stpf <- matrix(0, nrow = Time.step, ncol = d)
  v_stpf <- matrix(0, nrow = Time.step, ncol = d)
  for (t in 1:Time.step) {
    # dac
    tic()
    res_dac <- dac_car_lightweight(x_dac[, , t], y[t, ], sigmaX, sigmaY, Sigma.det, M)
    runtime <- toc()
    tRep_dac[j] <- tRep_dac[j] + runtime$toc - runtime$tic
    x_dac[, , t+1] <- res_dac[, 1:d]
    lZ_dac[t] <- res_dac[1, d+1]
    m_dac[t, ] <- colMeans(x_dac[, , t+1])
    v_dac[t, ] <- colVars(x_dac[, , t+1])

    # stpf
    tic()
    res_stpf <- stpf_car(x_stpf, y[t, ], sigmaX, sigmaY)
    runtime <- toc()
    tRep_stpf[j] <- tRep_stpf[j] + runtime$toc - runtime$tic
    x_stpf <- res_stpf[1][[1]]
    lZ_stpf[t] <- res_stpf[2][[1]]
    # reshape x
    tmp_x <- matrix(rbind(x_stpf), Nparticles*M, d)
    m_stpf[t, ] <- colMeans(tmp_x)
    v_stpf[t, ] <- colVars(tmp_x)
  }
  # Zrep[j] <- sum(lZ)
  se[, , j] <- (m_dac - t(true_means))^2
  # vse[, , j] <- (v - true_variances)^2
  # Zreplw[j] <- sum(lZlw)
  # selw[, , j] <- (mlw - t(true_means))^2
  # vselw[, , j] <- (vlw - true_variances)^2
}
mse <- apply(se, c(1,2), mean)
# vmse <- apply(vse, c(1,2), mean)
# mselw <- apply(selw, c(1,2), mean)
# vmselw <- apply(vselw, c(1,2), mean)
#
# # MSE
plot(1:Time.step, type = "l", rowMeans(mse), col = "blue", xlab=" ", ylab=" ", cex = 1.5)
# lines(1:Time.step, rowMeans(mselw), col = "red")
# legend(1, 0.001, legend = c("dac", "dac-lw"), col=c("red", "blue"), lty=1, cex=0.8)
#
# boxplot(Zrep)
# abline(h = true_ll, col = "red")


