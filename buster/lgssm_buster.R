# libraries
library(MASS)
library(FKF)
library(LaplacesDemon)
library(resample)
library(foreach)
library(doParallel)
source("utils.R")
source("dac_lgssm.R")
source("mult_resample.R")
source("lgssm_obs.R")

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
Nrep <- 100

registerDoParallel(10)
res <- foreach (j=1:Nrep, .packages= c('MASS', 'resample', 'tictoc'), .combine='rbind',
                .multicombine=TRUE, .inorder = FALSE,
                .init=list(list(), list(), vector(), list(), list(), vector())) %dopar% {
                  x <- array(0, dim = c(Nparticles, d, Time.step+1))
                  x[, , 1] <- mvrnorm(n = Nparticles, mu0, Sigma0)
                  m <- matrix(0, nrow = Time.step, ncol = d)
                  v <- matrix(0, nrow = Time.step, ncol = d)
                  xlw <- array(0, dim = c(Nparticles, d, Time.step+1))
                  xlw[, , 1] <- mvrnorm(n = Nparticles, mu0, Sigma0)
                  mlw <- matrix(0, nrow = Time.step, ncol = d)
                  vlw <- matrix(0, nrow = Time.step, ncol = d)
                  runt <- 0
                  runtlw <- 0
                  for (t in 1:Time.step) {
                    tic()
                    res_dac <- dac_lgssm(x[, , t], y[t, ], tau, lambda, sigmaY, Sigma.det)
                    runtime <- toc()
                    runt <- runt + runtime$toc - runtime$tic
                    x[, , t+1] <- res_dac[, 1:d]
                    m[t, ] <- colMeans(x[, , t+1])
                    v[t, ] <- colVars(x[, , t+1])
                    tic()
                    res_dac_lw <- dac_lgssm_lightweight(x[, , t], y[t, ], tau, lambda, sigmaY, Sigma.det, ceiling(sqrt(Nparticles)))
                    runtime <- toc()
                    runtlw <- runtlw + runtime$toc - runtime$tic
                    xlw[, , t+1] <- res_dac_lw[, 1:d]
                    mlw[t, ] <- colMeans(xlw[, , t+1])
                    vlw[t, ] <- colVars(xlw[, , t+1])
                  }
                  list((m - t(true_means))^2, (v - true_variances)^2, runt, (mlw - t(true_means))^2, (vlw - true_variances)^2, runtlw)
                }
mse <- apply(array(unlist(res[, 1]), dim = c(Time.step, d, Nrep)), c(1,2), mean)
vmse <- apply(array(unlist(res[, 2]), dim = c(Time.step, d, Nrep)), c(1,2), mean)
mselw <- apply(array(unlist(res[, 4]), dim = c(Time.step, d, Nrep)), c(1,2), mean)
vmselw <- apply(array(unlist(res[, 5]), dim = c(Time.step, d, Nrep)), c(1,2), mean)

df <- data.frame(cbind(mse, mselw))
write.csv(df, file = "mse.csv", row.names = FALSE)

df <- data.frame(cbind(vmse, vmselw))
write.csv(df, file = "vmse.csv", row.names = FALSE)

df <- data.frame(cbind(unlist(res[, 3]), unlist(res[, 6])))
write.csv(df, file = "runtime.csv", row.names = FALSE)
