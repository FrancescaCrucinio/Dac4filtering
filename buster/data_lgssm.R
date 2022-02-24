# devtools::load_all("/storage/u1693998/Dac4filtering")

# ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
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
Time.step <- 100

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

filename <- paste0("data/data_lgssm_d", d, "ID", ID)
df <- data.frame(rbind(y, true_means, true_variances))
df$type <- as.factor(rep(c("obs", "mean", "var"), each = Time.step))
write.csv(x=df, file=filename)
