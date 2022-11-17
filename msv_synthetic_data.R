# Load parameters
# first semester of 2013
set.seed(7777)
Time.step <- 100
d <- 4
# model parameters
phi <- 0.91
SigmaU <- matrix(0.3, nrow =d, ncol = d)
diag(SigmaU) <- 1
SigmaV <- matrix(0.2, nrow =d, ncol = d)
diag(SigmaV) <- 1
SigmaUV <- matrix(-0.1, nrow =d, ncol = d)
diag(SigmaUV) <- -0.2
Sigma0 <- SigmaU/(1-phi^2)
SigmaX <- SigmaU - SigmaUV %*% solve(SigmaV) %*% SigmaUV
# generate data
x <- matrix(0, nrow = Time.step, ncol = d)
y <- matrix(0, nrow = Time.step, ncol = d)
# initial state
x[1, ] <- mvrnorm(n = 1, rep(0, d), Sigma0)
cov_matrix <- diag(sqrt(exp(x[1, ]))) %*% SigmaV %*% diag(sqrt(exp(x[1, ])))
y[1, ] <- mvrnorm(n = 1, rep(0, d), cov_matrix)
# loop over time
for(i in 2:Time.step){
  mu_i <- phi*x[i-1,] + SigmaUV %*% solve(SigmaV) %*% diag(1/sqrt(exp(x[i-1,]))) %*% y[i-1, ]
  x[i, ] <- mu_i + mvrnorm(n = 1, rep(0, d), SigmaX)
  cov_matrix <- diag(sqrt(exp(x[i, ]))) %*% SigmaV %*% diag(sqrt(exp(x[i, ])))
  y[i, ] <- mvrnorm(n = 1, rep(0, d), cov_matrix)
}
write.csv(x=data.frame(y), file="data/synthetic_data_msv_nofactor_y", row.names = FALSE)
write.csv(x=data.frame(x), file="data/synthetic_data_msv_nofactor_x", row.names = FALSE)
