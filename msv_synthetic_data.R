# Load parameters
load("~/Documents/Dac4filtering/msv_parameters.RData")
# first semester of 2013
set.seed(1234)
Time.step <- 122
p <- 26
d <- 30
x <- matrix(0, nrow = Time.step+1, ncol = d)
y <- matrix(0, nrow = Time.step, ncol = p)
# initial state
x[1, ] <- mvrnorm(n = 1, mu, diag(Sigma))
# loop over time
for(i in 2:(Time.step+1)){
  x[i, ] <- mu + diag(Phi)%*%(x[i-1, ]-mu) + mvrnorm(n = 1, rep(0, d), diag(Sigma))
  cov_matrix <- Lambda %*% diag(exp(x[i, (p+1):d])) %*% t(Lambda) + diag(exp(x[i, 1:p]))
  y[i-1, ] <- mvrnorm(n = 1, rep(0, p), cov_matrix)
}
write.csv(x=data.frame(y), file="data/synthetic_data_msv_y", row.names = FALSE)
write.csv(x=data.frame(x), file="data/synthetic_data_msv_x", row.names = FALSE)
