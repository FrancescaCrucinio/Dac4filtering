# Load parameters
load("~/Documents/Dac4filtering/data/msv_parameters.RData")
# first semester of 2013
set.seed(1234)
Time.step <- 122
p <- 26
d <- 30
x <- matrix(0, nrow = Time.step+1, ncol = d)
y <- matrix(0, nrow = Time.step, ncol = p)
# initial state
x[1, ] <- mvrnorm(n = 1, mu, diag(Sigma^2))
# loop over time
for(i in 2:(Time.step+1)){
  x[i, ] <- mu + diag(Phi)%*%(x[i-1, ]-mu) + mvrnorm(n = 1, rep(0, d), diag(Sigma^2))
  cov_matrix <- Lambda %*% diag(exp(x[i, (p+1):d])) %*% t(Lambda) + diag(exp(x[i, 1:p]))
  y[i-1, ] <- mvrnorm(n = 1, rep(0, p), cov_matrix)
}
load("/Users/francescacrucinio/Documents/Dac4filtering/data/exrates.RData")
colnames(y) <- colnames(dat)
write.csv(x=data.frame(y), file="data/synthetic_data_msv_y", row.names = FALSE)
write.csv(x=data.frame(x), file="data/synthetic_data_msv_x", row.names = FALSE)

# run factorstochvol
res <- fsvsample(y, factors = 4, draws = 50000, zeromean = TRUE, thin = 10, quiet = TRUE, restric = "auto")
cov_n <- rowMeans(covmat(res)[, , , 1], dims = 2)
colnames(cov_n) <- colnames(dat)
rownames(cov_n) <- colnames(dat)
write.csv(x=data.frame(cov_n), file="data/synthetic_data_msv_cov_mat", row.names = FALSE)


