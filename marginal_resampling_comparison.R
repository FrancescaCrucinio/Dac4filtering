### Linear Gaussian SSM -- comparison of dac and dac with mixture reweighting (both lightweight and full cost)
ID <- 1
set.seed(1234*ID)

d <- 128
timeinterval <- 1
# get time interval extrema
ti_begin <- 1 + (timeinterval - 1)*10
ti_end <- timeinterval*10
Time.step <- 10

filename <- paste0("/storage/u1693998/data/data_lgssm_d", d, "ID", ID)
data <- read.csv(filename, row.names = 1)

# parameters
mu0 <- rep(0, times = d)
Sigma0 <- diag(x = 1, d, d)
tau <- 1
lambda <- 1
sigmaY <- 0.5^2
y <- data.matrix(data[ti_begin:ti_end, 1:d], rownames.force = NA)
true_means <- unname(data.matrix(data[data$type == "mean", 1:d]))
rm(data)

Nparticles <- 100

if(timeinterval>1){
  filename <- paste0("/storage/u1693998/results/dac_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval-1)
  data_dac <- read.table(filename, row.names = 1)
  history <- data.matrix(data_dac)
} else{
  # initial value
  history <- mvrnorm(n = Nparticles, mu0, Sigma0)
}

# lightweight
res_dac_light <- history
tic()
for (t in ti_begin:ti_end) {
  res_dac_light <- marginal_dac_lgssm_lightweight(res_dac_light, y[t-ti_begin+1, ], tau, lambda, sigmaY)
}
runtime <- toc()
se <- (colMeans(res_dac_light) - true_means[Time.step, ])^2
df <- data.frame(t(c(se, runtime)))
df$algo <- as.factor(c("light"))
filename <- paste0("/storage/u1693998/results/light_marginal_resampling_comparison_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
write.csv(x=df, file=filename)


# adaptive lightweight
res_dac_light_ada <- history
tic()
for (t in ti_begin:ti_end) {
  res_dac_light_ada <- marginal_dac_lgssm_lightweight(res_dac_light_ada, y[t-ti_begin+1, ], tau, lambda, sigmaY, adaptive = TRUE)
}
runtime <- toc()
se <- (colMeans(res_dac_light_ada) - true_means[Time.step, ])^2
df <- data.frame(t(c(se, runtime)))
df$algo <- as.factor(c("light-ada"))
filename <- paste0("/storage/u1693998/results/ada_light_marginal_resampling_comparison_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
write.csv(x=df, file=filename)

# mix
res_dac_mix <- history
tic()
for (t in ti_begin:ti_end) {
  res_dac_mix <- marginal_dac_lgssm_mix(res_dac_mix, y[t-ti_begin+1, ], tau, lambda, sigmaY)
}
runtime <- toc()
se <- (colMeans(res_dac_mix) - true_means[Time.step, ])^2
df <- data.frame(t(c(se, runtime)))
df$algo <- as.factor(c("mix"))
filename <- paste0("/storage/u1693998/results/mix_marginal_resampling_comparison_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
write.csv(x=df, file=filename)

# linear cost
res_dac_lc <- history
tic()
for (t in ti_begin:ti_end) {
  res_dac_lc <- marginal_dac_lgssm_lc(res_dac_lc, y[t-ti_begin+1, ], tau, lambda, sigmaY)
}
runtime <- toc()
se <- (colMeans(res_dac_lc) - true_means[Time.step, ])^2
df <- data.frame(t(c(se, runtime)))
df$algo <- as.factor(c("lc"))
filename <- paste0("/storage/u1693998/results/lc_marginal_resampling_comparison_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
write.csv(x=df, file=filename)

