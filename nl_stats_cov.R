devtools::load_all("/storage/u1693998/Dac4filtering")
### Linear Gaussian SSM -- comparison of dac, stpf, nsmc
# dimension
d <- 32
# parameters
sigmaX <- 1
nu <- 10
tau <- 1/4
delta <- 1
timeinterval <- 2
Time.step <- 100

Nparticles <- 100
M <- 100

data_dac <- array(0, dim = c(d, d, 50))
mse_dac <- array(0, dim = c(d, d, 50))
dac_time <- rep(0, times = 50)
for (i in 1:50){
  ground_truth <- unname(data.matrix(read.csv(paste0("/storage/u1693998/data/data_truth_nl_tau_", tau, "d", d, "ID", i),
                                              row.names = 1, nrows=d, skip=(10)*d)))
  runtime <- 0
  for (t in 1:timeinterval) {
    d1 <- read.csv(paste0("/storage/u1693998/results/results/dac_nl_cov_d", d, "N", Nparticles, "ID", i, "step", timeinterval), row.names = 1,
                   nrows=d, skip=(5-1)*d)
    runtime <- runtime + d1[1, d+2]
  }

  data_dac[, , i] <- unname(data.matrix(d1[, 1:d]))

  dac_time[i] <- runtime

  mse_dac[, , i] <- (ground_truth - data_dac[, , i])^2
}
var_dac <- apply(mse_dac, c(1,2), var)

mse_dac <- apply(mse_dac, c(1,2), mean)


df_dac <- data.frame(rbind(mse_dac, var_dac))
df_dac$type <- rep(c("mse", "var"), each = d)
df_dac$runtime <- mean(dac_time)
write.csv(x=df_dac, file=paste0("run_dac/results/stats_dac_nl_cov_d", d, "N", Nparticles))
