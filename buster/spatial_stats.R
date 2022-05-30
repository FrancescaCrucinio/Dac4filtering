### Spatial model stats
set.seed(1234*ID)
# dimension
d <- 32
# parameters
sigmaX <- 1
nu <- 10
tau_diag <- 1
tau <- -0.25
# number of time steps
Time.step <- 10
timeinterval <- 1
ti_begin <- 1 + (timeinterval - 1)*1
ti_end <- timeinterval*1

Nparticles <- 100
Nrep <- 50

data_dac <- array(0, dim = c(d, d, 50))
mse_dac <- array(0, dim = c(d, d, 50))
dac_time <- rep(0, times = 50)
for (i in 1:50){
  ground_truth <- unname(data.matrix(read.csv(paste0("/storage/u1693998/data/data_truth_spatial_tau_", -tau, "d", d, "ID", i),
                                              row.names = 1, nrows=d, skip=(Time.step)*d)))
  df <- read.csv(paste0("/storage/u1693998/results/dac_spatial_d", d, "N", Nparticles, "ID", i), row.names = 1,
                 nrows=d, skip=(Time.step-1)*d)
  data_dac[, , i] <- unname(data.matrix(df[, 1:d]))

  dac_time[i] <- df[1, d+2]

  mse_dac[, , i] <- (ground_truth - data_dac[, , i])^2
}

average_mse_dac <- apply(mse_dac, c(1,2), mean)

df_dac <- data.frame(average_mse_dac)
df_dac$runtime <- mean(dac_time)
write.csv(x=df_nsmc, file=paste0("run_dac/results/stats_dac_spatial_tau", -tau, "d", d, "N", Nparticles))

