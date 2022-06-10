### Spatial model stats
# dimension
d <- 32
# parameters
sigmaX <- 1
nu <- 10
tau_diag <- 1
tau <- -0.25
# number of time steps
Time.step <- 10

Nparticles <- 500
Nrep <- 50

data_dac <- array(0, dim = c(d, d, Nrep))
mse_dac <- array(0, dim = c(d, d, Nrep))
dac_time <- rep(0, times = Nrep)
ground_truth <- unname(data.matrix(read.csv(paste0("/storage/u1693998/data/corrected_data_truth_spatial_tau_", -tau, "d", d, "ID", 1),
                                            row.names = 1, nrows=d, skip=(Time.step)*d)))
for (i in 1:Nrep){
  df <- read.csv(paste0("/storage/u1693998/results/same_nt_corrected_means_dac_spatial_d", d, "N", Nparticles, "ID", i, "step", 1), row.names = 1,
                  nrows=d, skip=(Time.step-1)*d)
  # df2 <- read.csv(paste0("/storage/u1693998/results/nt_corrected_means_dac_spatial_d", d, "N", Nparticles, "ID", i, "step", 2), row.names = 1,
  #                 nrows=d)
  # df3 <- read.csv(paste0("/storage/u1693998/results/nt_corrected_means_dac_spatial_d", d, "N", Nparticles, "ID", i, "step", 3), row.names = 1,
  #                 nrows=d, skip=(4-1)*d)
  # df4 <- read.csv(paste0("/storage/u1693998/results/corrected_means_dac_spatial_d", d, "N", Nparticles, "ID", i, "step", 4), row.names = 1,
  #                 nrows=d, skip=(2-1)*d)
  # df5 <- read.csv(paste0("/storage/u1693998/results/corrected_means_dac_spatial_d", d, "N", Nparticles, "ID", i, "step", 5), row.names = 1,
  #                 nrows=d, skip=(2-1)*d)
  data_dac[, , i] <- unname(data.matrix(df[, 1:d]))

  dac_time[i] <- df[1, d+2]
  # + df3[1, d+2] + df4[1, d+2] + df5[1, d+2]

  mse_dac[, , i] <- (ground_truth - data_dac[, , i])^2
  print(paste(i))
}
average_mse_dac <- apply(mse_dac, c(1,2), mean)

df_dac <- data.frame(average_mse_dac)
df_dac$runtime <- mean(dac_time)
write.csv(x=df_dac, file=paste0("run_dac/results/same_nt_corrected_stats_dac_spatial_tau", -tau, "d", d, "N", Nparticles))

