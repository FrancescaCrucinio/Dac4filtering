devtools::load_all("/storage/u1693998/Dac4filtering")
### Linear Gaussian SSM -- comparison of dac, stpf, nsmc
# dimension
d <- 32
# parameters
sigmaX <- 1
nu <- 10
tau <- 1/4
delta <- 1
Time.step <- 10

Nparticles <- 1000
M <- 100

data_nsmc <- array(0, dim = c(d, d, 50))
data_stpf <- array(0, dim = c(d, d, 50))
mse_nsmc <- array(0, dim = c(d, d, 50))
mse_stpf <- array(0, dim = c(d, d, 50))
nsmc_time <- rep(0, times = 50)
stpf_time <- rep(0, times = 50)
for (i in 1:50){
  ground_truth <- unname(data.matrix(read.csv(paste0("/storage/u1693998/data/data_truth_nl_tau_", tau, "d", d, "ID", i),
                                              row.names = 1, nrows=d, skip=(Time.step)*d)))
  d1 <- read.csv(paste0("/storage/u1693998/results/results/stpf_nl_iid_d", d, "N", Nparticles, "ID", i), row.names = 1,
                 nrows=d, skip=(Time.step-1)*d)
  d2 <- read.csv(paste0("/storage/u1693998/results/results/nsmc_nl_iid_d", d, "N", Nparticles, "ID", i), row.names = 1,
                nrows=d, skip=(Time.step-1)*d)
  data_stpf[, , i] <- unname(data.matrix(d1[, 1:d]))
  data_nsmc[, , i] <- unname(data.matrix(d2[, 1:d]))

  stpf_time[i] <- d1[1, d+2]
  nsmc_time[i] <- d2[1, d+2]

  mse_nsmc[, , i] <- (ground_truth - data_nsmc[, , i])^2
  mse_stpf[, , i] <- (ground_truth - data_stpf[, , i])^2
}


var_nsmc <- apply(mse_nsmc, c(1,2), var)
var_stpf <- apply(mse_stpf, c(1,2), var)

mse_nsmc <- apply(mse_nsmc, c(1,2), mean)
mse_stpf <- apply(mse_stpf, c(1,2), mean)

df_nsmc <- data.frame(rbind(mse_nsmc, var_nsmc))
df_nsmc$type <- rep(c("mse", "var"), each = d)
df_nsmc$runtime <- mean(nsmc_time)
write.csv(x=df_nsmc, file=paste0("run_dac/results/stats_nsmc_nl_iid_d", d, "N", Nparticles))

df_stpf <- data.frame(rbind(mse_stpf, var_stpf))
df_stpf$type <- rep(c("mse", "var"), each = d)
df_stpf$runtime <- mean(stpf_time)
write.csv(x=df_stpf, file=paste0("run_dac/results/stats_stpf_nl_iid_d", d, "N", Nparticles))
