devtools::load_all("/storage/u1693998/Dac4filtering")
### Linear Gaussian SSM -- comparison of dac, stpf, nsmc
# dimension
d <- 32
# parameters
sigmaX <- 1
nu <- 10
tau <- 1/4
delta <- 1
Time.step <- 100

Nparticles <- 100
M <- 100

data_nsmc <- array(0, dim = c(d, d, 50))
data_stpf <- array(0, dim = c(d, d, 50))
mse_nsmc <- array(0, dim = c(d, d, 50))
mse_stpf <- array(0, dim = c(d, d, 50))
for (ID in 1:50){
  ground_truth <- unname(data.matrix(read.csv(paste0("/storage/u1693998/data/data_truth_nl_tau_", tau, "d", d, "ID", ID),
                                              row.names = 1, nrows=d, skip=(Time.step-1)*d)))
  data_stpf[, , i] <- unname(data.matrix(read.csv(paste0("/storage/u1693998/results/results/stpf_nl_iid_d", d, "N", Nparticles, "ID", ID), row.names = 1,
                                           nrows=d, skip=(Time.step-1)*d)))
  data_nsmc[, , i] <- unname(data.matrix(read.csv(paste0("/storage/u1693998/results/results/nsmc_nl_iid_d", d, "N", Nparticles, "ID", ID), row.names = 1,
                                           nrows=d, skip=(Time.step-1)*d)))

  mse_nsmc[, , i] <- (ground_truth - data_nsmc[, , i])^2
  mse_stpf[, , i] <- (ground_truth - data_stpf[, , i])^2
}

var_nsmc <- apply(data_nsmc, c(1,2), var)
var_stpf <- apply(data_stpf, c(1,2), var)

