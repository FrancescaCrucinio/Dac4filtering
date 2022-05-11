### Linear Gaussian SSM -- comparison of dac, stpf, nsmc
ID <- 1
set.seed(1234*ID)
# dimension
d <- 2048
# initial state
mu0 <- rep(0, times = d)
Sigma0 <- diag(x = 1, d, d)
# parameters
tau <- 1
lambda <- 1
sigmaY <- 0.5^2
Time.step <- 100

Nparticles <- 1000
M <- 100

# get data
filename <- paste0("data/data_lgssm_d", d, "ID", ID)
true_means <- unname(data.matrix(read.csv(filename, row.names = 1, nrows=Time.step, skip=Time.step)))[, 1:d]
true_variances <- unname(data.matrix(read.csv(filename, row.names = 1, nrows=Time.step, skip=2*Time.step)))[, 1:d]

marginals <- matrix(0, nrow = 10^5, ncol = d)
for(i in 1:d){
  marginals[, i] <- rnorm(10^5, mean = true_means[Time.step, i], sd = sqrt(true_variances[Time.step, i]))
}

df <- data.frame()
data_dac <- read.csv(paste0("/data/results/dac_lgssm_d", d, "N", Nparticles, "ID", ID), row.names = 1)
data_nsmc <- read.csv(paste0("/data/results/nsmc_lgssm_d", d, "N", Nparticles, "ID", ID), row.names = 1)
data_stpf <- read.csv(paste0("/data/results/stpf_lgssm_d", d, "N", Nparticles, "ID", ID), row.names = 1)

rse_dac <- (unname(matrix(data_dac$V1, ncol = d, byrow = TRUE)) - true_means)^2/true_variances
rse_nsmc <- (unname(matrix(data_nsmc$V1, ncol = d, byrow = TRUE))  - true_means)^2/true_variances
rse_stpf <- (unname(matrix(data_stpf$V1, ncol = d, byrow = TRUE))  - true_means)^2/true_variances


last_t_nsmc <- unname(data.matrix(read.table(paste0("/data/results/t100_nsmc_lgssm_d", d, "N", Nparticles, "ID", ID), row.names = 1)))
last_t_dac <- unname(data.matrix(read.table(paste0("/data/results/t100_dac_lgssm_d", d, "N", Nparticles, "ID", ID), row.names = 1)))
last_t_stpf <- unname(data.matrix(read.table(paste0("/data/results/t100_stpf_lgssm_d", d, "N", Nparticles, "ID", ID), row.names = 1)))

res_dac_ks <- apply(rbind(last_t_dac, marginals), ks_dist, N = Nparticles, MARGIN = 2)
res_dac_w1 <- apply(rbind(last_t_dac, marginals), w1_dist, N = Nparticles, MARGIN = 2)

res_nsmc_ks <- apply(rbind(last_t_nsmc, marginals), ks_dist, N = Nparticles, MARGIN = 2)
res_nsmc_w1 <- apply(rbind(last_t_nsmc, marginals), w1_dist, N = Nparticles, MARGIN = 2)

res_stpf_ks <- apply(rbind(last_t_stpf, marginals), ks_dist, N = Nparticles, MARGIN = 2)
res_stpf_w1 <- apply(rbind(last_t_stpf, marginals), w1_dist, N = Nparticles, MARGIN = 2)

df1 <- data.frame(cbind(t(rse_dac), res_dac_w1, res_dac_ks, rep(sum(data_dac$runtime), times = d)))
colnames(df1)[(Time.step+1):(Time.step+2)] <- c("V101", "V102")
df2 <- data.frame(rbind(df, cbind(t(rse_nsmc), res_nsmc_w1, res_nsmc_ks, rep(sum(data_nsmc$runtime), times = d))))
colnames(df2)[(Time.step+1):(Time.step+2)] <- c("V101", "V102")
df3 <- data.frame(rbind(df, cbind(t(rse_stpf), res_stpf_w1, res_stpf_ks, rep(sum(data_stpf$runtime), times = d))))
colnames(df2)[(Time.step+1):(Time.step+2)] <- c("V101", "V102")

df <- rbind(df1, df2, df3)
df$algo <- as.factor(rep(c("dac-light", "nsmc", "stpf"), each = d))


filename <- paste0("run_dac/results/lgssm_d", d, "N", Nparticles, "ID", ID)
write.csv(x=df, file=filename)
