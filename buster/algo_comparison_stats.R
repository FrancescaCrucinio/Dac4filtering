devtools::load_all("/storage/u1693998/Dac4filtering")
### Linear Gaussian SSM -- comparison of dac, stpf, nsmc
ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
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


filename <- paste0("/storage/u1693998/data/data_lgssm_d", d, "ID", ID)
true_means <- unname(data.matrix(read.csv(filename, row.names = 1, nrows=Time.step, skip=Time.step)))[, 1:d]
true_variances <- unname(data.matrix(read.csv(filename, row.names = 1, nrows=Time.step, skip=2*Time.step)))[, 1:d]

marginals <- matrix(0, nrow = 10^5, ncol = d)
for(i in 1:d){
  marginals[, i] <- rnorm(10^5, mean = true_means[Time.step, i], sd = sqrt(true_variances[Time.step, i]))
}

df <- data.frame()
data_light <- read.csv(paste0("/storage/u1693998/results/results/light_lgssm_d", d, "N", Nparticles, "ID", ID), row.names = 1)
data_nsmc <- read.csv(paste0("/storage/u1693998/results/results/nsmc_lgssm_d", d, "N", Nparticles, "ID", ID), row.names = 1)

rse_dac_light <- (unname(matrix(data_light$V1, ncol = d, byrow = TRUE)) - true_means)^2/true_variances
rse_nsmc <- (unname(matrix(data_nsmc$V1, ncol = d, byrow = TRUE))  - true_means)^2/true_variances

last_t_nsmc <- unname(data.matrix(read.table(paste0("/storage/u1693998/results/t100_nsmc_lgssm_d", d, "N", Nparticles, "ID", ID), row.names = 1)))
last_t_light <- unname(data.matrix(read.table(paste0("/storage/u1693998/results/t100_light_lgssm_d", d, "N", Nparticles, "ID", ID), row.names = 1)))


res_dac_light_ks <- apply(rbind(last_t_light, marginals), ks_dist, N = Nparticles, MARGIN = 2)
res_dac_light_w1 <- apply(rbind(last_t_light, marginals), w1_dist, N = Nparticles, MARGIN = 2)

res_nsmc_ks <- apply(rbind(last_t_nsmc, marginals), ks_dist, N = Nparticles, MARGIN = 2)
res_nsmc_w1 <- apply(rbind(last_t_nsmc, marginals), w1_dist, N = Nparticles, MARGIN = 2)

df1 <- data.frame(cbind(t(rse_dac_light), res_dac_light_w1, res_dac_light_ks, rep(data_light$runtime[1], times = d)))
colnames(df1)[(Time.step+1):(Time.step+2)] <- c("V101", "V102")
df2 <- data.frame(rbind(df, cbind(t(rse_nsmc), res_nsmc_w1, res_nsmc_ks, rep(data_nsmc$runtime[1], times = d))))
colnames(df2)[(Time.step+1):(Time.step+2)] <- c("V101", "V102")

df <- rbind(df1, df2)
df$algo <- as.factor(rep(c("dac-light", "nsmc"), each = d))


filename <- paste0("run_dac/results/lgssm_d", d, "N", Nparticles, "ID", ID)
write.csv(x=df, file=filename)
