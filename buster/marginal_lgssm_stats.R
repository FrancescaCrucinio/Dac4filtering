### Linear Gaussian SSM -- comparison of dac, stpf, nsmc
devtools::load_all("/storage/u1693998/Dac4filtering")

ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))

set.seed(1234*ID)
# dimension
d <- 32
# initial state
mu0 <- rep(0, times = d)
Sigma0 <- diag(x = 1, d, d)
# parameters
tau <- 1
lambda <- 1
sigmaY <- 0.5^2
Time.step <- 100
timeinterval <- 1

Nparticles <- 100
M <- 100

# get data
filename <- paste0("/storage/u1693998/data/data_lgssm_d", d, "ID", ID)
true_means <- unname(data.matrix(read.csv(filename, row.names = 1, nrows=Time.step, skip=Time.step)))[, 1:d]
true_variances <- unname(data.matrix(read.csv(filename, row.names = 1, nrows=Time.step, skip=2*Time.step)))[, 1:d]

marginals <- matrix(0, nrow = 10^5, ncol = d)
for(i in 1:d){
  marginals[, i] <- rnorm(10^5, mean = true_means[Time.step, i], sd = sqrt(true_variances[Time.step, i]))
}

df <- data.frame()
data_dac <- read.csv(paste0("marginal_dac_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval))

rse_dac <- (unname(data.matrix(data_dac[1:Time.step, 2:(d+1)])) - true_means)^2/true_variances

last_t_dac <- unname(data.matrix(read.table(paste0("/storage/u1693998/results/dac_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval), row.names = 1)))

res_dac_ks <- apply(rbind(last_t_dac, marginals), ks_dist, N = Nparticles, MARGIN = 2)
res_dac_w1 <- apply(rbind(last_t_dac, marginals), w1_dist, N = Nparticles, MARGIN = 2)

df <- data.frame(cbind(t(rse_dac), res_dac_w1, res_dac_ks, rep(sum(data_dac$runtime), times = d)))
colnames(df)[(Time.step+1):(Time.step+2)] <- c("V101", "V102")

df$algo <- as.factor(rep("dac-light-marginal", each = d))


filename <- paste0("run_dac/results/marginal_lgssm_d", d, "N", Nparticles, "ID", ID)
write.csv(x=df, file=filename)
