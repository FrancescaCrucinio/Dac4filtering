devtools::load_all("/storage/u1693998/Dac4filtering")
### Linear Gaussian SSM -- comparison of dac, stpf, nsmc
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

filename <- paste0("/storage/u1693998/data/data_lgssm_d", d, "ID", ID)
# filename <- paste0("data/data_lgssm_d", d, "ID", ID)
data <- read.csv(filename, row.names = 1)
y <- unname(data.matrix(data[1:Time.step, 1:d], rownames.force = NA))
true_means <- unname(data.matrix(data[data$type == "mean", 1:d]))
true_variances <- unname(data.matrix(data[data$type == "var", 1:d]))
# samples from marginals at last time step
marginals <- matrix(0, nrow = 10^5, ncol = d)
for(i in 1:d){
  marginals[, i] <- rnorm(10^5, mean = true_means[Time.step, i], sd = sqrt(true_variances[Time.step, i]))
}


Nparticles <- 100
M <- 100
df <- data.frame()

x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
# dac (lightweight adaptive)
tic()
res_dac <- dac_time_lgssm_crossover(tau, lambda, sigmaY, Nparticles, x0, y, method = "adaptive", marginals = marginals)
runtime <- toc()
rse_dac <- (res_dac$m - true_means)^2/true_variances
df <- data.frame(rbind(df, cbind(t(rse_dac), res_dac$w1, res_dac$ks, rep(runtime, times = d))))
# dac (lightweight)
tic()
res_dac_light <- dac_time_lgssm_crossover(tau, lambda, sigmaY, Nparticles, x0, y, method = "light", marginals = marginals)
runtime <- toc()
rse_dac_light <- (res_dac_light$m - true_means)^2/true_variances
df <- data.frame(rbind(df, cbind(t(rse_dac_light), res_dac_light$w1, res_dac_light$ks, rep(runtime, times = d))))
# nsmc
tic()
res_nsmc <- nsmc_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, M = M, marginals = marginals)
runtime <- toc()
rse_nsmc <- (res_nsmc$m - true_means)^2/true_variances
df <- data.frame(rbind(df, cbind(t(rse_nsmc), res_nsmc$w1, res_nsmc$ks, rep(runtime, times = d))))
# stpf
x0 <- array(mvrnorm(n = Nparticles*M, mu0, Sigma0), dim = c(Nparticles, M, d))
tic()
res_stpf <- stpf_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, M = M, marginals = marginals)
runtime <- toc()
rse_stpf <- (res_stpf$m - true_means)^2/true_variances
df <- data.frame(rbind(df, cbind(t(rse_stpf), res_stpf$w1, res_stpf$ks, rep(runtime, times = d))))

df$algo <- as.factor(rep(c("dac-ada", "dac-light", "nsmc", "stpf"), each = d))


filename <- paste0("run_dac/results/lgssm_d", d, "N", Nparticles, "ID", ID)
write.csv(x=df, file=filename)
