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
df_light <- data.frame()
df_nsmc <- data.frame()

res_dac_light <- mvrnorm(n = Nparticles, mu0, Sigma0)
res_nsmc <- res_dac_light
tic()
for (t in 1:Time.step) {
  y <- unname(data.matrix(read.csv(paste0("/storage/u1693998/data/data_lgssm_d", d, "ID", ID), row.names = 1, nrows=1, skip=t-1)))[1:d]
  res_dac_light <- dac_lgssm_lightweight_crossover(res_dac_light, y, tau, lambda, sigmaY)

  df_light <- data.frame(rbind(df_light, cbind(colMeans(res_dac_light), t)))
}
runtime <- toc()
df_light$runtime <- runtime
# nsmc
tic()
for (t in 1:Time.step) {
  y <- unname(data.matrix(read.csv(paste0("/storage/u1693998/data/data_lgssm_d", d, "ID", ID), row.names = 1, nrows=1, skip=t-1)))[1:d]
  res_nsmc <- nsmc_lgssm(res_nsmc, y, tau, lambda, sigmaY, M)
  df_nsmc <- data.frame(rbind(df_nsmc, cbind(colMeans(res_nsmc), t)))
}
runtime <- toc()
df_nsmc$runtime <- runtime

filename <- paste0("/storage/u1693998/results/t100_light_lgssm_d", d, "N", Nparticles, "ID", ID)
write.table(res_dac_light, file = filename, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
filename <- paste0("/storage/u1693998/results/t100_nsmc_lgssm_d", d, "N", Nparticles, "ID", ID)
write.table(res_nsmc, file = filename, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

write.csv(x=df_nsmc, file=paste0("run_dac/results/nsmc_lgssm_d", d, "N", Nparticles, "ID", ID))
write.csv(x=df_light, file=paste0("run_dac/results/light_lgssm_d", d, "N", Nparticles, "ID", ID))
