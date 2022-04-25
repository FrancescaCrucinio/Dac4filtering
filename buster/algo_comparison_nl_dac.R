devtools::load_all("/storage/u1693998/Dac4filtering")
### Linear Gaussian SSM -- comparison of dac, stpf, nsmc
ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
set.seed(1234*ID)
# dimension
d <- 32
# parameters
sigmaX <- 1
nu <- 10
tau <- -1/4
delta <- 1
Time.step <- 100

Nparticles <- 100
df_dac <- data.frame()

# initial distribution
res_dac <- sqrt(sigmaX)*array(rnorm(Nparticles*d^2), dim = c(d, d, Nparticles))

# dac
obs_old <- matrix(0, ncol = d, nrow = d)
tic()
for (t in 1:Time.step){
  y <- unname(data.matrix(read.csv(paste0("/storage/u1693998/data/data_cov_nl_tau_", -tau, "d", d, "ID", ID),
                                   row.names = 1, nrows=d, skip=(t-1)*d)))
  res_dac <- dac_nl_lightweight(res_dac, y, sigmaX, nu, covariance = TRUE, tempering = FALSE, obs_old = obs_old, tau = tau)
  df_dac <- rbind(df_dac, cbind(apply(res_dac, c(1,2), mean), rep(t, times = d)))
  obs_old <- y
}
runtime <- toc()
df_dac$runtime <- runtime

write.csv(x=df_dac, file=paste0("/storage/u1693998/results/results/dac_nl_cov_d", d, "N", Nparticles, "ID", ID))
