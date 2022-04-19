devtools::load_all("/storage/u1693998/Dac4filtering")
### Linear Gaussian SSM -- comparison of dac, stpf, nsmc
ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
set.seed(1234*ID)
# dimension
d <- 32
# parameters
sigmaX <- 1
nu <- 10
tau <- 1/4
delta <- 1
Time.step <- 100

Nparticles <- 1000
M <- 100
df_nsmc <- data.frame()
df_stpf <- data.frame()

# initial distribution
res_nsmc <- sqrt(sigmaX)*array(rnorm(Nparticles*d^2), dim = c(d, d, Nparticles))
res_stpf <- sqrt(sigmaX)*array(rnorm(Nparticles*M*d^2), dim = c(d, d, Nparticles, M))

# nsmc
tic()
for (t in 1:Time.step){
  y <- unname(data.matrix(read.csv(paste0("/storage/u1693998/data/data_iid_nl_tau_", tau, "d", d, "ID", ID),
                                   row.names = 1, nrows=d, skip=(t-1)*d)))
  res_nsmc <- nsmc_nl(res_nsmc, y, nu, sigmaX, M)
  df_nsmc <- rbind(df_nsmc, cbind(apply(res_nsmc, c(1,2), mean), rep(t, times = d)))
}
runtime <- toc()
df_nsmc$runtime <- runtime
# stpf
tic()
for (t in 1:Time.step){
  y <- unname(data.matrix(read.csv(paste0("/storage/u1693998/data/data_iid_nl_tau_", tau, "d", d, "ID", ID),
                                   row.names = 1, nrows=d, skip=(t-1)*d)))
  res_stpf <- stpf_nl(res_stpf, y, nu, sigmaX)
  df_stpf <- rbind(df_stpf, cbind(apply(res_stpf, c(1,2), mean), rep(t, times = d)))
}
runtime <- toc()
df_stpf$runtime <- runtime

write.csv(x=df_nsmc, file=paste0("/storage/u1693998/results/results/nsmc_nl_iid_d", d, "N", Nparticles, "ID", ID))
write.csv(x=df_stpf, file=paste0("/storage/u1693998/results/results/stpf_nl_iid_d", d, "N", Nparticles, "ID", ID))
