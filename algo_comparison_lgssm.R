### Linear Gaussian SSM -- comparison of dac, stpf, nsmc
ID <- 1
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

Nparticles <- 100
M <- 100

# observations
filename <- paste0("data/data_lgssm_d", d, "ID", ID)

df_nsmc <- data.frame()
df_stpf <- data.frame()
df_dac <- data.frame()

# initial distribution
res_nsmc <- mvrnorm(n = Nparticles, mu0, Sigma0)
res_stpf <- array(mvrnorm(n = Nparticles*M, mu0, Sigma0), dim = c(Nparticles, M, d))
res_dac <- res_nsmc

for (t in 1:Time.step){
  y <- unname(data.matrix(read.csv(filename, row.names = 1, nrows = 1, skip = t-1)))[1:d]
  # dac (lightweight)
  tic()
  res_dac <- marginal_dac_lgssm_lightweight(res_dac, y, tau, lambda, sigmaY)
  rt <- toc()
  df_dac <- data.frame(rbind(df_dac, cbind(colMeans(res_dac), t, rt)))

  # nsmc
  tic()
  res_nsmc <- nsmc_lgssm(res_nsmc, y, tau, lambda, sigmaY, M)
  rt <- toc()
  df_nsmc <- data.frame(rbind(df_nsmc, cbind(colMeans(res_nsmc), t, rt)))

  # stpf
  tic()
  res_stpf <- stpf_lgssm(res_stpf, y, tau, lambda, sigmaY)
  rt <- toc()
  df_stpf <- data.frame(rbind(df_stpf, cbind(colMeans(res_stpf, dims = 2), t, rt)))
}

# save particles at last time step for W1 and KS
filename <- paste0("/data/results/t100_dac_lgssm_d", d, "N", Nparticles, "ID", ID)
write.table(res_dac, file = filename, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
filename <- paste0("/data/results/t100_nsmc_lgssm_d", d, "N", Nparticles, "ID", ID)
write.table(res_nsmc, file = filename, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
filename <- paste0("/data/results/t100_stpf_lgssm_d", d, "N", Nparticles, "ID", ID)
write.table(matrix(res_stpf, ncol = d, nrow = Nparticles*M), file = filename, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

# save means
write.csv(x=df_dac, file=paste0("/data/results/dac_lgssm_d", d, "N", Nparticles, "ID", ID))
write.csv(x=df_nsmc, file=paste0("/data/results/nsmc_lgssm_d", d, "N", Nparticles, "ID", ID))
write.csv(x=df_stpf, file=paste0("/data/results/stpf_lgssm_d", d, "N", Nparticles, "ID", ID))

