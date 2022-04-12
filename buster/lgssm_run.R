# devtools::load_all("/storage/u1693998/Dac4filtering")

ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
d <- 2048
timeinterval <- 1
# get time interval extrema
ti_begin <- 1 + (timeinterval - 1)*10
ti_end <- timeinterval*10
Time.step <- 100
set.seed(1234*ID)

filename <- paste0("/storage/u1693998/data/data_lgssm_d", d, "ID", ID)
# filename <- paste0("data/data_lgssm_d", d, "ID", ID)
data <- read.csv(filename, row.names = 1)

# parameters
mu0 <- rep(0, times = d)
Sigma0 <- diag(x = 1, d, d)
tau <- 1
lambda <- 1
sigmaY <- 0.5^2
y <- data.matrix(data[ti_begin:ti_end, 1:d], rownames.force = NA)
rm(data)
Nparticles <- 1000
M <- 100

history <- array(0, dim = c(Nparticles, d, 2))
if(timeinterval>1){
  filename <- paste0("/storage/u1693998/results/dac_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval-1)
  # filename <- paste0("data/dac_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval-1)
  data_dac <- read.table(filename, row.names = 1)
  history[, , 2] <- matrix(0, ncol = d, nrow = Nparticles)
  history[, , 1] <- data.matrix(data_dac)
  filename <- paste0("/storage/u1693998/results/nsmc_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval-1)
  # filename <- paste0("data/nsmc_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval-1)
  data_nsmc <- read.table(filename, row.names = 1)
  res_nsmc <- data.matrix(data_nsmc)
  filename <- paste0("/storage/u1693998/results/stpf_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval-1)
  # filename <- paste0("data/stpf_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval-1)
  #data_stpf <- read.table(filename, row.names = 1)
  #res_stpf <- array(c(data.matrix(data_stpf)), dim = c(Nparticles, M, d))
} else{
  # initial value
  res_nsmc <- mvrnorm(n = Nparticles, mu0, Sigma0)
  #res_stpf <- array(mvrnorm(n = Nparticles*M, mu0, Sigma0), dim = c(Nparticles, M, d))
  history[, , 1] <- res_nsmc
}

m_dac <- matrix(0, nrow = 10, ncol = d)
v_dac <- matrix(0, nrow = 10, ncol = d)
runtime_dac <- rep(0, times = 10)
m_nsmc <- matrix(0, nrow = 10, ncol = d)
v_nsmc <- matrix(0, nrow = 10, ncol = d)
runtime_nsmc <- rep(0, times = 10)
#m_stpf <- matrix(0, nrow = 10, ncol = d)
#v_stpf <- matrix(0, nrow = 10, ncol = d)
#runtime_stpf <- rep(0, times = 10)
for (t in ti_begin:ti_end) {
  # dac (lightweight)
  tic()
  res_dac <- dac_lgssm_lightweight_crossover(history, y[t-ti_begin+1, ], tau, lambda, sigmaY)
  history[, , 2] <- history[, , 1]
  history[, , 1] <- res_dac
  rt <- toc()
  # save means and variances
  m_dac[t - (timeinterval-1)*10, ] <- colMeans(res_dac)
  v_dac[t - (timeinterval-1)*10, ] <- colVars(res_dac)
  runtime_dac[t - (timeinterval-1)*10] <- rt

  # nsmc
  tic()
  res_nsmc <- nsmc_lgssm(res_nsmc, y[t-ti_begin+1, ], tau, lambda, sigmaY, M)
  rt <- toc()
  m_nsmc[t - (timeinterval-1)*10, ] <- colMeans(res_nsmc)
  v_nsmc[t - (timeinterval-1)*10, ] <- colVars(res_nsmc)
  runtime_nsmc[t - (timeinterval-1)*10] <- rt

  # stpf
  # tic()
  # res_stpf <- stpf_lgssm(res_stpf, y[t, ], tau, lambda, sigmaY)
  # rt <- toc()
  # m_stpf[t - (timeinterval-1)*10, ] <- colMeans(res_stpf, dims = 2)
  # v_stpf[t - (timeinterval-1)*10, ] <- colMeans(res_stpf^2, dims = 2) - m_stpf[t - (timeinterval-1)*10, ]^2
  # runtime_stpf[t - (timeinterval-1)*10] <- rt
  # save last time step to file
  if(t == ti_end){
    filename <- paste0("/storage/u1693998/results/dac_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
    # filename <- paste0("data/dac_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
    write.table(res_dac, file = filename, append = FALSE, sep = " ", dec = ".",
                row.names = TRUE, col.names = TRUE)
    filename <- paste0("/storage/u1693998/results/nsmc_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
    # filename <- paste0("data/nsmc_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
    write.table(res_nsmc, file = filename, append = FALSE, sep = " ", dec = ".",
                row.names = TRUE, col.names = TRUE)
    # filename <- paste0("/storage/u1693998/results/stpf_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
    # # filename <- paste0("data/stpf_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
    # write.table(matrix(res_stpf, nrow = Nparticles*M, ncol = d), file = filename, append = FALSE, sep = " ", dec = ".",
                # row.names = TRUE, col.names = TRUE)
  }
}

filename <- paste0("run_dac/results/dac_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
df_dac <- data.frame(rbind(m_dac, v_dac))
df_dac$type <- as.factor(rep(c("mean", "var"), each = 10))
df_dac$runtime <- rep(runtime_dac, times = 2)
write.table(x=df_dac, file=filename, append = TRUE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)
filename <- paste0("run_dac/results/nsmc_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
df_nsmc <- data.frame(rbind(m_nsmc, v_nsmc))
df_nsmc$type <- as.factor(rep(c("mean", "var"), each = 10))
df_nsmc$runtime <- rep(runtime_nsmc, each = 2)
write.table(x=df_nsmc, file=filename, append = TRUE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)
# filename <- paste0("/storage/u1693998/results/stpf_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
# # filename <- paste0("data/results_stpf_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
# df_stpf <- data.frame(rbind(m_stpf, v_stpf))
# df_stpf$type <- as.factor(rep(c("mean", "var"), each = 10))
# df_stpf$runtime <- rep(runtime_stpf, each = 2)
# write.table(x=df_stpf, file=filename, append = TRUE, sep = " ", dec = ".",
#             row.names = FALSE, col.names = FALSE)


