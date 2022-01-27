# devtools::load_all("/storage/u1693998/Dac4filtering")

ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
ID <- 2
d <- 32
timeinterval <- 10
# get time interval extrema
ti_begin <- 1 + (timeinterval - 1)*10
ti_end <- timeinterval*10
Time.step <- 100
set.seed(1234*ID)

filename <- paste0("data/data_lgssm_d", d, "ID", ID)
data <- read.csv(filename, row.names = 1)

# parameters
mu0 <- rep(0, times = d)
Sigma0 <- diag(x = 1, d, d)
tau <- 1
lambda <- 1
sigmaY <- 0.5^2
y <- data.matrix(data[1:Time.step, 1:d], rownames.force = NA)

Nparticles <- 1000
M <- 2*d

history <- array(0, dim = c(Nparticles, d, 2))
if(timeinterval>1){
  filename <- paste0("data/dac_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval-1)
  data_dac <- read.csv(filename, row.names = 1)
  history[, , 2] <- data.matrix(data_dac[data_dac$timestep == ti_begin-2, 1:d])
  history[, , 1] <- data.matrix(data_dac[data_dac$timestep == ti_begin-1, 1:d])
  filename <- paste0("data/nsmc_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval-1)
  data_nsmc <- read.csv(filename, row.names = 1)
  res_nsmc <- data.matrix(data_nsmc[data_nsmc$timestep == ti_begin-1, 1:d])
  filename <- paste0("data/stpf_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval-1)
  data_stpf <- read.csv(filename, row.names = 1)
  res_stpf <- array(c(data.matrix(data_stpf[data_stpf$timestep == ti_begin-1, 1:d])), dim = c(Nparticles, M, d))
} else{
  # initial value
  res_nsmc <- mvrnorm(n = Nparticles, mu0, Sigma0)
  res_stpf <- array(mvrnorm(n = Nparticles*M, mu0, Sigma0), dim = c(Nparticles, M, d))
  history[, , 1] <- res_nsmc
}

df_dac <- data.frame()
df_nsmc <- data.frame()
df_stpf <- data.frame()

for (t in ti_begin:ti_end) {
  # dac (lightweight)
  tic()
  res_dac <- dac_lgssm_lightweight_crossover(history, y[t, ], tau, lambda, sigmaY, "adaptive")
  history[, , 2] <- history[, , 1]
  history[, , 1] <- res_dac
  rt <- toc()
  # save data
  timestep <- rep(t, times = Nparticles)
  runtime <- rep(rt$toc[[1]] - rt$tic[[1]], times = Nparticles)
  df_dac <- rbind(df_dac, cbind(res_dac, timestep, runtime))

  # nsmc
  tic()
  res_nsmc <- nsmc_lgssm(res_nsmc, y[t, ], tau, lambda, sigmaY, M)
  rt <- toc()
  runtime <- rep(rt$toc[[1]] - rt$tic[[1]], times = Nparticles)
  df_nsmc <- rbind(df_nsmc, cbind(res_nsmc, timestep, runtime))

  # stpf
  tic()
  res_stpf <- stpf_lgssm(res_stpf, y[t, ], tau, lambda, sigmaY)
  rt <- toc()
  runtime <- rep(rt$toc[[1]] - rt$tic[[1]], times = Nparticles*M)
  timestep <- rep(t, times = Nparticles*M)
  island <- rep(1:M, each = Nparticles)
  df_stpf <- rbind(df_stpf, cbind(matrix(res_stpf, nrow = Nparticles*M, ncol = d), timestep, runtime, island))
  print(paste(t))
}

filename <- paste0("data/dac_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
write.csv(x=df_dac, file=filename)
filename <- paste0("data/nsmc_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
write.csv(x=df_nsmc, file=filename)
filename <- paste0("data/stpf_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
write.csv(x=df_stpf, file=filename)


