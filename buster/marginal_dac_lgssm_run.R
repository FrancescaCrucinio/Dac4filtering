devtools::load_all("/storage/u1693998/Dac4filtering")

ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
d <- 256
timeinterval <- 1
# get time interval extrema
ti_begin <- 1 + (timeinterval - 1)*100
ti_end <- timeinterval*100
Time.step <- 100
set.seed(1234*ID)

filename <- paste0("/storage/u1693998/data/data_lgssm_d", d, "ID", ID)
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

if(timeinterval>1){
  filename <- paste0("/storage/u1693998/results/dac_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval-1)
  data_dac <- read.table(filename, row.names = 1)
  history <- data.matrix(data_dac)
} else{
  # initial value
  history <- mvrnorm(n = Nparticles, mu0, Sigma0)
}

m_dac <- matrix(0, nrow = 100, ncol = d)
v_dac <- matrix(0, nrow = 100, ncol = d)
runtime_dac <- rep(0, times = 100)
for (t in ti_begin:ti_end) {
  # dac (lightweight)
  tic()
  res_dac <- marginal_dac_lgssm_lightweight(history, y[t-ti_begin+1, ], tau, lambda, sigmaY)
  history <- res_dac
  rt <- toc()
  # save means and variances
  m_dac[t - (timeinterval-1)*100, ] <- colMeans(res_dac)
  v_dac[t - (timeinterval-1)*100, ] <- colVars(res_dac)
  runtime_dac[t - (timeinterval-1)*100] <- rt

  # save last time step to file
  if(t == ti_end){
    filename <- paste0("/storage/u1693998/results/dac_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
    write.table(res_dac, file = filename, append = FALSE, sep = " ", dec = ".",
                row.names = TRUE, col.names = TRUE)
  }
}

filename <- paste0("/storage/u1693998/results/marginal_dac_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", timeinterval)
df_dac <- data.frame(rbind(m_dac, v_dac))
df_dac$type <- as.factor(rep(c("mean", "var"), each = 100))
df_dac$runtime <- rep(runtime_dac, times = 2)
write.csv(x=df_dac, file=filename)
