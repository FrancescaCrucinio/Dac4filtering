### Linear Gaussian SSM -- comparison of dac, stpf, nsmc
set.seed(1234)
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

df_nsmc <- data.frame()
df_stpf <- data.frame()
df_dac <- data.frame()

m_dac <- matrix(0, nrow = Time.step, ncol = d)
v_dac <- matrix(0, nrow = Time.step, ncol = d)
runtime_dac <- rep(0, times = Time.step)
m_nsmc <- matrix(0, nrow = Time.step, ncol = d)
v_nsmc <- matrix(0, nrow = Time.step, ncol = d)
runtime_nsmc <- rep(0, times = Time.step)
m_stpf <- matrix(0, nrow = Time.step, ncol = d)
v_stpf <- matrix(0, nrow = Time.step, ncol = d)
runtime_stpf <- rep(0, times = Time.step)

# initial distribution
res_nsmc <- mvrnorm(n = Nparticles, mu0, Sigma0)
res_stpf <- array(mvrnorm(n = Nparticles*M, mu0, Sigma0), dim = c(Nparticles, M, d))
res_dac <- res_nsmc

for (t in 1:Time.step){
  y <- unname(data.matrix(read.csv(filename, row.names = 1, nrows = 1, skip = t-1)))[1:d]
  # dac (lightweight)
  tic()
  res_dac <- dac_lgssm_lightweight_crossover(res_dac, y, tau, lambda, sigmaY)
  rt <- toc()
  # save means and variances
  m_dac[t, ] <- colMeans(res_dac)
  v_dac[t, ] <- colVars(res_dac)
  runtime_dac[t] <- rt

  # nsmc
  tic()
  res_nsmc <- nsmc_lgssm(res_nsmc, y, tau, lambda, sigmaY, M)
  rt <- toc()
  m_nsmc[t, ] <- colMeans(res_nsmc)
  v_nsmc[t, ] <- colVars(res_nsmc)
  runtime_nsmc[t] <- rt

  # stpf
  tic()
  res_stpf <- stpf_lgssm(res_stpf, y, tau, lambda, sigmaY)
  rt <- toc()
  m_stpf[t, ] <- colMeans(res_stpf, dims = 2)
  v_stpf[t, ] <- colMeans(res_stpf^2, dims = 2) - m_stpf[t - (timeinterval-1)*10, ]^2
  runtime_stpf[t] <- rt
}


filename <- paste0("data/lgssm_tempering/dac_lgssm_d", d, "N", Nparticles, "ID", ID)
df_dac <- data.frame(rbind(m_dac, v_dac))
df_dac$type <- as.factor(rep(c("mean", "var"), each = Time.step))
df_dac$runtime <- rep(runtime_dac, times = 2)
write.table(x=df_dac, file=filename, append = TRUE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)
filename <- paste0("data/lgssm_tempering/nsmc_lgssm_d", d, "N", Nparticles, "ID", ID)
df_nsmc <- data.frame(rbind(m_nsmc, v_nsmc))
df_nsmc$type <- as.factor(rep(c("mean", "var"), each = Time.step))
df_nsmc$runtime <- rep(runtime_nsmc, each = 2)
write.table(x=df_nsmc, file=filename, append = TRUE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)
filename <- paste0("data/lgssm_tempering/stpf_lgssm_d", d, "N", Nparticles, "ID", ID)
df_stpf <- data.frame(rbind(m_stpf, v_stpf))
df_stpf$type <- as.factor(rep(c("mean", "var"), each = Time.step))
df_stpf$runtime <- rep(runtime_stpf, each = 2)
write.table(x=df_stpf, file=filename, append = TRUE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)

