devtools::load_all("/storage/u1693998/Dac4filtering")
ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
set.seed(1234*ID)
# get parameters
d <- 128
# model parameters
phi <- 0.1
SigmaU <- matrix(0.3, nrow =d, ncol = d)
diag(SigmaU) <- 1
SigmaV <- matrix(0.2, nrow =d, ncol = d)
diag(SigmaV) <- 1
SigmaUV <- matrix(-0.1, nrow =d, ncol = d)
diag(SigmaUV) <- -0.2
Sigma0 <- SigmaU/(1-phi^2)
SigmaX <- SigmaU - SigmaUV %*% solve(SigmaV) %*% SigmaUV

# number of time steps
Time.step <- 100
# number of particles
Nparticles <- 10

m_dac <- matrix(0, nrow = Time.step, ncol = d)
v_dac <- matrix(0, nrow = Time.step, ncol = d)
q1_dac <- matrix(0, nrow = Time.step, ncol = d)
q2_dac <- matrix(0, nrow = Time.step, ncol = d)
q3_dac <- matrix(0, nrow = Time.step, ncol = d)
min_dac <- matrix(0, nrow = Time.step, ncol = d)
max_dac <- matrix(0, nrow = Time.step, ncol = d)
tic()
t <- 1
y <- c(data.matrix(read.csv(file="/storage/u1693998/data/synthetic_data_msv_y", nrows = 1, skip = t-1, header = FALSE)))
res_dac <- marginal_dac_msv_first_step(y, SigmaV, Sigma0, Nparticles)
m_dac[t, ] <- colMeans(res_dac)
v_dac[t, ] <- colVars(res_dac)
q1_dac[t, ] <- apply(res_dac, 2, quantile, probs = 0.25)
q2_dac[t, ] <- apply(res_dac, 2, quantile, probs = 0.5)
q3_dac[t, ] <- apply(res_dac, 2, quantile, probs = 0.75)
min_dac[t, ] <- apply(res_dac, 2, min)
max_dac[t, ] <- apply(res_dac, 2, max)
y_past <- y
print(paste(t))
for (t in 2:Time.step) {
  y <- c(data.matrix(read.csv(file="/storage/u1693998/data/synthetic_data_msv_y", nrows = 1, skip = t-1, header = FALSE)))
  res_dac <- marginal_dac_msv(res_dac, y, y_past, SigmaV, SigmaUV, SigmaX, phi, adaptive = FALSE)
  m_dac[t, ] <- colMeans(res_dac)
  v_dac[t, ] <- colVars(res_dac)
  q1_dac[t, ] <- apply(res_dac, 2, quantile, probs = 0.25)
  q2_dac[t, ] <- apply(res_dac, 2, quantile, probs = 0.5)
  q3_dac[t, ] <- apply(res_dac, 2, quantile, probs = 0.75)
  min_dac[t, ] <- apply(res_dac, 2, min)
  max_dac[t, ] <- apply(res_dac, 2, max)
  y_past <- y
  print(paste(t))
}
runtime <- toc()
df <- data.frame(rbind(m_dac, v_dac, min_dac, q1_dac, q2_dac, q3_dac, mac_dac, exp_m_dac))
df$runtime <- runtime
df$type <- rep(c("mean", "var", "min", "q1", "q2", "q3", "max"), each = Time.step)
write.csv(x=df, file=paste0("/storage/u1693998/results/means_dac_msv_syn_d", d, "N", Nparticles, "ID", ID))
