devtools::load_all("/storage/u1693998/Dac4filtering")
ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
set.seed(1234*ID)
# get parameters
load("/storage/u1693998/data/msv_parameters.RData")
# dimension
p <- 26
r <- 4
d <- p + r

# number of time steps
Time.step <- 122
# number of particles
Nparticles <- 100

m_dac <- matrix(0, nrow = Time.step, ncol = d)
exp_m_dac <- matrix(0, nrow = Time.step, ncol = d)
v_dac <- matrix(0, nrow = Time.step, ncol = d)
q1_dac <- matrix(0, nrow = Time.step, ncol = d)
q2_dac <- matrix(0, nrow = Time.step, ncol = d)
q3_dac <- matrix(0, nrow = Time.step, ncol = d)
min_dac <- matrix(0, nrow = Time.step, ncol = d)
max_dac <- matrix(0, nrow = Time.step, ncol = d)
res_dac <- mvrnorm(n = Nparticles, mu, diag(Sigma^2))
tic()
for (t in 1:Time.step) {
  y <- unname(data.matrix(read.csv(file="/storage/u1693998/data/synthetic_data_msv_y", nrows = 1, skip = t-1)))
  res_dac <- marginal_dac_msv(res_dac, y, mu, Phi, Lambda, Sigma, adaptive = FALSE)
  m_dac[t, ] <- colMeans(res_dac)
  exp_m_dac[t, ] <- colMeans(exp(res_dac))
  v_dac[t, ] <- colVars(res_dac)
  q1_dac[t, ] <- apply(res_dac, 2, quantile, probs = 0.25)
  q2_dac[t, ] <- apply(res_dac, 2, quantile, probs = 0.5)
  q3_dac[t, ] <- apply(res_dac, 2, quantile, probs = 0.75)
  min_dac[t, ] <- apply(res_dac, 2, min)
  max_dac[t, ] <- apply(res_dac, 2, max)
  print(paste(t))
}
runtime <- toc()
df <- data.frame(rbind(m_dac, v_dac, min_dac, q1_dac, q2_dac, q3_dac, mac_dac, exp_m_dac))
df$runtime <- runtime
df$type <- rep(c("mean", "var", "min", "q1", "q2", "q3", "max", "exp_m"), each = Time.step)
write.csv(x=df, file=paste0("/storage/u1693998/results/means_dac_msv_syn_d", d, "N", Nparticles, "ID", ID))
