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
Nparticles <- 50
M <- 20
m_stpf <- matrix(0, nrow = Time.step, ncol = d)
v_stpf <- matrix(0, nrow = Time.step, ncol = d)
q1_stpf <- matrix(0, nrow = Time.step, ncol = d)
q2_stpf <- matrix(0, nrow = Time.step, ncol = d)
q3_stpf <- matrix(0, nrow = Time.step, ncol = d)
min_stpf <- matrix(0, nrow = Time.step, ncol = d)
max_stpf <- matrix(0, nrow = Time.step, ncol = d)
runtime <- rep(0, times = Time.step)
tic()
t <- 1
y <- c(data.matrix(read.csv(file="/storage/u1693998/data/synthetic_data_msv_y", nrows = 1, skip = t-1, header = FALSE)))
res_stpf <-  stpf_msv_first_step(Nparticles, M, y, Sigma0, SigmaV)
m_stpf[t, ] <- apply(res_stpf, 3, mean)
v_stpf[t, ] <- apply(res_stpf, 3, function(x) {var(c(x))})
q1_stpf[t, ] <- apply(res_stpf, 3, quantile, probs = 0.25)
q2_stpf[t, ] <- apply(res_stpf, 3, quantile, probs = 0.5)
q3_stpf[t, ] <- apply(res_stpf, 3, quantile, probs = 0.75)
min_stpf[t, ] <- apply(res_stpf, 3, min)
max_stpf[t, ] <- apply(res_stpf, 3, max)
y_past <- y
print(paste(t))
runtime[t] <- toc()
for (t in 2:Time.step) {
  tic()
  y <- c(data.matrix(read.csv(file="/storage/u1693998/data/synthetic_data_msv_y", nrows = 1, skip = t-1, header = FALSE)))
  res_stpf <- stpf_msv(res_stpf, y, y_past, phi, SigmaUV, SigmaV, SigmaX)
  m_stpf[t, ] <- apply(res_stpf, 3, mean)
  v_stpf[t, ] <- apply(res_stpf, 3, function(x) {var(c(x))})
  q1_stpf[t, ] <- apply(res_stpf, 3, quantile, probs = 0.25)
  q2_stpf[t, ] <- apply(res_stpf, 3, quantile, probs = 0.5)
  q3_stpf[t, ] <- apply(res_stpf, 3, quantile, probs = 0.75)
  min_stpf[t, ] <- apply(res_stpf, 3, min)
  max_stpf[t, ] <- apply(res_stpf, 3, max)
  y_past <- y
  print(paste(t))
  runtime[t] <- toc()
}
df <- data.frame(rbind(m_stpf, v_stpf, min_stpf, q1_stpf, q2_stpf, q3_stpf, max_stpf))
df$runtime <- sum(runtime)
df$type <- rep(c("mean", "var", "min", "q1", "q2", "q3", "max"), each = Time.step)
write.csv(x=df, file=paste0("/storage/u1693998/results/means_stpf_msv_syn_d", d, "N", Nparticles, "ID", ID))
