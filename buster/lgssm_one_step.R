# devtools::load_all("/storage/u1693998/Dac4filtering")

ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
set.seed(1234*ID)
d <- 128
timeinterval <- 1
Time.step <- 10

filename <- paste0("/storage/u1693998/data/data_lgssm_d", d, "ID", ID)
# filename <- paste0("data/data_lgssm_d", d, "ID", ID)
data <- read.csv(filename, row.names = 1)
y <- unname(data.matrix(data[1:Time.step, 1:d], rownames.force = NA))
true_means <- unname(data.matrix(data[(Time.step+1):(2*Time.step), 1:d], rownames.force = NA))
# parameters
mu0 <- rep(0, times = d)
Sigma0 <- diag(x = 1, d, d)
tau <- 1
lambda <- 1
sigmaY <- 0.5^2

Nparticles <- 10000

if(timeinterval == 1){
  history <- array(0, dim = c(Nparticles, d, 2))
  x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
}

df <- data.frame()
# dac
if(timeinterval > 1){
  filename <- paste0("/storage/u1693998/results/nc_dac_lgssm_d", d, "N", Nparticles, "ID", ID, "step", timeinterval-1)
  history <- array(0, dim = c(Nparticles, d, 2))
  history[, , 2] <- matrix(0, ncol = d, nrow = Nparticles)
  history[, , 1] <- unname(data.matrix(read.table(filename, row.names = 1)))
}
tic()
res_dac <- dac_lgssm_lc_crossover(history, y[timeinterval, ], tau, lambda, sigmaY)
runtime <- toc()
m_dac <- (colMeans(res_dac) - true_means[timeinterval, ])^2
df <- data.frame(rbind(df, t(c(m_dac, runtime))))
filename <- paste0("/storage/u1693998/results/nc_dac_lgssm_d", d, "N", Nparticles, "ID", ID, "step", timeinterval)
write.table(res_dac, file = filename, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
# mix
#if(timeinterval > 1){
#  filename <- paste0("/storage/u1693998/results/mix_dac_lgssm_d", d, "N", Nparticles, "ID", ID, "step", timeinterval-1)
#  history <- array(0, dim = c(Nparticles, d, 2))
#  history[, , 2] <- matrix(0, ncol = d, nrow = Nparticles)
#  history[, , 1] <- unname(data.matrix(read.table(filename, row.names = 1)))
#}
#tic()
#res_dac_mix <- dac_lgssm_crossover(history, y[timeinterval, ], tau, lambda, sigmaY)
#runtime <- toc()
#m_dac_mix <- (colMeans(res_dac_mix) - true_means[timeinterval, ])^2
#df <- data.frame(rbind(df, t(c(m_dac_mix, runtime))))
#filename <- paste0("/storage/u1693998/results/mix_dac_lgssm_d", d, "N", Nparticles, "ID", ID, "step", timeinterval)
#write.table(res_dac_mix, file = filename, append = FALSE, sep = " ", dec = ".",
#            row.names = TRUE, col.names = TRUE)
# light
if(timeinterval > 1){
  filename <- paste0("/storage/u1693998/results/nc_light_dac_lgssm_d", d, "N", Nparticles, "ID", ID, "step", timeinterval-1)
  history <- array(0, dim = c(Nparticles, d, 2))
  history[, , 2] <- matrix(0, ncol = d, nrow = Nparticles)
  history[, , 1] <- unname(data.matrix(read.table(filename, row.names = 1)))
}
tic()
res_dac_light <- dac_lgssm_lightweight_crossover(history, y[timeinterval, ], tau, lambda, sigmaY)
runtime <- toc()
m_dac_light <- (colMeans(res_dac_light) - true_means[timeinterval, ])^2
df <- data.frame(rbind(df, t(c(m_dac_light, runtime))))
filename <- paste0("/storage/u1693998/results/nc_light_dac_lgssm_d", d, "N", Nparticles, "ID", ID, "step", timeinterval)
write.table(res_dac_light, file = filename, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
# adaptive light
if(timeinterval > 1){
  filename <- paste0("/storage/u1693998/results/nc_v2_dac_lgssm_d", d, "N", Nparticles, "ID", ID, "step", timeinterval-1)
  history <- array(0, dim = c(Nparticles, d, 2))
  history[, , 2] <- matrix(0, ncol = d, nrow = Nparticles)
  history[, , 1] <- unname(data.matrix(read.table(filename, row.names = 1)))
}
tic()
res_dac_v2 <- dac_lgssm_lightweight_crossover(history, y[timeinterval, ], tau, lambda, sigmaY, "adaptive")
runtime <- toc()
m_dac_v2 <- (colMeans(res_dac_v2) - true_means[timeinterval, ])^2
df <- data.frame(rbind(df, t(c(m_dac_v2, runtime))))
filename <- paste0("/storage/u1693998/results/nc_v2_dac_lgssm_d", d, "N", Nparticles, "ID", ID, "step", timeinterval)
write.table(res_dac_v2, file = filename, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

df$algo <- as.factor(c("dac", "light", "light_ada_temp"))
df$mutation <- as.factor(c(1,1,1))

filename <- paste0("run_dac/results/crossover_resampling_comparison_d", d, "N", Nparticles, "ID", ID, "step", timeinterval)
write.csv(x=df, file=filename)

