set.seed(1234)
y <- as.matrix(read.csv(file="data/synthetic_data_msv_y"))
true_x <- read.csv(file="data/synthetic_data_msv_x")
d <- dim(y)[2]
Time.step <- 100
# parameters
phi <- 0.1
Phi <- diag(rep(phi, times = d))
SigmaU <- matrix(0.3, nrow =d, ncol = d)
diag(SigmaU) <- 1
SigmaV <- matrix(0.2, nrow =d, ncol = d)
diag(SigmaV) <- 1
SigmaUV <- matrix(-0.1, nrow =d, ncol = d)
diag(SigmaUV) <- -0.2
Sigma0 <- SigmaU/(1-phi^2)
Sigma0_inv <- solve(Sigma0)
SigmaX <- SigmaU - SigmaUV %*% solve(SigmaV) %*% SigmaUV
# BPF
Nparticles <- 10000
res_bpf <- array(0, dim = c(Time.step, Nparticles, d))
tic()
lW <- rep(0, times = Nparticles)
for (i in 1:Nparticles) {
  res_bpf[1, i, ] <- mvrnorm(n = 1, rep(0, d), Sigma0)
  obs_covariance <- diag(sqrt(exp(res_bpf[1, i, ]))) %*% SigmaV %*% diag(sqrt(exp(res_bpf[1, i, ])))
  lW[i] <- -0.5*(y[1, ] %*% solve(obs_covariance) %*% y[1, ])
}
max.lW <- max(lW)
W <- exp(lW - max(lW))
W <- W/sum(W)
ancestors <- stratified_resample(W, Nparticles)
res_bpf[1, , ] <- res_bpf[1, ancestors, ]

for (t in 2:Time.step) {
  res_bpf[t, , ] <- msv_bpf(res_bpf[t-1, , ], y[t-1, ], y[t, ], SigmaV, SigmaUV, SigmaX, phi)
  print(paste(t))
}
runtime <- toc()

# STPF
Nparticles <- 50
M <- 20
res_stpf <- array(0, dim = c(Time.step, Nparticles, M, d))
tic()
res_stpf[1, , , ] <- stpf_msv_first_step(Nparticles, M, y[1, ], Sigma0, SigmaV)
for (t in 2:Time.step) {
  res_stpf[t, , , ] <- stpf_msv(res_stpf[t-1, , , ], y[t, ], y[t-1, ], phi, SigmaUV, SigmaV, SigmaX)
  print(paste(t))
}
runtime <- toc()

# DAC
Nparticles <- 100
res_dac <- array(0, dim = c(Time.step, Nparticles, d))
tic()
res_dac[1, , ] <- marginal_dac_msv_first_step(y[1, ], SigmaV, Sigma0, Nparticles, adaptive = FALSE)
toc()
tic()
for (t in 2:Time.step) {
  res_dac[t, , ] <- marginal_dac_msv(res_dac[t-1, , ], y[t, ], y[t-1, ], SigmaV, SigmaUV, SigmaX, phi, adaptive = FALSE)
  print(paste(t))
}
runtime <- toc()

# plot
dim <- 4
stpf_mean <- apply(res_stpf, c(1,4), mean)[, dim]
stpf_q1 <- apply(res_stpf, c(1,4), quantile, 0.25)[, dim]
stpf_q3 <- apply(res_stpf, c(1,4), quantile, 0.75)[, dim]
bpf_mean <- apply(res_bpf, c(1,3), mean)[, dim]
bpf_q1 <- apply(res_bpf, c(1,3), quantile, 0.25)[, dim]
bpf_q3 <- apply(res_bpf, c(1,3), quantile, 0.75)[, dim]
dac_mean <- apply(res_dac, c(1, 3), mean)[, dim]
dac_q1 <- apply(res_dac, c(1,3 ), quantile, 0.25)[, dim]
dac_q3 <- apply(res_dac, c(1, 3), quantile, 0.75)[, dim]

cbPalette <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
pdf("msv_comparison_d4.pdf", width = 8, height = 6)
plot(1:Time.step, true_x[1:Time.step, dim], type = "l", col = "black", ylim = c(-3, 4), lwd = 3, xlab = "Time", ylab = "Mean")
lines(1:Time.step, stpf_mean[1:Time.step], type = "l", col = cbPalette[1], lty = "dashed", lwd = 2)
# lines(1:Time.step, stpf_q1[1:Time.step], type = "l", col = cbPalette[1], lty = "dashed")
# lines(1:Time.step, stpf_q3[1:Time.step], type = "l", col = cbPalette[1], lty = "dashed")
lines(1:Time.step, bpf_mean[1:Time.step], type = "l", col = cbPalette[2], lty = "dashed", lwd = 2)
# lines(1:Time.step, bpf_q1[1:Time.step], type = "l", col = cbPalette[2], lty = "dashed")
# lines(1:Time.step, bpf_q3[1:Time.step], type = "l", col = cbPalette[2], lty = "dashed")
lines(1:Time.step, dac_mean, type = "l", col = cbPalette[5], lty = "dashed", lwd = 2)
# lines(1:Time.step, dac_q1, type = "l", col = cbPalette[5], lty = "dashed")
# lines(1:Time.step, dac_q3, type = "l", col = cbPalette[5], lty = "dashed")
legend("bottom", legend=c("truth", "stpf", "bpf", "dac"),
       col=c("black", cbPalette[1], cbPalette[2], cbPalette[5]), lty=c(1, 2, 2, 2), cex=1.1, ncol = 4, lwd = c(3, 2, 2, 2 ))
dev.off()

