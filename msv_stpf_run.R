# devtools::load_all("/storage/u1693998/Dac4filtering")
# ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
ID <- 1
# get data
y <- as.matrix(read.csv(file="data/synthetic_data_msv_nofactor_y"))
true_x <- read.csv(file="data/synthetic_data_msv_nofactor_x")
d <- dim(y)[2]
Time.step <- dim(y)[1]
Time.step <- 100
# parameters
phi <- 0.91
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
# number of particles
Nparticles <- 50
M <- 20

set.seed(1234*ID)
res_stpf <- array(0, dim = c(Time.step, Nparticles, M, d))
tic()
res_stpf[1, , , ] <- stpf_msv_first_step(Nparticles, M, y[1, ], Sigma0, SigmaV)
for (t in 2:Time.step) {
  res_stpf[t, , , ] <- stpf_msv(res_stpf[t-1, , , ], y[t, ], y[t-1, ], phi, SigmaUV, SigmaV, SigmaX)
  print(paste(t))
}
runtime <- toc()


dim <- 1
stpf_mean <- apply(res_stpf, c(1,4), mean)[, dim]
stpf_q1 <- apply(res_stpf, c(1,4), quantile, 0.25)[, dim]
stpf_q3 <- apply(res_stpf, c(1,4), quantile, 0.75)[, dim]
plot(1:Time.step, true_x[1:Time.step, dim], type = "l", col = "black", ylim = c(-10, 10))
lines(1:Time.step, stpf_mean, type = "l", col = "red")
lines(1:Time.step, stpf_q1, type = "l", col = "red", lty = "dashed")
lines(1:Time.step, stpf_q3, type = "l", col = "red", lty = "dashed")

