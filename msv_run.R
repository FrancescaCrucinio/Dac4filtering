# devtools::load_all("/storage/u1693998/Dac4filtering")
# ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
ID <- 1
# get data
y <- as.matrix(read.csv(file="data/synthetic_data_msv_nofactor_y"))
true_x <- read.csv(file="data/synthetic_data_msv_nofactor_x")
d <- dim(y)[2]
Time.step <- dim(y)[1]
Time.step <- 10
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
Nparticles <- 100

set.seed(1234*ID)
res_dac <- array(0, dim = c(Time.step, Nparticles, d))
tic()
res_dac[1, , ] <- marginal_dac_msv_first_step(y[1, ], SigmaV, Sigma0, Nparticles)
toc()
tic()
for (t in 2:Time.step) {
  res_dac[t, , ] <- marginal_dac_msv(res_dac[t-1, , ], y[t, ], y[t-1, ], SigmaU, SigmaV, SigmaUV, SigmaX, phi, adaptive = FALSE)
  print(paste(t))
}
runtime <- toc()

means.along <- function(a, i) {
  n <- length(dim(a))
  b <- aperm(a, c(seq_len(n)[-i], i))
  rowMeans(b, dims = n - 1)
}

dim <- 1
plot(1:Time.step, true_x[1:Time.step, dim], type = "l", col = "black", ylim = c(min(true_x[1:Time.step, dim], means.along(res_dac, 2)[, dim]),
                                                                             max(true_x[1:Time.step, dim], means.along(res_dac, 2)[, dim])))
lines(1:Time.step, means.along(res_dac, 2)[, dim], type = "l", col = "red")


means.along(res_dac, 2)[1,]
4.510652 1.879582 2.589651 3.163248
