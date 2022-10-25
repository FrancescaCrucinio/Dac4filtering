# devtools::load_all("/storage/u1693998/Dac4filtering")
# ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
ID <- 1
set.seed(1234*ID)
# get parameters
load("~/Documents/Dac4filtering/data/msv_parameters.RData")
# get data
y <- as.matrix(read.csv(file="data/synthetic_data_msv_y"))
true_x <- read.csv(file="data/synthetic_data_msv_x")
# dimension
p <- 26
r <- 4
d <- p + r

# number of time steps
Time.step <- 10
# number of particles
Nparticles <- 500

res_dac <- array(0, dim = c(Time.step+1, Nparticles, d))
res_dac[1, , ] <- mvrnorm(n = Nparticles, mu, diag(Sigma^2))
tic()
for (t in 1:Time.step) {
  res_dac[t+1, , ] <- marginal_dac_msv(res_dac[t, , ], y[t, ], mu, Phi, Lambda, Sigma, adaptive = FALSE)
  print(paste(t, t, t ,t, t, t))
}
toc()

means.along <- function(a, i) {
  n <- length(dim(a))
  b <- aperm(a, c(seq_len(n)[-i], i))
  rowMeans(b, dims = n - 1)
}
(means.along(res_dac, 2) - true_x[1:11, ])^2/true_x[1:11, ]^2

plot(zoo(true_x[1:11, 2]))
lines((means.along(res_dac, 2))[1:11, 2], col = "red", type = "l")

