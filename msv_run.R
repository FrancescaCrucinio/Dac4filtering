# devtools::load_all("/storage/u1693998/Dac4filtering")
# ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
ID <- 1
set.seed(1234*ID)
# get parameters
load("~/Documents/Dac4filtering/data/msv_parameters.RData")
# get data
y <- as.matrix(read.csv(file="data/synthetic_data_msv_y"))
true_x <- as.matrix(read.csv(file="data/synthetic_data_msv_x"))
# dimension
p <- 26
r <- 4
d <- p + r

# number of time steps
Time.step <- 1
# number of particles
Nparticles <- 1000

res_dac <- mvrnorm(n = Nparticles, mu, diag(Sigma^2))
tic()
for (t in 1:Time.step) {
  res_dac <- marginal_dac_msv(res_dac, y[t, ], mu, Phi, Lambda, Sigma, adaptive = FALSE)
  print(paste(t))
}
toc()


mean((colMeans(res_dac) - true_x[2, ])^2)
