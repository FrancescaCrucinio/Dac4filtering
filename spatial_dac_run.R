### Spatial model
ID <- 1
set.seed(1234*ID)
# dimension
d <- 32
# parameters
sigmaX <- 1
nu <- 10
tau_diag <- 1
tau <- -0.25
# number of time steps
Time.step <- 10
timeinterval <- 1
ti_begin <- 1 + (timeinterval - 1)*10
ti_end <- timeinterval*10
Nparticles <- 5000
df_dac <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df_dac) <- c("mean", "var", "first_q", "second_q", "third_q", "time")
# initial state
res_dac <- sqrt(sigmaX)*array(rnorm(Nparticles*d^2), dim = c(d, d, Nparticles))
# observations
filename <- paste0("data/data_spatial_tau_", -tau, "d", d, "ID", 1)

tic()
for (t in ti_begin:ti_end){
  print(paste(t))
  y <- unname(data.matrix(read.csv(filename, row.names = 1, nrows=d, skip=(t-1)*d)))
  res_dac <- marginal_dac_spatial(res_dac, y, sigmaX, nu, tau, adaptive = TRUE)
  df_dac <- rbind(df_dac, cbind(c(apply(res_dac, c(1,2), mean)), c(apply(res_dac, c(1,2), var)), c(apply(res_dac, c(1,2), quantile, probs = 0.25)),
	c(apply(res_dac, c(1,2), quantile, probs = 0.5)), c(apply(res_dac, c(1,2), quantile, probs = 0.75)), rep(t, times = d^2)))
  print(paste(t))
}
runtime <- toc()
df_dac$runtime <- runtime
filename <- paste0("/data/results/dac_spatial_d", d, "N", Nparticles, "ID", ID, "step", timeinterval)
write.table(data.frame(matrix(res_dac, ncol = Nparticles)), file = filename, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

write.csv(x=df_dac, file=paste0("/data/results/means_dac_spatial_d", d, "N", Nparticles, "ID", ID, "step", timeinterval))
