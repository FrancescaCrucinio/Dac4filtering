devtools::load_all("/storage/u1693998/Dac4filtering")
### Linear Gaussian SSM -- comparison of dac, stpf, nsmc
ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
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
ti_begin <- 1 + (timeinterval - 1)*1
ti_end <- timeinterval*1

Nparticles <- 10000
df_dac <- data.frame()

if(timeinterval == 1){
  res_dac <- sqrt(sigmaX)*array(rnorm(Nparticles*d^2), dim = c(d, d, Nparticles))
} else {
  filename <- paste0("/storage/u1693998/results/dac_spatial_d", d, "N", Nparticles, "ID", ID, "step", timeinterval-1)
  res_dac <- unname(data.matrix(read.table(filename, row.names = 1)))
  dim(res_dac) <- c(d, d, Nparticles)
}

tic()
for (t in ti_begin:ti_end){
  y <- unname(data.matrix(read.csv(paste0("data/data_spatial_tau_", -tau, "d", d, "ID", ID),
                                   row.names = 1, nrows=d, skip=(t-1)*d)))
  res_dac <- marginal_dac_spatial(res_dac, y, sigmaX, nu, tau, adaptive = TRUE)
  df_dac <- rbind(df_dac, cbind(apply(res_dac, c(1,2), mean), rep(t, times = d)))
  print(paste(t))
}
runtime <- toc()
df_dac$runtime <- runtime
filename <- paste0("/storage/u1693998/results/dac_spatial_d", d, "N", Nparticles, "ID", ID, "step", timeinterval)
write.table(data.frame(matrix(res_dac, ncol = Nparticles)), file = filename, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

write.csv(x=df_dac, file=paste0("/storage/u1693998/results/means_dac_spatial_d", d, "N", Nparticles, "ID", ID, "step", timeinterval))
