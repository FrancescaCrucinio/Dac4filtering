### Spatial model stats
# dimension
d <- 4
# parameters
sigmaX <- 1
nu <- 10
tau_diag <- 1
tau <- -0.25
# number of time steps
Time.step <- 10

Nparticles <- 5000
Nrep <- 50
df <- data.frame()
for (i in 1:50){
  dfnew <- read.csv(paste0("/storage/u1693998/results/means_adaptive_dac_spatial_d", d, "N", Nparticles, "ID", i, "step", 1), row.names = 1)
  colnames(dfnew) <- c("mean", "var", "first_q", "median", "third_q", "t", "runtime")
  dfnew$run <- i
  dfnew$dim <- 1:(d^2)
  df <- rbind(df, dfnew)
}
df$N <- Nparticles
write.csv(x=df, file=paste0("run_dac/results/new_stats_dac_spatial_tau", -tau, "d", d, "N", Nparticles))
