### Linear Gaussian SSM -- comparison of dac, stpf, nsmc
ID <- 1
set.seed(1234*ID)
# dimension
d <- 8
# parameters
sigmaX <- 1
nu <- 10
tau_diag <- 1
tau <- -0.25
y.error.prec <- matrix(0, nrow = d^2, ncol = d^2)
diag(y.error.prec) <- tau_diag
vertical_neighbours <- ((0:d^2) * (d^2 + 1) + d+1)
horizontal_neighbours <- (0:(d^2-1)) * (d^2+1) + 2
y.error.prec[vertical_neighbours[(vertical_neighbours <= d^4)]] <- tau
y.error.prec[horizontal_neighbours[mod(horizontal_neighbours, 2) == 0]] <- tau
y.error.prec[upper.tri(y.error.prec)] <- t(y.error.prec)[upper.tri(y.error.prec)]
# number of time steps
Time.step <- 10
# data
timeinterval <- 1
ti_begin <- 1 + (timeinterval - 1)*10
ti_end <- timeinterval*10
Nparticles <- 100
# initial state
res_bpf <- matrix(rnorm(Nparticles*d^2, sd = sqrt(sigmaX)), nrow = Nparticles, ncol = d^2)
df_bpf <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df_bpf) <- c("mean", "var", "first_q", "second_q", "third_q", "time")

tic()
for (t in ti_begin:ti_end){
  print(paste(t))
  y <- unname(data.matrix(read.csv(paste0("/data/data_spatial_tau_", -tau, "d", d, "ID", 1),
                                   row.names = 1, nrows=d, skip=(t-1)*d)))
  res_bpf <- spatial_bpf(res_bpf, sigmaX, y.error.prec, c(y), Nparticles)
  df_bpf <- rbind(df_bpf, cbind(colMeans(res_bpf), colVars(res_bpf), apply(res_bpf, 2, quantile, probs = 0.25),
                                apply(res_bpf, 2, quantile, probs = 0.5), apply(res_bpf, 2, quantile, probs = 0.75),
                                rep(t, times = d^2)))
  print(paste(t))
}
runtime <- toc()
df_bpf$runtime <- runtime
write.csv(x=df_bpf, file=paste0("/data/results/means_bpf_spatial_d", d, "N", Nparticles, "ID", ID, "step", timeinterval))

