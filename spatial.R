set.seed(1234*1)
d <- 2
sigmaX <- 1
nu <- 10
tau_diag <- 1
tau <- -0.25
Time.step <- 1
spatial_data <- spatial_obs(d, sigmaX, nu, tau, tau_diag, Time.step)
y.error.prec <- matrix(c(1, -0.25, -0.25, 0, -0.25, 1, 0, -0.25, -0.25, 0, 1, -0.25, 0, -0.25, -0.25, 1), nrow = 4)
Nparticles <- 1000
history <- sqrt(sigmaX)*array(rnorm(Nparticles*d^2), dim = c(d, d, Nparticles))
res_dac <- history
tic()
for (t in 1:Time.step){
  res_dac <- marginal_dac_spatial(res_dac, spatial_data$y[, , t], sigmaX, nu, tau, tau_diag, adaptive = TRUE)
  # history <- res_dac
  print(paste(t))
}
toc()
res_bpf <- matrix(history, nrow = Nparticles, byrow = TRUE)
tic()
for (t in 1:Time.step){
  res_bpf <- spatial_bpf(res_bpf, sigmaX, nu, y.error.prec, c(spatial_data$y[, , t]), Nparticles)
  print(paste(t))
}
toc()
colMeans(res_bpf)
apply(res_dac, c(1,2), mean)
