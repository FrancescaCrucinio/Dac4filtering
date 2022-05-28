set.seed(1234*1)
d <- 16
sigmaX <- 1
nu <- 10
tau_diag <- 1
tau <- -0.25
Time.step <- 3
spatial_data <- spatial_obs(d, sigmaX, nu, tau, tau_diag, Time.step)
y <- spatial_data$y
Nparticles <- 1000
# initial state
history <- sqrt(sigmaX)*array(rnorm(Nparticles*d^2), dim = c(d, d, Nparticles))
tic()
for (t in 1:Time.step){
  res_dac <- marginal_dac_spatial(history, y[, , t], sigmaX, nu, tau, adaptive = TRUE)
  history <- res_dac
  print(paste(t))
}
toc()
boxplot(c((apply(res_dac, c(1,2), mean) - spatial_data$x[, , Time.step+1])^2))
image((apply(res_dac, c(1,2), mean) - spatial_data$x[, , Time.step+1])^2)
mean((apply(res_dac, c(1,2), mean) - spatial_data$x[, , Time.step+1])^2)
