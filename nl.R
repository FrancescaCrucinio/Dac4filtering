set.seed(1234)
d <- 4
sigmaX <- 1
nu <- 10
delta <- 1
Time.step <- 20
y.error.prec <- matrix(0, nrow = d^2, ncol = d^2)
diag(y.error.prec) <- 1
diag(y.error.prec[-1, ]) <- 1/4
vertical_neighbours <- ((0:d^2) * (d^2 + 1) + d+1)
y.error.prec[vertical_neighbours[(vertical_neighbours <= d^4)]] <- 1/4
y.error.prec[upper.tri(y.error.prec)] = t(y.error.prec)[upper.tri(y.error.prec)]

nl_data <- nl_obs(d, sigmaX, nu, delta, y.error.prec, Time.step)
obs <- matrix(unlist(nl_data$y), ncol = d, nrow = d)

Nparticles <- 100
M <- 100
# initial state
history <- sqrt(sigmaX)*array(rnorm(Nparticles*d^2), dim = c(d, d, Nparticles))
xOld <- sqrt(sigmaX)*array(rnorm(Nparticles*M*d^2), dim = c(d, d, Nparticles, M))
xOld2 <- history
tic()
for (t in 1:Time.step){
  print(paste(t))
  res <- dac_nl_lightweight(history, obs, sigmaX, nu, covariance = TRUE)
  history <- res
}
toc()
tic()
for (t in 1:Time.step){
  print(paste(t))
  res_adaptive <- dac_nl_adaptive_lightweight(history, obs, sigmaX, nu, covariance = FALSE)
  history <- res_adaptive
}
toc()
tic()
for (t in 1:Time.step){
  print(paste(t))
  res_stpf <- stpf_nl(xOld, obs, nu, sigmaX)
  xOld <- res_stpf
}
toc()
tic()
for (t in 1:Time.step){
  print(paste(t))
  res_nsmc <- nsmc_nl(xOld2, obs, nu, sigmaX, M)
  xOld2 <- res_nsmc
}
toc()
apply(res, c(1,2), mean)
apply(res_stpf, c(1, 2), mean)

mean((apply(res, c(1,2), mean) - nl_data$x[, , Time.step+1])^2)
mean((apply(res_adaptive, c(1,2), mean) - nl_data$x[, , Time.step+1])^2)
mean((apply(res_nsmc, c(1,2), mean) - nl_data$x[, , Time.step+1])^2)
mean((apply(res_stpf, c(1, 2), mean) - nl_data$x[, , Time.step+1])^2)
