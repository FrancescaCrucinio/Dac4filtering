set.seed(1234)
d <- 4
sigmaX <- 1
nu <- 10
tau <- 1/4
delta <- 1
Time.step <- 1
y.error.prec <- matrix(0, nrow = d^2, ncol = d^2)
diag(y.error.prec) <- 1
diag(y.error.prec[-1, ]) <- tau
vertical_neighbours <- ((0:d^2) * (d^2 + 1) + d+1)
y.error.prec[vertical_neighbours[(vertical_neighbours <= d^4)]] <- tau
y.error.prec[upper.tri(y.error.prec)] <- t(y.error.prec)[upper.tri(y.error.prec)]
y.error.var <- inv(y.error.prec)
y.error.var[upper.tri(y.error.var)] <- t(y.error.var)[upper.tri(y.error.var)]

nl_data <- nl_obs(d, sigmaX, nu, delta, y.error.var, Time.step)
y <- nl_data$yiid
y_cov <- nl_data$y
Nparticles <- 100
M <- 100
# initial state
history <- sqrt(sigmaX)*array(rnorm(Nparticles*d^2), dim = c(d, d, Nparticles))
xOld <- sqrt(sigmaX)*array(rnorm(Nparticles*M*d^2), dim = c(d, d, Nparticles, M))
xOld2 <- history
tic()
for (t in 1:Time.step){
  print(paste(t))
  res_dac <- dac_nl_lightweight(history, y_cov[, , t], sigmaX, nu, covariance = FALSE, tempering = FALSE)
  # history <- res_dac
}
toc()
tic()
for (t in 1:Time.step){
  print(paste(t))
  res_dac_tempering <- dac_nl_lightweight(history, y_cov[, , t], sigmaX, nu, M = NULL, covariance = TRUE,
                                          obs_old = matrix(0, nrow = d, ncol = d), tau = tau)
  # history <- res_dac_tempering
}
toc()
tic()
for (t in 1:Time.step){
  print(paste(t))
  res_stpf <- stpf_nl(xOld, y_cov[, , t], nu, sigmaX)
  xOld <- res_stpf
}
toc()
tic()
for (t in 1:Time.step){
  print(paste(t))
  res_nsmc <- nsmc_nl(xOld2, y_cov[, , t], nu, sigmaX, M)
  xOld2 <- res_nsmc
}
toc()

mean((apply(res_dac, c(1,2), mean) - nl_data$x[, , Time.step+1])^2)
mean((apply(res_dac_tempering, c(1,2), mean) - nl_data$x[, , Time.step+1])^2)
mean((apply(res_nsmc, c(1,2), mean) - nl_data$x[, , Time.step+1])^2)
mean((apply(res_stpf, c(1, 2), mean) - nl_data$x[, , Time.step+1])^2)
