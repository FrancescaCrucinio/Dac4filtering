set.seed(1234*5)
d <- 8
sigmaX <- 1
nu <- 10
tau <- -1/4
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
Nparticles <- 10
M <- 100
# initial state
history <- sqrt(sigmaX)*array(rnorm(Nparticles*d^2), dim = c(d, d, Nparticles))
history_dac_tempering <- history_dac
history_nsmc <- sqrt(sigmaX)*array(rnorm(Nparticles*d^2), dim = c(d, d, Nparticles))
history_stpf <- sqrt(sigmaX)*array(rnorm(Nparticles*M*d^2), dim = c(d, d, Nparticles, M))
tic()
for (t in 1:Time.step){
  res_dac <- marginal_dac_nl_lightweight(history_dac, y_cov[, , t], sigmaX, nu, adaptive = TRUE, covariance = TRUE, tau)
  history_dac <- res_dac
  print(paste(t))
}
toc()
# tic()
# obs_old <- matrix(0, ncol = d, nrow = d)
# for (t in 1:Time.step){
#   print(paste(t))
#   res_dac_tempering <- dac_nl_lightweight(history_dac_tempering, y_cov[, , t], sigmaX, nu, covariance = TRUE, tempering = FALSE,
#                                           obs_old = obs_old, tau = tau)
#   history_dac_tempering <- res_dac_tempering
#   obs_old <- y_cov[, , t]
# }
# toc()
tic()
for (t in 1:Time.step){
  print(paste(t))
  res_nsmc <- nsmc_nl(history_nsmc, y_cov[, , t], nu, sigmaX, M)
  history_nsmc <- res_nsmc
}
toc()
tic()
for (t in 1:Time.step){
  print(paste(t))
  res_stpf <- stpf_nl(history_stpf, y_cov[, , t], nu, sigmaX)
  history_stpf <- res_stpf
}
toc()
mean((apply(res_dac, c(1,2), mean) - nl_data$x[, , Time.step+1])^2, col = grey(seq(0, 1, length = 256)))
# mean((apply(res_dac_tempering, c(1,2), mean) - nl_data$x[, , Time.step+1])^2, col = grey(seq(0, 1, length = 256)))
mean((apply(res_nsmc, c(1,2), mean) - nl_data$x[, , Time.step+1])^2, col = grey(seq(0, 1, length = 256)))
mean((apply(res_stpf, c(1, 2), mean) - nl_data$x[, , Time.step+1])^2, col = grey(seq(0, 1, length = 256)))

