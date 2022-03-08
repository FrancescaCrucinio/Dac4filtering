set.seed(1234)
d <- 4
sigmaX <- 1
nu <- 10
delta <- 1
Time.step <- 10
y.error.prec <- matrix(0, nrow = d^2, ncol = d^2)
diag(y.error.prec) <- 1
# diag(y.error.prec[-1, ]) <- 1/4
# vertical_neighbours <- ((0:d^2) * (d^2 + 1) + d+1)
# y.error.prec[vertical_neighbours[(vertical_neighbours <= d^4)]] <- 1/4
# y.error.prec[upper.tri(y.error.prec)] = t(y.error.prec)[upper.tri(y.error.prec)]

nl_data <- nl_obs(d, sigmaX, nu, delta, y.error.prec, Time.step)
obs <- matrix(unlist(nl_data$y), ncol = d, nrow = d)

Nparticles <- 100
# initial state
history <- sqrt(sigmaX)*array(rnorm(Nparticles*d^2), dim = c(d, d, Nparticles))
for (t in 1:Time.step){
  print(paste(t))
  res <- dac_nl_lightweight(history, obs, sigmaX, nu)
  history <- res
}
apply(res, c(1,2), mean)
mean((apply(res, c(1,2), mean) - nl_data$x[, , Time.step+1])^2)
