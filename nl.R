d <- 2
sigmaX <- 1
nu <- 10
delta <- 1
Time.step <- 1

nl_data <- nl_obs(d, sigmaX, nu, delta, Time.step)
obs <- matrix(unlist(nl_data$y), ncol = d, nrow = d)

Nparticles <- 1000
# initial state
history <- array(0, dim = c(Nparticles, d, d, 2))
history[, , , 1] <- sqrt(sigmaX)*array(rnorm(Nparticles*d^2), dim = c(Nparticles, d, d))
res <- dac_nl_lightweight(history, obs, sigmaX, nu)
colMeans(res)
