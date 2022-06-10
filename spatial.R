set.seed(1234*11)
d <- 4
sigmaX <- 1
nu <- 10
tau_diag <- 1
tau <- -0.25
Time.step <- 10
spatial_data <- spatial_obs(d, sigmaX, nu, tau, tau_diag, Time.step)
N <- 10000
history <- sqrt(sigmaX)*array(rnorm(N*d^2), dim = c(d, d, N))
tic()
for (t in 1:Time.step){
  res_dac <- marginal_dac_spatial(history, spatial_data$y[, , t], sigmaX, nu, tau, tau_diag, adaptive = TRUE)
  history <- res_dac
  print(paste(t))
}
toc()
df <- data.frame(mean = c(apply(res_dac, c(1,2), mean)), variance = c(apply(res_dac^2, c(1,2), mean) - apply(res_dac, c(1,2), mean)^2),
                 first_q = c(apply(res_dac, c(1,2), quantile, probs = 0.25)), third_q = c(apply(res_dac, c(1,2), quantile, probs = 0.75)))
write_csv(x=df, file = "d4N10000results")
df <- read.csv("d4N10000results")
mean((spatial_data$x[, , Time.step+1]- df$mean)^2)
Nparticles <- c(100, 500, 1000)
Nrep <- 10
mse_gt <- matrix(0, nrow = 3, ncol = Nrep)
mse_filter <- matrix(0, nrow = 3, ncol = Nrep)
var_filter <- matrix(0, nrow = 3, ncol = Nrep)
first_q_filter <- matrix(0, nrow = 3, ncol = Nrep)
third_q_filter <- matrix(0, nrow = 3, ncol = Nrep)
for (k in 1:Nrep) {
  for (i in 1:length(Nparticles)) {
    N <- Nparticles[i]
    # initial state
    history <- sqrt(sigmaX)*array(rnorm(N*d^2), dim = c(d, d, N))
    tic()
    for (t in 1:Time.step){
      res_dac <- marginal_dac_spatial(history, spatial_data$y[, , t], sigmaX, nu, tau, tau_diag, adaptive = TRUE)
      history <- res_dac
      print(paste(t))
    }
    toc()
    mse_gt[i, k] <- mean((c(apply(res_dac, c(1,2), mean)) - spatial_data$x[, , Time.step+1])^2)
    mse_filter[i, k] <- mean((c(apply(res_dac, c(1,2), mean)) - df$mean)^2)
    var_filter[i, k] <- mean((c(apply(res_dac^2, c(1,2), mean) - apply(res_dac, c(1,2), mean)^2) - df$var)^2)
    first_q_filter[i, k] <-mean((c(apply(res_dac, c(1,2), quantile, probs = 0.25)) - df$first_q)^2)
    third_q_filter[i, k] <- mean((c(apply(res_dac, c(1,2), quantile, probs = 0.75)) - df$third_q)^2)
  }
}
rowMeans(mse_gt)
rowMeans(mse_filter)
rowMeans(var_filter)
rowMeans(first_q_filter)
rowMeans(third_q_filter)

image((apply(res_dac, c(1,2), mean) - spatial_data$x[, , Time.step+1])^2)
(apply(res_dac, c(1,2), function(x) length(unique(x))))

