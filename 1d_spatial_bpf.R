set.seed(1234*11)
sigmaX <- 1
nu <- 10
tau_diag <- 1
Time.step <- 100
spatial_data <- spatial_obs1D(sigmaX, nu, tau_diag, Time.step)

# BPF
Nparticles <- c(10, 50, 100, 500, 1000, 5000, 10000, 20000, 30000)
mse <- matrix(0, ncol = length(Nparticles), nrow = Time.step+1)
var <- matrix(0, ncol = length(Nparticles), nrow = Time.step+1)
first_q <- matrix(0, ncol = length(Nparticles), nrow = Time.step+1)
third_q <- matrix(0, ncol = length(Nparticles), nrow = Time.step+1)
for (i in 1:length(Nparticles)) {
  N <- Nparticles[i]
  x <- matrix(0, ncol = N, nrow = Time.step+1)
  x[1, ] <- rnorm(N, sd = sqrt(sigmaX))
  for (t in 2:(Time.step+1)){
    obs <- spatial_data$y[t-1]
    x[t, ] <- x[t-1, ] + sqrt(sigmaX)*rnorm(N)
    lW <- -0.5*(nu+1)*log(1+(x[t, ] - obs)^2/nu)
    max.lW <- max(lW)
    W <- exp(lW - max(lW))
    W <- W/sum(W)
    ancestors <- stratified_resample(W, N)
    x[t, ] <- x[t, ancestors]
  }
  mse[, i] <- (rowMeans(x) - spatial_data$x)^2
  var[, i] <- rowMeans(x^2) - rowMeans(x)^2
  first_q[, i] <- apply(x, 1, quantile, probs = 0.25)
  third_q[, i] <- apply(x, 1, quantile, probs = 0.75)
  print(paste(i))
}
plot(1:(Time.step+1), (var[, 1] - var[, length(Nparticles)])^2, type = "l", log='y')
lines(1:(Time.step+1), (var[, 2] - var[, length(Nparticles)])^2, type = "l", col = "blue")
lines(1:(Time.step+1), (var[, 3] - var[, length(Nparticles)])^2, type = "l", col = "red")
lines(1:(Time.step+1), (var[, 4] - var[, length(Nparticles)])^2, type = "l", col = "green")
lines(1:(Time.step+1), (var[, 5] - var[, length(Nparticles)])^2, type = "l", col = "yellow")
lines(1:(Time.step+1), (var[, 6] - var[, length(Nparticles)])^2, type = "l", col = "orange")
lines(1:(Time.step+1), (var[, 7] - var[, length(Nparticles)])^2, type = "l", col = "brown")
colMeans(mse)
colMeans((var[, 1:(length(Nparticles)-1)] - var[, length(Nparticles)])^2)
colMeans((first_q[, 1:(length(Nparticles)-1)] - first_q[, length(Nparticles)])^2)
colMeans((third_q[, 1:(length(Nparticles)-1)] - third_q[, length(Nparticles)])^2)

