


x1 <- seq(-10, 10, length=100)
x2 <- x1

set.seed(1234*11)
d <- 32
sigmaX <- 1
nu <- 10
tau_diag <- 1
tau <- -0.25
Time.step <- 1
spatial_data <- spatial_obs(d, sigmaX, nu, tau, tau_diag, Time.step)
y1 <- spatial_data$y[1,1,]
y2 <- spatial_data$y[1,2,]

z <- matrix(0, ncol = 100, nrow = 100)
for (i in 1:100) {
  for (j in 1:100) {
    sum_over_neighbours_obs_left <- (x1[i]-y1)^2
    sum_over_neighbours_obs_right <- (x2[j]-y2)^2
    sum_over_neighbours_obs_merged <- (x1[i]-y1)^2+(x2[j]-y2)^2-0.5*(y1-x1[i])*(y2-x2[j])
    z[i,j] <- - 0.5*(nu+2)*log(1+abs(sum_over_neighbours_obs_merged)/nu)
    + 0.5*(nu+1)*(log(1+abs(sum_over_neighbours_obs_left)/nu) + log(1+abs(sum_over_neighbours_obs_right)/nu))
  }
}
image(z)
axis(1,x1) # add a new x-axis
axis(2,x2, pos=0)
