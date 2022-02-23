# nl_obs <- function(d, sigmaX, nu, delta, rx, ry, Time.Step){
#   y.error.prec <- matrix(0, nrow = d^2, ncol = d^2)
#   diag(y.error.prec) <- 1
#   if(ry > 3) ry <- 3
#
#   x <- array(0, dims = c(d, d, Time.step+1))
#   y <- array(0, dims = c(d, d, Time.step))
#   # initial state
#   x[, , 1] <- sqrt(sigmaX)*matrix(rnorm(d^2), nrow = 2)
#
#   for (1 in 1:Time.Step) {
#     for (row in 1:d) {
#       for(col in 1:d){
#
#
#       }
#     }
#
#   }
#   return(x, y)
# }
#
#
# find_x_neighbours <- function(rx, d){
#   if(rx > 3) rx <- 3
# }
# # dimensions of y
# p <- nrow(y.error.var)
#
#
# y <- matrix(0, nrow = Time.step, ncol = p)
# # initial state
# x <- mvrnorm(n = 1, mu0, Sigma0)
# # loop over time
# for(i in 1:Time.step){
#   x <- t(rmvnp(1, c(x.coeff%*%x), x.error.prec))
#   y[i, ] <- mvrnorm(n = 1, y.coeff%*%x, y.error.var)
# }
# return(y)
