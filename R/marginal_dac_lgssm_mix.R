# Mixture resampling DaC for linear Gaussian SSM
marginal_dac_lgssm_mix <- function(history, obs, tau, lambda, sigmaY){
  # dimension and number of particles
  d <- ncol(history)
  Nparticles <- nrow(history)
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # number of variables
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  lW <- matrix(0, nrow = Nparticles, ncol = d)
  for (i in 1:nchild^nlevels){
    if (i == 1) {
      x[, i] <- 0.5*history[, i] + rnorm(Nparticles)/sqrt(tau)
    } else {
      x[, i] <- 0.5*tau*history[, i]/(tau+lambda) + rnorm(Nparticles)/sqrt(tau+lambda)
    }
    lW[, i] <- -0.5*(obs[i] - x[, i])^2/sigmaY - 0.5*log(2*pi*sigmaY)
    max.lW <- max(lW[, i])
    W <- exp(lW[, i] - max.lW)
  }

  # loop over tree levels excluding leaves
  for (u in 1:nlevels){
    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    # number of variables in each node
    nvNew <- nchild^u

    # updated particles/normalizing constant
    xNew <- matrix(0, nrow = Nparticles, ncol = d)
    for (i in 1:nodes){

      # mixture weights
      lWmix <- matrix(0, ncol = Nparticles, nrow = Nparticles)
      ci <- child_indices(i, nvNew)

      if(u == 1){
        for (n1 in 1:Nparticles) {
          for (n2 in 1:Nparticles) {
            mx <- x[n1, ]
            mx[(ci[1]+nv):ci[2]] <- x[n2, (ci[1]+nv):ci[2]]
            integral_per_dimension <- matrix(0, nrow = Nparticles, ncol = 2*nv)
            for (j in 2:(2*nv)) {
              dimension <- (ci[1]:ci[2])[j]
              integral_per_dimension[, j] <- -0.5*(tau+lambda)*(mx[dimension] - 0.5*tau*history[, dimension]/(tau+lambda))^2 -
                0.5*lambda*tau*(mx[dimension-1] * history[, dimension])/(tau+lambda)
            }
            if(ci[1] == 1){
              integral_per_dimension[, 1] <- -0.5*tau*(mx[1] - 0.5*history[, 1])^2
            } else {
              integral_per_dimension[, 1] <- -0.5*(tau+lambda)*(mx[ci[1]] - 0.5*tau*history[, ci[1]]/(tau+lambda))^2
            }
            integral_merged <- mean(exp(rowSums(integral_per_dimension)))
            integral_left <- mean(exp(rowSums(integral_per_dimension[, 1:nv, drop = FALSE])))
            integral_right <- mean(exp(rowSums(integral_per_dimension[, (nv+1):(2*nv), drop = FALSE])))
            lWmix[n1, n2] <- lW[n1, (nchild*(i-1)+1)] + lW[n2, i*nchild] -
                              0.5*lambda * (lambda *x[n1, (ci[1]+nv-1)]^2/(tau+lambda) -
                              2*x[n1, (ci[1]+nv-1)] * x[n2, (ci[1]+nv)]) +
                              log(integral_merged) - log(integral_left) - log(integral_right)
          }
        }
        max.lWmix <- max(lWmix)
        Wmix <- exp(lWmix - max.lWmix)
      } else {
        for (n1 in 1:Nparticles) {
          for (n2 in 1:Nparticles) {
            mx <- x[n1, ]
            mx[(ci[1]+nv):ci[2]] <- x[n2, (ci[1]+nv):ci[2]]
            integral_per_dimension <- matrix(0, nrow = Nparticles, ncol = 2*nv)
            for (j in 2:(2*nv)) {
              dimension <- (ci[1]:ci[2])[j]
              integral_per_dimension[, j] <- -0.5*(tau+lambda)*(mx[dimension] - 0.5*tau*history[, dimension]/(tau+lambda))^2 -
                0.5*lambda*tau*(mx[dimension-1] * history[, dimension])/(tau+lambda)
            }
            if(ci[1] == 1){
              integral_per_dimension[, 1] <- -0.5*tau*(mx[1] - 0.5*history[, 1])^2
            } else {
              integral_per_dimension[, 1] <- -0.5*(tau+lambda)*(mx[ci[1]] - 0.5*tau*history[, ci[1]]/(tau+lambda))^2
            }
            integral_merged <- mean(exp(rowSums(integral_per_dimension)))
            integral_left <- mean(exp(rowSums(integral_per_dimension[, 1:nv, drop = FALSE])))
            integral_right <- mean(exp(rowSums(integral_per_dimension[, (nv+1):(2*nv), drop = FALSE])))
            lWmix[n1, n2] <-  -0.5*lambda * (lambda *x[n1, (ci[1]+nv-1)]^2/(tau+lambda) -
                              2*x[n1, (ci[1]+nv-1)] * x[n2, (ci[1]+nv)]) +
                              log(integral_merged) - log(integral_left) - log(integral_right)
          }
        }
        max.lWmix <- max(lWmix)
        Wmix <- exp(lWmix - max.lWmix)
      }
      Wmix <- Wmix/sum(Wmix)
      # resampling
      indices <- stratified_resample(Wmix, Nparticles)
      # get row/column indices
      res <- rc_indices(indices, Nparticles)
      for(n in 1:Nparticles) {
        # update particles
        xNew[n, ci[1]:ci[2]] <- c(x[res[n, 1], ci[1]:(ci[1]+nv-1)], x[res[n, 2], (ci[1]+nv):ci[2]])
      }
    }
    x <- xNew
    nv <- nvNew
  }
  return(x)
}
