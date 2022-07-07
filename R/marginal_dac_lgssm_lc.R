# Linear cost DaC for linear Gaussian SSM
marginal_dac_lgssm_lc <- function(history, obs, tau, lambda, sigmaY){
  # dimension and number of particles
  d <- ncol(history)
  Nparticles <- nrow(history)
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # number of variables
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  W <- matrix(0, nrow = Nparticles, ncol = d)
  for (i in 1:nchild^nlevels){
    if (i == 1) {
      x[, i] <- 0.5*history[, i] + rnorm(Nparticles)/sqrt(tau)
    } else {
      x[, i] <- 0.5*tau*history[, i]/(tau+lambda) + rnorm(Nparticles)/sqrt(tau+lambda)
    }
    lW <- -0.5*(obs[i] - x[, i])^2/sigmaY - 0.5*log(2*pi*sigmaY)
    max.lW <- max(lW)
    W[, i] <- exp(lW - max.lW)
    W[, i] <- W[, i]/sum(W[, i])
  }

  # loop over tree levels excluding leaves
  for (u in 1:nlevels){
    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    # number of variables in each node
    nvNew <- nchild^u

    # updated particles/normalizing constant
    xNew <- matrix(0, nrow = Nparticles, ncol = d)
    WNew <- matrix(0, nrow = Nparticles, ncol = nodes)

    for (i in 1:nodes){
      # get children indices
      ci <- child_indices(i, nvNew)
      # resample on each children
      # child 1
      indices1 <- stratified_resample(W[, nchild*(i-1)+1], Nparticles)
      # child 2 (with random permutation)
      indices2 <- sample(stratified_resample(W[, nchild*i], Nparticles))
      # weights
      integral_merged <- rep(0, times = Nparticles)
      integral_left <- rep(0, times = Nparticles)
      integral_right <- rep(0, times = Nparticles)
      for (n in 1:Nparticles) {
        mx <- x[indices1[n], ]
        mx[(ci[1]+nv):ci[2]] <- x[indices2[n], (ci[1]+nv):ci[2]]
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
        integral_merged[n] <- mean(exp(rowSums(integral_per_dimension)))
        integral_left[n] <- mean(exp(rowSums(integral_per_dimension[, 1:nv, drop = FALSE])))
        integral_right[n] <- mean(exp(rowSums(integral_per_dimension[, (nv+1):(2*nv), drop = FALSE])))
      }
      lWmix <- - 0.5*lambda * (lambda *x[indices1, (ci[1]+nv-1)]^2/(tau+lambda) -
                                        2*x[indices1, (ci[1]+nv-1)] * x[indices2, (ci[1]+nv)]) +
                  integral_merged - integral_left - integral_right
      max.lWmix <- max(lWmix)
      WNew[, i] <- exp(lWmix - max.lWmix)
      # update particles
      xNew[, ci[1]:ci[2]] <- cbind(x[indices1, ci[1]:(ci[1]+nv-1)], x[indices2, (ci[1]+nv):ci[2]])
    }
    x <- xNew
    nv <- nvNew
    W <- WNew
  }
  return(x)
}
