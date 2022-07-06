
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
            lWmix[n1, n2] <- lW[n1, (nchild*(i-1)+1)] + lW[n2, i*nchild] -
              -0.5*lambda * (lambda *x[n1, (ci[1]+nv-1)]^2/(tau+lambda) -
                               2*x[n1, (ci[1]+nv-1)] * x[n2, (ci[1]+nv)]) +
              log(mean(exp(-0.5*(x[n1, (ci[1]+nv-1)]-history[, ci[1]+nv])^2/(tau+lambda) -
                             0.5*lambda*tau*history[, ci[1]+nv] * x[n1, (ci[1]+nv-1)]/(tau+lambda)))) -
              log(mean(exp(-0.5*(x[n1, (ci[1]+nv-1)]-history[, ci[1]+nv])^2/(tau+lambda))))
          }
        }
        max.lWmix <- max(lWmix)
        Wmix <- exp(lWmix - max.lWmix)
      } else {
        for (n1 in 1:Nparticles) {
          for (n2 in 1:Nparticles) {
            lWmix[n1, n2] <- -0.5*lambda * (lambda *x[n1, (ci[1]+nv-1)]^2/(tau+lambda) -
                                              2*x[n1, (ci[1]+nv-1)] * x[n2, (ci[1]+nv)]) +
              log(mean(exp(-0.5*(x[n1, (ci[1]+nv-1)]-history[, ci[1]+nv])^2/(tau+lambda) -
                             0.5*lambda*tau*history[, ci[1]+nv] * x[n1, (ci[1]+nv-1)]/(tau+lambda)))) -
              log(mean(exp(-0.5*(x[n1, (ci[1]+nv-1)]-history[, ci[1]+nv])^2/(tau+lambda))))
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
