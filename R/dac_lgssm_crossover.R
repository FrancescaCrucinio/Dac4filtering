dac_lgssm_lc_crossover <- function(history, obs, tau, lambda, sigmaY){
  # dimension and number of particles
  d <- ncol(history[, , 1])
  Nparticles <- nrow(history[, , 1])
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # leaves
  # number of variables
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  W <- matrix(0, nrow = Nparticles, ncol = d)
  for (i in 1:nchild^nlevels){
    if (i == 1) {
      x[, i] <- 0.5*history[, i, 1] + rnorm(Nparticles)/sqrt(tau)
    } else {
      x[, i] <- 0.5*tau*history[, i, 1]/(tau+lambda) + rnorm(Nparticles)/sqrt(tau+lambda)
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
    xOldNew <- matrix(0, nrow = Nparticles, ncol = d)
    WNew <- matrix(0, nrow = Nparticles, ncol = nodes)

    for (i in 1:nodes){
      # get children indices
      ci <- child_indices(i, nvNew)
      # resample on each children
      # child 1
      indices1 <- stratified_resample(W[, nchild*(i-1)+1], Nparticles)
      # child 2 (with random permutation)
      indices2 <- sample(stratified_resample(W[, nchild*i], Nparticles))
      # update history
      history[, ci[1]:(ci[1]+nv-1)] <- history[indices1, ci[1]:(ci[1]+nv-1)]
      history[, (ci[1]+nv):ci[2]] <- history[indices2, (ci[1]+nv):ci[2]]
      # update particles
      x[, ci[1]:ci[2]] <- cbind(x[indices1, ci[1]:(ci[1]+nv-1)], x[indices2, (ci[1]+nv):ci[2]])
      # mutation
      # sample crossover point
      crossover_point <- sample.int(d)
      if(crossover_point < ci[2] & crossover_point >= ci[1]){
        # merge the x on the two children
        history[, crossover_point:ci[2], ] <-
      }
      # weights
      lW <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 -
                             2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*history[, (ci[1]+nv), 1])
      )
      max.lW <- max(lW)
      W[, i] <- exp(lW - max.lW)
    }
    nv <- nvNew
  }
  return(x)
}
