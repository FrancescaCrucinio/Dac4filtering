dac_car <- function(xOld, obs, sigmaX, sigmaY, Sigma.det){
  # dimension and number of particles
  d <- ncol(xOld)
  Nparticles <- nrow(xOld)
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # function returning index of first and last child of node i in level u

  # leaves
  # number of variables
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  lW <- matrix(0, nrow = Nparticles, ncol = d)
  lZ <- rep(0, times = d)
  for (i in 1:nchild^nlevels){
    # propose
    x[, i] <- rowSums(xOld[, i:d, drop = FALSE])/d + sqrt(sigmaX) * rnorm(Nparticles)
    # weights
    lW[, i] <- -0.5*(obs[i] - x[, i])^2/sigmaY - 0.5*log(2*pi*sigmaY)
    max.lW <- max(lW[, i])
    W <- exp(lW[, i] - max.lW)
    lZ[i] <- log(mean(W)) + max.lW
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
    lZNew <- rep(0, times = nodes)

    for (i in 1:nodes) {
      # mixture weights
      lWmix <- matrix(0, ncol = Nparticles, nrow = Nparticles)
      ci <- child_indices(i, nvNew)
      # all children are leaves
      if(u == 1){
        for (n1 in 1:Nparticles) {
          for (n2 in 1:Nparticles) {
            lWmix[n1, n2] <- lW[n1, (nchild*(i-1)+1)] + lW[n2, i*nchild] +
              x[n1, ci[1]:(ci[1]+nv-1)]*x[n2, (ci[1]+nv):ci[2]]/(d*sigmaX) -
              (x[n1, ci[1]:(ci[1]+nv-1)]/d)^2/(2*sigmaX) - 
              x[n1, ci[1]] * sum(xOld[n2, ci[1]:d])/(d^2*sigmaX)
          }
        }
        max.lWmix <- max(lWmix)
        Wmix <- exp(lWmix - max.lWmix)
        lZNew[i] <- log(mean(Wmix)) + max.lWmix -
          0.5*Sigma.det[[u]][nchild*(i-1)+1] - 0.5*Sigma.det[[u]][i*nchild] + 0.5*Sigma.det[[u+1]][i]
      } else {
        for (n1 in 1:Nparticles) {
          for (n2 in 1:Nparticles) {
            # merge the two children nodes
            mx <- c(x[n1, ci[1]:(ci[1]+nv-1)], x[n2, (ci[1]+nv):ci[2]])
            # get last term in mixture weights
            tmp <- 0
            for (i1 in ci[1]:(ci[1]+nv-1)){
              tmp <- tmp + x[n1, i1] * sum(xOld[n2, i1:d])
            }
            tmp <- nv*tmp/(d^2*sigmaX)
            lWmix[n1, n2] <- sum(x[n1, ci[1]:(ci[1]+nv-1)])*sum(x[n2, (ci[1]+nv):ci[2]])/(d*sigmaX) -
              (sum(cumsum(mx[1:(nvNew-1)]/d)^2) - sum(cumsum(x[n1, ci[1]:(ci[1]+nv-1)][1:(nv-1)]/d)^2) -
                 sum(cumsum(x[n2, (ci[1]+nv):ci[2]][1:(nv-1)]/d)^2))/(2*sigmaX) - tmp
          }
        }
        max.lWmix <- max(lWmix)
        Wmix <- exp(lWmix - max.lWmix)
        lZNew[i] <- log(mean(Wmix)) + max.lWmix + lZ[(nchild*(i-1)+1)] + lZ[i*nchild] -
          0.5*Sigma.det[[u]][nchild*(i-1)+1] - 0.5*Sigma.det[[u]][i*nchild] + 0.5*Sigma.det[[u+1]][i]
      }
      Wmix <- Wmix/sum(Wmix)
      # resampling
      res <- mult_resample(Wmix, Nparticles)
      for(n in 1:Nparticles) {
        # update particles
        xNew[n, ci[1]:ci[2]] <- c(x[res[n, 1], ci[1]:(ci[1]+nv-1)], x[res[n, 2], (ci[1]+nv):ci[2]])
        # update xOld
        xOldNew[n, ci[1]:ci[2]] <- c(xOld[res[n, 1], ci[1]:(ci[1]+nv-1)], xOld[res[n, 2], (ci[1]+nv):ci[2]])
      }
    }
    x <- xNew
    xOld <- xOldNew
    lZ <- lZNew
    nv <- nvNew
  }
  return(cbind(x, lZ))
}
