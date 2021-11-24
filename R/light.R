light <- function(i, u, nv, ci, W, Nparticles, m, lambda, tau, x, xOld, indicesOld){
  # binary tree
  nchild <- 2
  # resample on each children
  if(u == 1){
    # child 1
    indices1 <- stratified_resample(W[, nchild*(i-1)+1], m*Nparticles)
    # child 2 (with random permutation)
    indices2 <- sample(stratified_resample(W[, nchild*i], m*Nparticles))
  }
  else{
    # child 1
    indices1 <- sample.int(Nparticles, size = m*Nparticles, replace = TRUE)
    # child 2
    indices2 <- sample.int(Nparticles, size = m*Nparticles, replace = TRUE)
  }
  # mixture weights
  lWmix <- -0.5*lambda * (lambda *x[indices1, (ci[1]+nv-1)]^2 -
                            2*x[indices1, (ci[1]+nv-1)] * (x[indices2, (ci[1]+nv)] - 0.5*tau*xOld[indicesOld[indices2, (ci[1]+nv)], (ci[1]+nv)])
  )
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # lZNew[i] <- lZ[(nchild*(i-1)+1)] + lZ[i*nchild] + log(mean(Wmix)) + max.lWmix -
  # 0.5*Sigma.det[[u]][nchild*(i-1)+1] - 0.5*Sigma.det[[u]][i*nchild] + 0.5*Sigma.det[[u+1]][i]
  # resampling the new population
  indices <- mult_resample(Wmix/sum(Wmix), Nparticles)
  return(cbind(indices1[indices], indices2[indices]))
}
