adaptive_light <- function(ess_target, i, u, nv, ci, lW, Nparticles, lambda, tau, x, xOld){
  # binary tree
  nchild <- 2
  # mixture weights
  if(u == 1){
    lWmix <- lW[, (nchild*(i-1)+1)] + lW[, i*nchild] -
      0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOld[, (ci[1]+nv)])
    )
  } else{
    lWmix <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOld[, (ci[1]+nv)])
      )
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # build ESS
  ess_s <- sum(Wmix)
  ess_ss <- sum(Wmix^2)
  ess <- ess_s^2/ess_ss
  ess <- 100
  # first permutation
  permutation <- 1:Nparticles
  m <- 1
  while (ess < ess_target) {
    m <- m+1
    new_perm <- sample.int(Nparticles)
    if(u == 1){
      lWmix_perm <- lW[, (nchild*(i-1)+1)] + lW[new_perm, i*nchild] -
        0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[new_perm, (ci[1]+nv)] - 0.5*tau*xOld[new_perm, (ci[1]+nv)])
        )
    } else{
      lWmix_perm <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[new_perm, (ci[1]+nv)] - 0.5*tau*xOld[new_perm, (ci[1]+nv)])
      )
    }
    permutation <- c(permutation, new_perm)
    max.lWmix <- max(lWmix_perm)
    Wmix <- exp(lWmix_perm - max.lWmix)
    # build ESS
    ess_s <- ess_s + sum(Wmix)
    ess_ss <- ess_ss + sum(Wmix^2)
    ess <- ess_s^2/ess_ss
    lWmix <- c(lWmix, lWmix_perm)
  }
  print(paste(m, "ESS", ess))
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(cbind(rep(1:Nparticles, time = m)[indices], permutation[indices]))
}
