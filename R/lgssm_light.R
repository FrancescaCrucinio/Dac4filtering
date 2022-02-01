lgssm_light <- function(i, u, nv, ci, W, Nparticles, m, lambda, tau, x, xOld){
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
  lWmix <- -0.5*lambda * (lambda *x[indices1, (ci[1]+nv-1)]^2/(tau+lambda) -
                            2*x[indices1, (ci[1]+nv-1)] * (x[indices2, (ci[1]+nv)] - 0.5*tau*xOld[indices2]/(tau+lambda))
  )
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(cbind(indices1[indices], indices2[indices]))
}
# xOldv is the ci[1]+nv component of xOld (main code)
lgssm_adaptive_light <- function(ess_target, i, u, nv, ci, lW, Nparticles, lambda, tau, x, xOldv){
  # binary tree
  nchild <- 2
  # mixture weights
  if(u == 1){
    lWmix <- lW[, (nchild*(i-1)+1)] + lW[, i*nchild] -
      0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOldv/(tau+lambda))
      )
  } else{
    lWmix <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOldv/(tau+lambda))
    )
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # build ESS
  ess_s <- sum(Wmix)
  ess_ss <- sum(Wmix^2)
  ess <- ess_s^2/ess_ss
  # first permutation
  permutation <- 1:Nparticles
  m <- 1
  while (ess < ess_target & m <= ceiling(sqrt(Nparticles))) {
    m <- m+1
    new_perm <- sample.int(Nparticles)
    if(u == 1){
      lWmix_perm <- lW[, (nchild*(i-1)+1)] + lW[new_perm, i*nchild] -
        0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[new_perm, (ci[1]+nv)] - 0.5*tau*xOldv[new_perm])
        )
    } else{
      lWmix_perm <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[new_perm, (ci[1]+nv)] - 0.5*tau*xOldv[new_perm])
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
  write.table(data.frame("u" = u, "m" = m), file = "data/adaptive_lgssm.csv", sep = ",", append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(cbind(rep(1:Nparticles, time = m)[indices], permutation[indices]))
}
