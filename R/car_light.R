car_light <- function(i, u, nv, nvNew, ci, W, Nparticles, m, sigmaX, x, xOld){
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
  lWmix <- rep(0, times = m*Nparticles)
  for (n in 1:(m*Nparticles)) {
    # merge the two children nodes
    mx <- c(x[indices1[n], ci[1]:(ci[1]+nv-1)], x[indices2[n], (ci[1]+nv):ci[2]])
    # get last term in mixture weights
    tmp <- 0
    for (h in 1:nv){
      tmp <- tmp + (nv-h+1)*xOld[indices2[n], d-h+1]
    }
    tmp <- sum(x[indices1[n], ci[1]:(ci[1]+nv-1)])*tmp/(d^2*sigmaX)
    lWmix[n] <- sum(x[indices1[n], ci[1]:(ci[1]+nv-1)])*sum(x[indices2[n], (ci[1]+nv):ci[2]])/(d*sigmaX) -
      (sum(cumsum(mx[1:(nvNew-1)]/d)^2) - sum(cumsum(x[indices1[n], ci[1]:(ci[1]+nv-1)][seq(length.out = (nv-1))]/d)^2) -
         sum(cumsum(x[indices2[n], (ci[1]+nv):ci[2]][seq(length.out = (nv-1))]/d)^2))/(2*sigmaX) - tmp
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(cbind(indices1[indices], indices2[indices]))
}
car_adaptive_light <- function(ess_target, i, u, nv, ci, lW, Nparticles, lambda, tau, x, xOld){
  # binary tree
  nchild <- 2
  # mixture weights
  lWmix <- rep(0, times = Nparticles)
  for (n in 1:Nparticles) {
    # merge the two children nodes
    mx <- c(x[indices1[n], ci[1]:(ci[1]+nv-1)], x[indices2[n], (ci[1]+nv):ci[2]])
    # get last term in mixture weights
    tmp <- 0
    for (h in 1:nv){
      tmp <- tmp + (nv-h+1)*xOld[indices2[n], d-h+1]
    }
    tmp <- sum(x[indices1[n], ci[1]:(ci[1]+nv-1)])*tmp/(d^2*sigmaX)
    lWmix[n] <- sum(x[indices1[n], ci[1]:(ci[1]+nv-1)])*sum(x[indices2[n], (ci[1]+nv):ci[2]])/(d*sigmaX) -
      (sum(cumsum(mx[1:(nvNew-1)]/d)^2) - sum(cumsum(x[indices1[n], ci[1]:(ci[1]+nv-1)][seq(length.out = (nv-1))]/d)^2) -
         sum(cumsum(x[indices2[n], (ci[1]+nv):ci[2]][seq(length.out = (nv-1))]/d)^2))/(2*sigmaX) - tmp
  }
  if(u == 1){
    lWmix <- lWmix + lW[, (nchild*(i-1)+1)] + lW[, i*nchild]
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
    # mixture weights
    lWmix_perm <- rep(0, times = Nparticles)
    for (n in 1:Nparticles) {
      # merge the two children nodes
      mx <- c(x[indices1[n], ci[1]:(ci[1]+nv-1)], x[indices2[n], (ci[1]+nv):ci[2]])
      # get last term in mixture weights
      tmp <- 0
      for (h in 1:nv){
        tmp <- tmp + (nv-h+1)*xOld[new_perm[n], d-h+1]
      }
      tmp <- sum(x[indices1[n], ci[1]:(ci[1]+nv-1)])*tmp/(d^2*sigmaX)
      lWmix[n] <- sum(x[indices1[n], ci[1]:(ci[1]+nv-1)])*sum(x[new_perm[n], (ci[1]+nv):ci[2]])/(d*sigmaX) -
        (sum(cumsum(mx[1:(nvNew-1)]/d)^2) - sum(cumsum(x[indices1[n], ci[1]:(ci[1]+nv-1)][seq(length.out = (nv-1))]/d)^2) -
           sum(cumsum(x[indices2[n], (ci[1]+nv):ci[2]][seq(length.out = (nv-1))]/d)^2))/(2*sigmaX) - tmp
    }
    if(u == 1){
      lWmix_perm <- lWmix_perm + lW[, (nchild*(i-1)+1)] + lW[new_perm, i*nchild]
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
  write.table(data.frame("u" = u, "m" = m), file = "adaptive_lgssm.csv", sep = ",", append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(cbind(rep(1:Nparticles, time = m)[indices], permutation[indices]))
}

