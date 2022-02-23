car_adaptive_light <- function(ess_target, i, u, nv, nvNew, ci, lW, Nparticles, sigmaX, x, xOld, historyIndex)
{
  # binary tree
  nchild <- 2
  # mixture weights
  lWmix <- rep(0, times = Nparticles)
  for (n in 1:Nparticles) {
    # get ancestor
    right_ancestor_coordinates <- cbind(historyIndex[n, (ci[1]+nv):d, nchild*i], (ci[1]+nv):d)
    right_ancestor <- xOld[right_ancestor_coordinates]
    # merge the two children nodes
    mx <- x[n, ci[1]:ci[2]]
    lWmix[n] <- sum(x[n, ci[1]:(ci[1]+nv-1)])*sum(x[n, (ci[1]+nv):ci[2]])/(d*sigmaX) -
      (sum(cumsum(mx[1:(nvNew-1)]/d)^2) - sum(cumsum(x[n, ci[1]:(ci[1]+nv-1)][seq(length.out = (nv-1))]/d)^2) -
         sum(cumsum(x[n, (ci[1]+nv):ci[2]][seq(length.out = (nv-1))]/d)^2))/(2*sigmaX) -
      sum(x[n, ci[1]:(ci[1]+nv-1)])*sum(cumsum(rev(right_ancestor)))/(d^2*sigmaX)
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
  while (ess < ess_target & m<=ceiling(sqrt(Nparticles))) {
    m <- m+1
    new_perm <- sample.int(Nparticles)
    # mixture weights
    lWmix_perm <- rep(0, times = Nparticles)
    for (n in 1:Nparticles) {
      # merge the two children nodes
      mx <- c(x[n, ci[1]:(ci[1]+nv-1)], x[new_perm[n], (ci[1]+nv):ci[2]])
      lWmix_perm[n] <- sum(x[n, ci[1]:(ci[1]+nv-1)])*sum(x[new_perm[n], (ci[1]+nv):ci[2]])/(d*sigmaX) -
        (sum(cumsum(mx[1:(nvNew-1)]/d)^2) - sum(cumsum(x[n, ci[1]:(ci[1]+nv-1)][seq(length.out = (nv-1))]/d)^2) -
           sum(cumsum(x[new_perm[n], (ci[1]+nv):ci[2]][seq(length.out = (nv-1))]/d)^2))/(2*sigmaX) -
        sum(x[n, ci[1]:(ci[1]+nv-1)])*sum(cumsum(rev(xOld[new_perm[n], (ci[1]+nv):d])))/(d^2*sigmaX)
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
  # write.table(data.frame("u" = u, "m" = m), file = "data/adaptive_car.csv", sep = ",", append = TRUE, quote = FALSE,
  #             col.names = FALSE, row.names = FALSE)
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(rep(1:Nparticles, times = m)[indices], permutation[indices]),
              "target_reached" = (ess >= ess_target), "resampled_particles_lW" = lWmix[indices]))
}
