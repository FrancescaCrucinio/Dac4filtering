# Lightweight resampling for linear Gaussian SSM
marginal_lgssm_light <- function(i, u, nv, ci, W, Nparticles, theta, lambda, tau, x, history, memory = FALSE){
  nchild <- 2
  # resample on each children
  if(u == 1){
    # child 1
    indices1 <- stratified_resample(W[, nchild*(i-1)+1], theta*Nparticles)
    # child 2 (with random permutation)
    indices2 <- sample(stratified_resample(W[, nchild*i], theta*Nparticles))
  }
  else{
    # child 1
    indices1 <- sample.int(Nparticles, size = theta*Nparticles, replace = TRUE)
    # child 2
    indices2 <- sample.int(Nparticles, size = theta*Nparticles, replace = TRUE)
  }
  # mixture weights
  if(memory){
    lWmix <- rep(0, times = theta*Nparticles)
    for (n in 1:(theta*Nparticles)) {
      lWmix[n] <- -0.5*lambda * (lambda *x[indices1[n], (ci[1]+nv-1)]^2/(tau+lambda) -
                                   2*x[indices1[n], (ci[1]+nv-1)] * x[indices2[n], (ci[1]+nv)])+
        log(mean(exp(-0.5*(x[indices1[n], (ci[1]+nv-1)]-history[, ci[1]+nv])^2/(tau+lambda) -
                       0.5*lambda*tau*history[, ci[1]+nv] * x[indices1[n], (ci[1]+nv-1)]/(tau+lambda)))) -
        log(mean(exp(-0.5*(x[indices1[n], (ci[1]+nv-1)]-history[, ci[1]+nv])^2/(tau+lambda))))
    }
  } else {
    lWmix <- -0.5*lambda * (lambda *x[indices1, (ci[1]+nv-1)]^2/(tau+lambda) -
                              2*x[indices1, (ci[1]+nv-1)] * x[indices2, (ci[1]+nv)]) +
      log(rowMeans(exp(-0.5*outer(x[indices1, (ci[1]+nv-1)], history[, ci[1]+nv], "-")^2/(tau+lambda) -
                         0.5*lambda*tau*outer(x[indices1, (ci[1]+nv-1)], history[, ci[1]+nv], "*")/(tau+lambda)))) -
      log(rowMeans(exp(-0.5*outer(x[indices1, (ci[1]+nv-1)], history[, ci[1]+nv], "-")^2/(tau+lambda))))
  }
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(indices1[indices], indices2[indices]), "resampled_particles_lW" = lWmix[indices]))
}


# Adaptive Lightweight resampling for linear Gaussian SSM
marginal_lgssm_light_adaptive <- function(ess_target, i, u, nv, ci, lW, Nparticles, lambda, tau, x, history){
  # binary tree
  nchild <- 2
  # normal_density <- function(z, coeff_mat, prec_mat, past) {
  #   normal_mean <- -sweep(0.5*coeff_mat %*% past, 1, z)
  #   out <- exp(-0.5*apply(normal_mean, 2, FUN = function(z) {t(z)%*%prec_mat%*%z}))
  #   return(out)
  # }

  # mixture weights
  # integral_left <- apply(x[, ci[1]:(ci[1]+nv-1), drop = FALSE], 1, FUN = normal_density,
  #                    x.coeff[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1)], x.error.prec[ci[1]:(ci[1]+nv-1), ci[1]:(ci[1]+nv-1)],
  #                    t(history[, ci[1]:(ci[1]+nv-1)]))
  # integral_right <- apply(x[, (ci[1]+nv):ci[2], drop = FALSE], 1, FUN = normal_density,
  #                       x.coeff[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2]], x.error.prec[(ci[1]+nv):ci[2],(ci[1]+nv):ci[2]],
  #                       t(history[, (ci[1]+nv):ci[2]]))
  # integral_merged <- apply(x[, ci[1]:ci[2], drop = FALSE], 1, FUN = normal_density,
  #                      x.coeff[ci[1]:ci[2], ci[1]:ci[2]], x.error.prec[ci[1]:ci[2],ci[1]:ci[2]],
  #                      t(history[, ci[1]:ci[2]]))
  integral_part <- log(rowMeans(exp(-0.5*outer(x[, (ci[1]+nv-1)], history[, ci[1]+nv], "-")^2/(tau+lambda) -
                                      0.5*lambda*tau*outer(x[, (ci[1]+nv-1)], history[, ci[1]+nv], "*")/(tau+lambda))))  -
            log(rowMeans(exp(-0.5*outer(x[, (ci[1]+nv-1)], history[, ci[1]+nv], "-")^2/(tau+lambda))))
  if(u == 1){
    lWmix <- lW[, (nchild*(i-1)+1)] + lW[, i*nchild] - 0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) -
                                                                      2*x[, (ci[1]+nv-1)] * x[, (ci[1]+nv)])
  } else{
    lWmix <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) -
                              2*x[, (ci[1]+nv-1)] * x[, (ci[1]+nv)])
  }
  # if(u == 1){
  #   lWmix <- lW[, (nchild*(i-1)+1)] + lW[, i*nchild] +
  #     log(colSums(integral_merged)) - log(sum(integral_left)) - log(sum(integral_right))
  # } else{
  #   lWmix <- log(sum(integral_merged)) - log(sum(integral_left)) - log(sum(integral_right))
  # }
  lWmix <- lWmix + integral_part
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # build ESS
  ess_s <- sum(Wmix)
  ess_ss <- sum(Wmix^2)
  ess <- ess_s^2/ess_ss
  # first permutation
  permutation <- 1:Nparticles
  theta <- 1
  while (ess < ess_target & theta <= ceiling(sqrt(Nparticles))) {
    theta <- theta+1
    new_perm <- sample.int(Nparticles)
    # integral_right <- apply(x[new_perm, (ci[1]+nv):ci[2], drop = FALSE], 1, FUN = normal_density,
    #                         x.coeff[(ci[1]+nv):ci[2], (ci[1]+nv):ci[2]], x.error.prec[(ci[1]+nv):ci[2],(ci[1]+nv):ci[2]],
    #                         t(history[, (ci[1]+nv):ci[2]]))
    # integral_merged <- apply(cbind(x[, ci[1]:(ci[1]+nv-1)], x[new_perm, (ci[1]+nv):ci[2]]), 1, FUN = normal_density,
    #                          x.coeff[ci[1]:ci[2], ci[1]:ci[2]], x.error.prec[ci[1]:ci[2],ci[1]:ci[2]],
    #                          t(history[, ci[1]:ci[2]]))
    if(u == 1){
      lWmix_perm <- lW[, (nchild*(i-1)+1)] + lW[new_perm, i*nchild] -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) -
                                                                                     2*x[, (ci[1]+nv-1)] * x[new_perm, (ci[1]+nv)])
    } else{
      lWmix_perm <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2/(tau+lambda) -
                                     2*x[, (ci[1]+nv-1)] * x[new_perm, (ci[1]+nv)])
    }
    # if(u == 1){
    #   lWmix_perm <- lW[, (nchild*(i-1)+1)] + lW[new_perm, i*nchild] +
    #     log(sum(integral_merged)) - log(sum(integral_left)) - log(sum(integral_right))
    # } else{
    #   lWmix_perm <- log(sum(integral_merged)) - log(sum(integral_left)) - log(sum(integral_right))
    # }
    lWmix_perm <- lWmix_perm + integral_part
    permutation <- c(permutation, new_perm)
    max.lWmix <- max(lWmix_perm)
    Wmix <- exp(lWmix_perm - max.lWmix)
    # build ESS
    ess_s <- ess_s + sum(Wmix)
    ess_ss <- ess_ss + sum(Wmix^2)
    ess <- ess_s^2/ess_ss
    lWmix <- c(lWmix, lWmix_perm)
  }
  # write.table(data.frame("u" = u, "theta" = theta), file = "data/adaptive_lgssm.csv", sep = ",", append = TRUE, quote = FALSE,
  #             col.names = FALSE, row.names = FALSE)
  max.lWmix <- max(lWmix)
  Wmix <- exp(lWmix - max.lWmix)
  # resampling the new population
  indices <- stratified_resample(Wmix/sum(Wmix), Nparticles)
  return(list("resampled_indices" = cbind(rep(1:Nparticles, times = theta)[indices], permutation[indices]),
              "target_reached" = (ess >= ess_target), "resampled_particles_lW" = lWmix[indices]))
}
