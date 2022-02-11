lgssm_tempering <- function(ci, nv, lambda, tau, x, xOldv, lW, ess_target, ess_decay_threshold, mcmc_sd){
  Nparticles <- lenght(lW)
  current_alpha <- bisection_ess(lW, ess_target)
  lOmega <- rep(0, times = Nparticles)

  while (abs(current_alpha - 1) > 1e-03) {
    new_alpha <- bisection_cess(lOmega, lW, current_alpha, ess_decay_threshold)$root
    newlOmega <- lOmega + (new_alpha - current_alpha)*lW
    Omega <- exp(newlOmega - max(newlOmega))
    ess <- sum(Omega)^2/sum(Omega^2)
    if(ess < Nparticles / 2){
      resampled_indices <- stratified_resample(Omega/sum(Omega), Nparticles)
      newlOmega <- rep(0, times = Nparticles)
      xNew[, ci[1]:ci[2]] <- cbind(x[resampled_indices[, 1], ci[1]:(ci[1]+nv-1)], x[resampled_indices[, 2], (ci[1]+nv):ci[2]])
      # update lW
      lW <- -0.5*lambda * (lambda *x[, (ci[1]+nv-1)]^2 - 2*x[, (ci[1]+nv-1)] * (x[, (ci[1]+nv)] - 0.5*tau*xOldv)
      )
    }
    # MCMC move
    xNew[, ci[1]:ci[2]] <- xNew[, ci[1]:ci[2]] + mcmc_sd*rnorm(dim(x[, ci[1]:ci[2]]))
    current_alpha <- new_alpha
    lOmega <- newlOmega
    x <- xNew
  }
  return(list("x" = x, "resampled_indices" = resampled_indices))
}
