msv_bpf <- function(history, obs_past, obs_current, SigmaV, SigmaUV, SigmaX, phi){
  d <- ncol(history)
  Nparticles <- nrow(history)
  x <- matrix(0, nrow = Nparticles, ncol = d)
  lW <- rep(0, times = Nparticles)
  for (i in 1:Nparticles) {
    mu <- phi*history[i, ] + SigmaUV %*% solve(SigmaV) %*% diag(1/sqrt(exp(history[i, ]))) %*% obs_past
    x[i, ] <- mu + mvrnorm(n = 1, rep(0, d), SigmaX)
    obs_covariance <- diag(sqrt(exp(x[i, ]))) %*% SigmaV %*% diag(sqrt(exp(x[i, ])))
    lW[i] <- -0.5*(obs_current %*% solve(obs_covariance) %*% obs_current)
  }
  max.lW <- max(lW)
  W <- exp(lW - max(lW))
  W <- W/sum(W)
  ancestors <- stratified_resample(W, Nparticles)
  x <- x[ancestors, ]
  return(x)
}

