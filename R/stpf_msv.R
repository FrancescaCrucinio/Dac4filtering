# Space-time particle filter for MSV
# Space-time particle filter for linear Gaussian SSM
stpf_msv_first_step <- function(Nparticles, M, obs_current, Sigma0, SigmaV){
  # dimension, number islands and number of particles
  d <- dim(Sigma0)[1]

  x <- array(0, dim = c(Nparticles, M, d))
  lZ <- rep(0, times = Nparticles)
  # loop over islands
  for(i in 1:Nparticles){
    # loop over dimension
    j <- 1
    # get mean vector
    x[i, , j] <- sqrt(Sigma0[j, j])*rnorm(M)
    # weights
    lW <- -0.5*(obs_current[j])^2/(exp(x[i, , j])*SigmaV[j,j])
    max.lW <- max(lW)
    W <- exp(lW - max(lW))
    lZ[i] <- log(mean(W)) + max.lW
    # resampling
    W <- W/sum(W)
    ancestors <- stratified_resample(W, M)
    x[i, , j] <- x[i, ancestors, j]
    j <- 2
    # propose
    mu_conditional <- sapply(x[i, , 1:(j-1)], function(s){Sigma0[j, 1:(j-1)] %*% solve(Sigma0[1:(j-1), 1:(j-1)]) %*% as.matrix(s)})
    x[i, , j] <- mu_conditional + sqrt(Sigma0[j, j] - Sigma0[j, 1:(j-1)] %*% solve(Sigma0[1:(j-1), 1:(j-1)]) %*% Sigma0[j, 1:(j-1)])[1]*rnorm(M)
    # weights
    lW <- -0.5*apply(x[i, , 1:j], 1, function(s) {obs_current[1:j] %*% solve(diag(sqrt(exp(s))) %*% SigmaV[1:j, 1:j] %*% diag(sqrt(exp(s)))) %*% obs_current[1:j]}) +
      0.5*sapply(x[i, , 1:(j-1)], function(s) {obs_current[1:(j-1)] %*% solve(diag(as.matrix(sqrt(exp(s)))) %*% SigmaV[1:(j-1), 1:(j-1)] %*% diag(as.matrix(sqrt(exp(s))))) %*% obs_current[1:(j-1)]})
    max.lW <- max(lW)
    W <- exp(lW - max(lW))
    lZ[i] <- lZ[i] + log(mean(W)) + max.lW
    # resampling
    W <- W/sum(W)
    ancestors <- stratified_resample(W, M)
    x[i, , 1:j] <- x[i, ancestors, 1:j]
    for(j in 3:d){
      # propose
      mu_conditional <- apply(x[i, , 1:(j-1)], 1, function(s){Sigma0[j, 1:(j-1)] %*% solve(Sigma0[1:(j-1), 1:(j-1)]) %*% s})
      x[i, , j] <- mu_conditional + sqrt(Sigma0[j, j] - Sigma0[j, 1:(j-1)] %*% solve(Sigma0[1:(j-1), 1:(j-1)]) %*% Sigma0[j, 1:(j-1)])[1]*rnorm(M)
      # weights
      lW <- -0.5*apply(x[i, , 1:j], 1, function(s) {obs_current[1:j] %*% solve(diag(sqrt(exp(s))) %*% SigmaV[1:j, 1:j] %*% diag(sqrt(exp(s)))) %*% obs_current[1:j]}) +
        0.5*apply(x[i, , 1:(j-1)], 1, function(s) {obs_current[1:(j-1)] %*% solve(diag(sqrt(exp(s))) %*% SigmaV[1:(j-1), 1:(j-1)] %*% diag(sqrt(exp(s)))) %*% obs_current[1:(j-1)]})
      max.lW <- max(lW)
      W <- exp(lW - max(lW))
      lZ[i] <- lZ[i] + log(mean(W)) + max.lW
      # resampling
      W <- W/sum(W)
      ancestors <- stratified_resample(W, M)
      x[i, , 1:j] <- x[i, ancestors, 1:j]
    }
  }
  # resampling islands
  Wisland <- exp(lZ - max(lZ))
  Wisland <- Wisland/sum(Wisland)
  ancestors_island <- stratified_resample(Wisland, Nparticles)
  x <- x[ancestors_island, ,]
  return(x)
}

stpf_msv <- function(xOld, obs_current, obs_past, phi, SigmaUV, SigmaV, SigmaX){
  # dimension, number islands and number of particles
  d <- dim(xOld)[3]
  Nparticles <- nrow(xOld)
  M <- ncol(xOld)

  x <- array(0, dim = c(Nparticles, M, d))
  lZ <- rep(0, times = Nparticles)
  # loop over islands
  for(i in 1:Nparticles){
    # loop over dimension
    j <- 1
    # get mean vector
    mu <- matrix(0, nrow = M, ncol = d)
    for (k in 1:M) {
      mu[k, ] <- phi*xOld[i, k, ] + SigmaUV %*% solve(SigmaV) %*% diag(1/sqrt(exp(xOld[i, k, ]))) %*% obs_past
    }
    x[i, , j] <- mu[, j] + sqrt(SigmaX[j,j])*rnorm(M)
    # weights
    lW <- -0.5*(obs_current[j])^2/(exp(x[i, , j])*SigmaV[j,j])
    max.lW <- max(lW)
    W <- exp(lW - max(lW))
    lZ[i] <- log(mean(W)) + max.lW
    # resampling
    W <- W/sum(W)
    ancestors <- stratified_resample(W, M)
    xOld[i, , ] <- xOld[i, ancestors, ]
    x[i, , j] <- x[i, ancestors, j]
    j <- 2
    # propose
    mu_conditional <- rep(0, times = M)
    for (k in 1:M) {
      mu_conditional[k] <- mu[k, j] + SigmaX[j, 1:(j-1)] %*% solve(SigmaX[1:(j-1), 1:(j-1)]) %*% (x[i, k, 1:(j-1)] - mu[k, 1:(j-1)])
    }
    x[i, , j] <- mu_conditional + sqrt(SigmaX[j, j] - SigmaX[j, 1:(j-1)] %*% solve(SigmaX[1:(j-1), 1:(j-1)]) %*% SigmaX[j, 1:(j-1)])[1]*rnorm(M)
    # weights
    lW <- -0.5*apply(x[i, , 1:j], 1, function(s) {obs_current[1:j] %*% solve(diag(sqrt(exp(s))) %*% SigmaV[1:j, 1:j] %*% diag(sqrt(exp(s)))) %*% obs_current[1:j]}) +
      0.5*sapply(x[i, , 1:(j-1)], function(s) {obs_current[1:(j-1)] %*% solve(diag(as.matrix(sqrt(exp(s)))) %*% SigmaV[1:(j-1), 1:(j-1)] %*% diag(as.matrix(sqrt(exp(s))))) %*% obs_current[1:(j-1)]})
    max.lW <- max(lW)
    W <- exp(lW - max(lW))
    lZ[i] <- lZ[i] + log(mean(W)) + max.lW
    # resampling
    W <- W/sum(W)
    ancestors <- stratified_resample(W, M)
    xOld[i, , ] <- xOld[i, ancestors, ]
    x[i, , 1:j] <- x[i, ancestors, 1:j]
    for(j in 3:d){
      # propose
      mu_conditional <- rep(0, times = M)
      for (k in 1:M) {
        mu_conditional[k] <- mu[k, j] + SigmaX[j, 1:(j-1)] %*% solve(SigmaX[1:(j-1), 1:(j-1)]) %*% (x[i, k, 1:(j-1)] - mu[k, 1:(j-1)])
      }
      x[i, , j] <- mu_conditional + sqrt(SigmaX[j, j] - SigmaX[j, 1:(j-1)] %*% solve(SigmaX[1:(j-1), 1:(j-1)]) %*% SigmaX[j, 1:(j-1)])[1]*rnorm(M)
      # weights
      lW <- -0.5*apply(x[i, , 1:j], 1, function(s) {obs_current[1:j] %*% solve(diag(sqrt(exp(s))) %*% SigmaV[1:j, 1:j] %*% diag(sqrt(exp(s)))) %*% obs_current[1:j]}) +
        0.5*apply(x[i, , 1:(j-1)], 1, function(s) {obs_current[1:(j-1)] %*% solve(diag(sqrt(exp(s))) %*% SigmaV[1:(j-1), 1:(j-1)] %*% diag(sqrt(exp(s)))) %*% obs_current[1:(j-1)]})
      max.lW <- max(lW)
      W <- exp(lW - max(lW))
      lZ[i] <- lZ[i] + log(mean(W)) + max.lW
      # resampling
      W <- W/sum(W)
      ancestors <- stratified_resample(W, M)
      xOld[i, , ] <- xOld[i, ancestors, ]
      x[i, , 1:j] <- x[i, ancestors, 1:j]
    }
  }
  # resampling islands
  Wisland <- exp(lZ - max(lZ))
  Wisland <- Wisland/sum(Wisland)
  ancestors_island <- stratified_resample(Wisland, Nparticles)
  x <- x[ancestors_island, ,]
  return(x)
}
